/**
 * @file hash_aligner.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-12
 */

#include "misc/hash_aligner.h"
#include "basic/bit_operation.h"

#include <omp.h>

#include <deque>
#include <iostream>

#include "sequence/sequence.h"

using namespace std;

bool Compare(const HashAlignerRecord &x, const HashAlignerRecord &y)
{
    return x.match_length > y.match_length;
}

HashAligner::HashAligner(uint32_t kmer_size, uint32_t min_length, uint32_t step)
{
    kmer_size_ = kmer_size;
    min_length_ = min_length;
    step_ = step;

    //if (num_threads == 0)
    int num_threads = omp_get_max_threads();

    buffer_records.resize(num_threads);
    buffer_tables.resize(num_threads);
    for (int i = 0; i < num_threads; ++i)
    {
        buffer_records[i].reserve(10000);
        buffer_tables[i].reserve(10000);
    }
}

void HashAligner::Initialize(const deque<Sequence> &sequences)
{
    uint64_t sum = 0;
    for (unsigned i = 0; i < sequences.size(); ++i)
    {
        if (sequences[i].size() >= min_length_)
            sum += sequences[i].size();
    }
    hash_map_.reserve(sum/step_);

    sequences_.resize(sequences.size());
    sequence_words.resize(sequences.size());
    reverse_sequence_words.resize(sequences.size());
//#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)sequences.size(); ++i)
    {
        if (sequences[i].size() < min_length_)
            continue;

        InsertSequence(sequences[i], i);

        sequences_[i] = sequences[i];
        Convert(sequences_[i], sequence_words[i]);
        sequences_[i].ReverseComplement();
        Convert(sequences_[i], reverse_sequence_words[i]);
        sequences_[i].ReverseComplement();
    }
}

int HashAligner::AlignSequence(const Sequence &seq, deque<HashAlignerRecord> &records, int min_match, int max_records)
{
    int max_mismatch = seq.size() - min_match;

    vector<uint64_t> words;
    Convert(seq, words);

    vector<HashAlignerRecord> &tmp_records = buffer_records[omp_get_thread_num()];
    tmp_records.resize(0);

    HashSet<uint64_t> &table = buffer_tables[omp_get_thread_num()];
    table.clear();

    records.clear();

    Kmer kmer(kmer_size_);
    int length = 0;
    for (uint64_t i = 0; i < seq.size(); ++i)
    {
        kmer.ShiftAppend(seq[i]);
        length = (seq[4] < 4) ? length + 1 : 0;

        if (length < (int)kmer_size_)
            continue;

        Kmer key = kmer.unique_format();
        HashMap<Kmer, position_list_type>::iterator p = hash_map_.find(key);

        if (p == hash_map_.end())
            continue;

        position_list_type &position_list = p->second;
        for (position_iterator iter = position_list.begin(); iter != position_list.end(); ++iter)
        {
            uint64_t position = *iter;
            uint64_t id = position >> 33;
            bool is_reverse = (position >> 32) & 1;
            uint64_t offset = position & ((1ULL << 32) - 1);

            if (key != kmer)
            {
                is_reverse = !is_reverse;
                offset = sequences_[id].size() - (offset + kmer_size_);
            }

            int64_t from = int64_t(offset) - (i+1 - kmer_size_);
            int64_t to = from + seq.size();

            if (from >= 0 && to <= sequences_[id].size())
            {
                uint64_t value = (id << 33) | (uint64_t(is_reverse) << 32) | from;
                if (table.insert(value).second)
                {
                    HashAlignerRecord record;
                    record.query_from = 0;
                    record.query_to = seq.size();
                    record.query_length = seq.size();
                    record.ref_id = id;
                    record.ref_from = from;
                    record.ref_to = to;
                    record.ref_length = sequences_[id].size();
                    record.is_reverse = is_reverse;
                    record.match_length = 0;

                    Match(words, record, max_mismatch);
                    if (record.match_length >= min_match)
                    {
                        if ((int)records.size() >= max_records)
                        {
                            records.resize(0);
                            return 0;
                        }
                        else
                            records.push_back(record);
                    }
                }
            }
        }
    }

    return records.size();
} 

void HashAligner::Match(const vector<uint64_t> &words, HashAlignerRecord &record, int max_mismatch)
{
    vector<uint64_t> &ref_words = (!record.is_reverse ? 
            sequence_words[record.ref_id] : reverse_sequence_words[record.ref_id]);

    int mismatch = 0;
    int length = record.query_to - record.query_from;
    for (int offset = 0; offset < length; offset += 32)
    {
        int len = min(32, length - offset);

        uint64_t x = GetWord(words, record.query_from + offset, record.query_from + offset + len);
        uint64_t y = GetWord(ref_words, record.ref_from + offset, record.ref_from + offset + len);

        mismatch += bit_operation::BaseCount(x ^ y);

        if (mismatch > max_mismatch)
        {
            mismatch = length;
            break;
        }
    }

    record.match_length = length - mismatch;
}

void HashAligner::InsertSequence(const Sequence &seq, uint64_t id)
{
    Kmer kmer(kmer_size_);
    int length = 0;
    int last = -1;
    for (int64_t i = 0; i < (int64_t)seq.size(); ++i)
    {
        kmer.ShiftAppend(seq[i]);
        length = (seq[4] < 4) ? length + 1 : 0;

        if (length < (int)kmer_size_)
            continue;

        if (i - last < (int)step_)
            continue;
        last = i;

        Kmer key = kmer.unique_format();

        uint64_t position = ((id << 1) << 32) | (i+1 - kmer_size_);
        if (key != kmer)
            position = (((id << 1) + 1) << 32) | (seq.size() - (i+1));

        hash_map_[key].set_pool(pool_);
        hash_map_[key].push_front(position);
    }
}

void HashAligner::ExtendRecord(const Sequence &seq, HashAlignerRecord &record)
{
    Sequence a = seq;
    const Sequence &b = sequences_[record.ref_id];

    for (int strand = 0; strand < 2; ++strand)
    {
        int mismatch = 0;
        
        while (record.query_to < record.query_length && record.ref_to < record.ref_length)
        {
            int x = record.query_to;
            int y = record.ref_to;

            if (!record.is_reverse)
            {
                if (a[x] != b[y])
                    ++mismatch;
                else
                    ++record.match_length;
            }
            else
            {
                if (a[x] != 3 - b[b.size() - 1 - y])
                    ++mismatch;
                else
                    ++record.match_length;
            }

            if (mismatch > 3)
                break;

            ++record.query_to;
            ++record.ref_to;
        }

        a.ReverseComplement();
        record.ReverseComplement();
    }
}

void HashAligner::Match(const Sequence &seq, HashAlignerRecord &record)
{
    const Sequence &a = seq;
    const Sequence &b = sequences_[record.ref_id];

    int match = 0;
    for (int i = 0; record.query_from + i < record.query_to; ++i)
    {
        int x = record.query_from + i;
        int y = record.ref_from + i;

        if (!record.is_reverse)
        {
            if (a[x] == b[y])
                ++match;
        }
        else
        {
            if (a[x] == 3 - b[b.size() - 1 - y])
                ++match;
        }
    }

    record.match_length = match;
}

void HashAligner::Convert(const Sequence &seq, vector<uint64_t> &words)
{
    words.resize((seq.size() + 31) >> 5);
    fill(words.begin(), words.end(), 0);
    for (unsigned i = 0; i < seq.size(); ++i)
        words[i>>5] |= uint64_t(seq[i]) << ((i&31) << 1);
}


