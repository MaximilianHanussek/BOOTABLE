/**
 * @file validate_contigs_blat.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.3
 * @date 2011-09-06
 */

#include <algorithm>
#include <cstdio>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "misc/blat_record.h"
#include "misc/options_description.h"
#include "misc/sequence_io.h"
#include "misc/utils.h"
#include "sequence/sequence.h"

using namespace std;

const int MaxLine = (1 << 20);

char line[MaxLine];
char buf[MaxLine];

int main(int argc, char *argv[])
{
    int min_contig = 100;
    double similar = 0.95;
    
    OptionsDescription desc;
    desc.AddOption("min_contig", "", min_contig, "minimum contigs");
    desc.AddOption("similar", "", similar, "similarity");

    try
    {
        desc.Parse(argc, argv);
    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "validate_contigs_blat - validate contigs by blat." << endl;
        cerr << "Usage: validate_contigs_blat ref.fa contigs.fa." << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }

    deque<Sequence> refs;
    deque<string> ref_names;
    ReadSequence(argv[1], refs, ref_names);

    deque<Sequence> contigs;
    deque<string> contig_names;
    ReadSequence(argv[2], contigs, contig_names);

    vector<vector<bool> > flags(refs.size());
    map<string, int> dict;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        flags[i].resize(refs[i].size(), false);
        size_t index = ref_names[i].find(' ');
        if (index != string::npos)
            ref_names[i].resize(index);
        dict[ref_names[i]] = i;
    }

    int num_gaps = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        size_t index = contig_names[i].find(' ');
        if (index != string::npos)
            contig_names[i].resize(index);

        bool is_new_gap = true;
        for (unsigned j = 0; j < contigs[i].size(); ++j)
        {
            if (contigs[i][j] == 4)
            {
                if (is_new_gap)
                {
                    is_new_gap = false;
                    ++num_gaps;
                }
            }
            else
                is_new_gap = true;
        }
    }

    string blat_file = string(argv[2]) + ".blat";
    FILE *fblat = OpenFile(blat_file, "rb");

    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);
    fgets(line, MaxLine, fblat);

    map<string, int> valid_contigs;
    deque<int> valid_lengths;
    int64_t num_mismatch = 0;
    while (fgets(line, MaxLine, fblat) != NULL)
    {
        BlatRecord record;
        record.Parse(line);

        int ref_id = dict[record.ref_name];

        if (record.match_count > similar * record.query_length && record.query_length >= min_contig
                && record.match_count > similar * abs(record.ref_to - record.ref_from))
        {
            for (unsigned i = 0; i < record.blocks.size(); ++i)
            {
                BlatBlock block = record.blocks[i];
                for (unsigned j = block.ref_from; j < block.ref_from + block.size; ++j)
                    flags[ref_id][j] = true;
            }

//            valid_lengths.push_back(abs(record.query_to - record.query_from));
//            valid_contigs.insert(record.query_name);
            //if (valid_contigs.insert(record.query_name).second)
            if (valid_contigs.find(record.query_name) == valid_contigs.end())
            {
                valid_contigs[record.query_name] = record.mismatch_count;
                valid_lengths.push_back(record.query_to - record.query_from);
                //valid_lengths.push_back(record.query_length);
                //num_mismatch += record.mismatch_count;
            }
            else
            {
                valid_contigs[record.query_name] = min(record.mismatch_count, (int64_t)valid_contigs[record.query_name]);
            }
        }
    }

    for (map<string, int>::iterator p = valid_contigs.begin(); p != valid_contigs.end(); ++p)
    {
        num_mismatch += p->second;
    }

    long long count = 0;
    long long total = 0;
    for (unsigned k = 0; k < flags.size(); ++k)
    {
        for (unsigned i = 0; i < flags[k].size(); ++i)
        {
            if (flags[k][i])
                ++count;
            ++total;
        }
    }

    //valid_lengths.push_back(60000);
    sort(valid_lengths.begin(), valid_lengths.end());
    reverse(valid_lengths.begin(), valid_lengths.end());

    long long n50 = 0;
    long long sum = 0;
    long long n80 = 0;

    for (unsigned i = 0; i < valid_lengths.size(); ++i)
    {
        sum += valid_lengths[i];
        if (sum >= 0.5 * total && n50 == 0)
            n50 = valid_lengths[i];
        if (sum >= 0.8 * total && n80 == 0)
            n80 = valid_lengths[i];
    }
    cout << "total " << total << " " << sum << endl;

    long long maximum = 0;
    long long mean = 0;
    if (valid_lengths.size() > 0)
    {
        maximum = valid_lengths[0];
        mean = sum / valid_lengths.size();
    }

    long long sum_wrong = 0;
    long long num_wrong = 0;
    long long corret_contigs = 0;
    long long sum_corret = 0;
    int last_id = 0;
    int last_error = 0;
    FastaWriter error_writer(FormatString("%s.error.fa", argv[2]));
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() < min_contig)
            continue;

        if (valid_contigs.find(contig_names[i]) == valid_contigs.end())
        {
            ++num_wrong;
            sum_wrong += contigs[i].size();
            error_writer.Write(contigs[i], contig_names[i]);
        }
        else
        {
            last_id = i;
            last_error = sum_wrong;

            ++corret_contigs;
            sum_corret += contigs[i].size();
        }
    }

    printf("last id %d %d total contigs %d gaps %d\n", last_id, last_error, (int)(num_wrong + corret_contigs), num_gaps);
    printf("contigs: %lld N50: %lld coverage: %.2f%% max: %lld mean: %lld total: %lld/%lld N80: %lld\n",
            (long long)valid_contigs.size(), n50, count * 100.0 / total, maximum, mean, count, total, n80);
    printf("substitution error: %.4f%% wrong contigs: %lld %lld correct: %lld %lld %s\n", 
            num_mismatch * 100.0 /sum, num_wrong, sum_wrong, corret_contigs, sum_corret, argv[2]);

    deque<int> lengths;
    for (unsigned i = 0; i < refs.size(); ++i)
    {
        int last = 0;
        for (unsigned j = 0; j < refs[i].size(); ++j)
        {
            if (flags[i][j] == 0)
            {
                if (flags[i][last])
                {
                    lengths.push_back(j - last);
                    last = j;
                }
            }
            else
            {
                if (flags[i][last] == 0)
                    last = j;
            }
        }
    }
    sort(lengths.begin(), lengths.end());
    reverse(lengths.begin(), lengths.end());
    
    return 0;
}

