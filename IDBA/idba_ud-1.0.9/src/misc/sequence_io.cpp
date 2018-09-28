/**
 * @file sequence_io.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-06
 */

#include "misc/sequence_io.h"

#include <string>

#include "misc/utils.h"
#include "sequence/sequence.h"
#include "misc/sequence_reader.h"
#include "misc/sequence_writer.h"
#include "sequence/short_sequence.h"

using namespace std;


uint64_t ReadSequence(const string &filename, deque<Sequence> &sequences)
{
    FastaReader reader(filename);
    Sequence seq;
    while (reader.Read(seq))
        sequences.push_back(seq);
    return sequences.size();
}

uint64_t ReadSequence(const string &filename, deque<Sequence> &sequences, deque<string> &names)
{
    FastaReader reader(filename);
    Sequence seq;
    string name;
    while (reader.Read(seq, name))
    {
        sequences.push_back(seq);
        names.push_back(name);
    }
    return sequences.size();
}

uint64_t ReadSequence(const string &filename, deque<ShortSequence> &sequences)
{
    FastaReader reader(filename);
    Sequence seq;
    while (reader.Read(seq))
    {
        seq.TrimN();
        ShortSequence short_seq(seq);
        sequences.push_back(short_seq);
    }
    return sequences.size();
}

uint64_t WriteSequence(const string &filename, const deque<Sequence> &sequences, const string &prefix)
{
    FastaWriter writer(filename);
    for (uint64_t i = 0; i < sequences.size(); ++i)
        writer.Write(sequences[i], FormatString("%s_%d", prefix.c_str(), i));
    return sequences.size();
}

uint64_t WriteSequence(FastaWriter &writer, const deque<Sequence> &sequences, const string &prefix)
{
    for (uint64_t i = 0; i < sequences.size(); ++i)
        writer.Write(sequences[i], FormatString("%s_%d", prefix.c_str(), i));
    return sequences.size();
}

