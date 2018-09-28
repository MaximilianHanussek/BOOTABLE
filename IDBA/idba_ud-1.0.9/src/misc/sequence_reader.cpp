/**
 * @file sequence_reader.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-02
 */

#include "misc/sequence_reader.h"

#include <fstream>
#include <istream>
#include <stdexcept>
#include <string>

#include "sequence/sequence.h"


using namespace std;

SequenceReader::SequenceReader(std::istream &is, bool is_own)
{
    if (is.fail())
        throw logic_error("SequenceReader::SequenceReader() istream is invalid");

    is_ = &is;
    is_own_ = is_own;
}

bool SequenceReader::Read(Sequence &seq)
{
    string comment;
    string quality;
    return ReadRecord(seq, comment, quality);
}

bool SequenceReader::Read(Sequence &seq, string &comment)
{
    string quality;
    return ReadRecord(seq, comment, quality);
}

bool SequenceReader::Read(Sequence &seq, string &comment, string &quality)
{
    return ReadRecord(seq, comment, quality);
}

bool FastaReader::ReadRecord(Sequence &seq, string &comment, string &quality)
{
    return ReadFasta(*is_, seq, comment);
}

bool FastqReader::ReadRecord(Sequence &seq, string &comment, string &quality)
{
    return ReadFastq(*is_, seq, comment, quality);
}

