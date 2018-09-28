/**
 * @file assembly_utility.h
 * @brief Utility functions for assembly.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-06
 */

#ifndef __ASSEMBLY_ASSEMBLY_UTILITY_H_

#define __ASSEMBLY_ASSEMBLY_UTILITY_H_


#include <stdint.h>

#include <string>
#include <deque>
#include <vector>
#include <algorithm>

#include "graph/contig_info.h"
#include "misc/sequence_io.h"

class ContigGraph;
class HashGraph;
class Sequence;
class ShortSequence;
class Sequence;
class HashAlignerRecord;
class HashAligner;

struct AssemblyInfo
{
    std::deque<ShortSequence> reads;
    std::deque<Sequence> long_reads;
    std::vector<bool> read_flags;
    std::vector<bool> long_read_flags;

    void ClearStatus()
    {
        read_flags.resize(reads.size());
        long_read_flags.resize(long_reads.size());
        std::fill(read_flags.begin(), read_flags.end(), true);
        std::fill(long_read_flags.begin(), long_read_flags.end(), true);
    }

    int read_length() const;
};

void ReadInput(const std::string &read_file, const std::string &long_read_file, AssemblyInfo &assembly_info);
bool ReadHashAlignerRecords(FILE *fp, std::deque<HashAlignerRecord> &records);
void ReadHashAlignerRecordBlock(FILE *fp, std::vector<HashAlignerRecord> &records);
void WriteHashAlignerRecordBlock(FILE *fp, std::vector<HashAlignerRecord> &records);

uint64_t WriteContig(const std::string &filename, const std::deque<Sequence> &contigs, 
        const std::deque<ContigInfo> &contig_infos, const std::string &prefix = "", int min_contig = 0);
uint64_t WriteContig(FastaWriter &writer, const std::deque<Sequence> &contigs, 
        const std::deque<ContigInfo> &contig_infos, const std::string &prefix = "", int min_contig = 0);

void BuildKmerFile(AssemblyInfo &assembly_info, int kmer_size, int min_count, int prefix_length, const std::string &kmer_file);
void ReadKmerFile(const std::string &kmer_file, HashGraph &hash_graph);
void WriteKmerFile(const std::string &kmer_file, HashGraph &hash_graph);
void InsertInternalKmers(AssemblyInfo &assembly_info, HashGraph &hash_graph, int min_count = 0);
void IterateHashGraph(AssemblyInfo &assembly_info, int old_kmer_size, int new_kmer_size, int min_support, HashGraph &hash_graph);
void IterateContigGraph(AssemblyInfo &assembly_info, int old_kmer_size, int new_kmer_size, int min_support, ContigGraph &contig_graph);

void IterateHashGraph(AssemblyInfo &assembly_info, int new_kmer_size, int min_support, HashGraph &hash_graph, std::deque<Sequence> &old_contigs);

void InsertExistKmers(AssemblyInfo &assembly_info, HashGraph &hash_graph);
int64_t InsertIterativeKmers(const HashGraph &old_hash_graph, const Sequence &seq, HashGraph &hash_graph);

int64_t AlignReads(AssemblyInfo &assembly_info, HashAligner &hash_aligner, double similar, const std::string &align_file, bool is_all = true);
void CorrectReads(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs, 
        std::deque<ContigInfo> &contig_infos, const std::string &align_file, int max_mismatch);

#endif

