/**
 * @file local_assmebler.h
 * @brief Local Assembler which assembles paired-end reads locally.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.9
 * @date 2012-01-05
 */

#ifndef __ASSEMBLY_LOCAL_ASSEMBLER_H_

#define __ASSEMBLY_LOCAL_ASSEMBLER_H_

#include "assembly/assembly_utility.h"
#include "sequence/sequence.h"
#include "graph/hash_graph.h"
#include "misc/hash_aligner.h"

#include <omp.h>

#include <deque>
#include <vector>

class LocalAssembler;

struct LocalAssemblyTask
{
    int id;
    LocalAssembler *local_assembler;
    std::deque<Sequence> local_contigs;
};

class LocalAssembler
{
public:
    LocalAssembler() {}
    ~LocalAssembler();

    void Initialize(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs);

    AssemblyInfo &assembly_info() { return *assembly_info_; }
    std::deque<Sequence> &contigs() { return *contigs_; }
    std::vector<std::deque<uint64_t> > &in_reads() { return in_reads_; }
    std::vector<std::deque<uint64_t> > &out_reads() { return out_reads_; }

    int num_threads() const { return num_threads_; }
    void set_num_threads(int num_threads) { num_threads_ = num_threads; }

    int mink() const { return mink_; }
    void set_mink(int mink) { mink_ = mink; }

    int maxk() const { return maxk_; }
    void set_maxk(int maxk) { maxk_ = maxk; }
    
    int step() const { return step_; }
    void set_step(int step) { step_ = step; }

    int min_contig() const { return min_contig_; }
    void set_min_contig(int min_contig) { min_contig_ = min_contig; }

    void set_insert_distance(double median, double sd)
    { median_ = median; sd_ = sd; local_range_ = std::min(median*2, median + 3*sd); }
    double median() const { return median_; }
    double sd() const { return sd_; }
    double local_range() const { return local_range_; }

    int64_t Assemble(std::deque<Sequence> &contigs);

    void AddReadByHashAlignerRecord(HashAlignerRecord &record, int read_id)
    {
        if (record.ref_from >= record.ref_length - local_range_)
        {
            int offset = ((record.ref_length - record.ref_from) << 4) | (record.query_length - record.match_length);
            if (record.is_reverse == false)
                AddOutRead(record.ref_id, offset, read_id ^ 1);
            else
                AddInRead(record.ref_id, offset, read_id ^ 1);
        }
    }

    void AddInRead(int contig_id, int offset, int read_id)
    { 
        omp_set_lock(&locks_[contig_id]); 
        in_reads_[contig_id].push_back(((uint64_t)offset << 32) | read_id); 
        omp_unset_lock(&locks_[contig_id]); 
    }

    void AddOutRead(int contig_id, int offset, int read_id)
    { 
        omp_set_lock(&locks_[contig_id]); 
        out_reads_[contig_id].push_back(((uint64_t)offset << 32) | read_id); 
        omp_unset_lock(&locks_[contig_id]); 
    }

    static void *LocalAssembleThread(void *p);

    void LocalAssemble(Sequence &contig, std::deque<uint64_t> &local_reads, std::deque<Sequence> &local_contigs);
    void IterativeAssemble(AssemblyInfo &assembly_info, std::deque<Sequence> &contigs);

private:
    AssemblyInfo *assembly_info_;
    std::deque<Sequence> *contigs_;
    std::vector<omp_lock_t> locks_;
    std::vector<std::deque<uint64_t> > in_reads_;
    std::vector<std::deque<uint64_t> > out_reads_;
    int read_length_;
    int num_threads_;
    int mink_;
    int maxk_;
    int step_;
    int min_contig_;
    double median_;
    double sd_;
    double local_range_;
};

#endif

