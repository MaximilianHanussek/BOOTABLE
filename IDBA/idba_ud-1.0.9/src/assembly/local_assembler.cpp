/**
 * @file local_assembler.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.9
 * @date 2012-01-05
 */

#include "assembly/local_assembler.h"

#include <pthread.h>

#include <vector>
#include <deque>

#include "assembly/assembly_utility.h"
#include "basic/histgram.h"
#include "graph/contig_graph.h"
#include "graph/hash_graph.h"
#include "misc/hash_aligner.h"
#include "sequence/sequence.h"
#include "sequence/short_sequence.h"

using namespace std;

LocalAssembler::~LocalAssembler()
{
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)locks_.size(); ++i)
        omp_destroy_lock(&locks_[i]);
}

void LocalAssembler::Initialize(AssemblyInfo &assembly_info, deque<Sequence> &contigs)
{
    assembly_info_ = &assembly_info;
    contigs_ = &contigs;

    locks_.resize(contigs.size());
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)locks_.size(); ++i)
        omp_init_lock(&locks_[i]);
    in_reads_.resize(contigs.size());
    out_reads_.resize(contigs.size());

    read_length_ = assembly_info.read_length();
}

int64_t LocalAssembler::Assemble(deque<Sequence> &contigs)
{
    if (num_threads_ == 0)
        num_threads_ = omp_get_max_threads();

    contigs.clear();
    omp_set_num_threads(1);

    vector<pthread_t> threads(num_threads_);
    vector<LocalAssemblyTask> tasks(num_threads_);
    for (int i = 0; i < num_threads_; ++i)
    {
        tasks[i].id = i;
        tasks[i].local_assembler = this;
        pthread_create(&threads[i], NULL, LocalAssembleThread, (void *)&tasks[i]);
    }
    
    for (int i = 0; i < num_threads_; ++i)
    {
        pthread_join(threads[i], NULL);
        contigs.insert(contigs.end(), tasks[i].local_contigs.begin(), tasks[i].local_contigs.end());
    }

    omp_set_num_threads(num_threads_);

    return contigs.size();
}

void * LocalAssembler::LocalAssembleThread(void *p)
{
    omp_set_num_threads(1);

    LocalAssemblyTask &task = *(LocalAssemblyTask *)p;
    LocalAssembler &local_assembler = *task.local_assembler;
    deque<Sequence> &contigs = local_assembler.contigs();
    vector<deque<uint64_t> > &in_reads = local_assembler.in_reads();
    vector<deque<uint64_t> > &out_reads = local_assembler.out_reads();
    int num_threads = local_assembler.num_threads();
    int min_contig = local_assembler.min_contig();

    for (int64_t i = task.id; i < (int64_t)contigs.size(); i += num_threads)
    {
        if ((int)contigs[i].size() < min_contig)
            continue;

        {
            deque<Sequence> local_contigs;
            local_assembler.LocalAssemble(contigs[i], in_reads[i], local_contigs);
            contigs[i].ReverseComplement();
            task.local_contigs.insert(task.local_contigs.end(), local_contigs.begin(), local_contigs.end());
        }

        {
            deque<Sequence> local_contigs;
            local_assembler.LocalAssemble(contigs[i], out_reads[i], local_contigs);
            contigs[i].ReverseComplement();
            task.local_contigs.insert(task.local_contigs.end(), local_contigs.begin(), local_contigs.end());
        }
    }

    return p;
}

void LocalAssembler::LocalAssemble(Sequence &contig, std::deque<uint64_t> &local_reads, std::deque<Sequence> &local_contigs)
{
    deque<ShortSequence> &reads = assembly_info_->reads;

    int min_num_reads = median_ / read_length_;

    if ((int)local_reads.size() > min_num_reads)
    {
        sort(local_reads.begin(), local_reads.end());

        AssemblyInfo assembly_info;
        int last = -1;
        int count = 0;
        for (unsigned j = 0; j < local_reads.size(); ++j)
        {
            int offset = local_reads[j] >> 36;
            int index = local_reads[j] & ((1ULL << 32) - 1);

            count = (offset == last) ? count+1 : 1;
            last = offset;

            if (count <= 3)
                assembly_info.reads.push_back(reads[index]);
        }
        Sequence seq = contig;
        if ((int64_t)seq.size() > local_range_)
            seq.resize(local_range_);
        assembly_info.long_reads.push_back(seq);
        assembly_info.ClearStatus();

        IterativeAssemble(assembly_info, local_contigs);
    }
}

void LocalAssembler::IterativeAssemble(AssemblyInfo &assembly_info, deque<Sequence> &contigs)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    Sequence contig_tail = assembly_info.long_reads[0];
    assembly_info.long_reads.clear();
    assembly_info.ClearStatus();

    HashGraph hash_graph(mink_);

    hash_graph.reserve(local_range_*4);
    deque<ContigInfo> contig_infos;
    ContigGraph contig_graph;
    for (int kmer_size = mink_; kmer_size <= maxk_; kmer_size += step_)
    {
        int64_t sum = 0;
        hash_graph.set_kmer_size(kmer_size);
        hash_graph.clear();
        for (int64_t i = 0; i < (int64_t)reads.size(); ++i)
        {
            if ((int)reads[i].size() < kmer_size)
                continue;

            Sequence seq(reads[i]);
            hash_graph.InsertKmers(seq);
            sum += seq.size() - kmer_size + 1;
        }

        Histgram<int> histgram = hash_graph.coverage_histgram();
        double mean = histgram.percentile(1 - 1.0 * local_range_ / hash_graph.num_vertices());
        double threshold = mean;

        hash_graph.InsertKmers(contig_tail);
        for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
            hash_graph.InsertUncountKmers(contigs[i]);

        hash_graph.Assemble(contigs, contig_infos);
        contig_graph.set_kmer_size(kmer_size);
        contig_graph.Initialize(contigs, contig_infos);
        contig_graph.RemoveDeadEnd(kmer_size*2);
        contig_graph.RemoveBubble();
        contig_graph.IterateCoverage(kmer_size*2, 1, threshold);
        contig_graph.Assemble(contigs, contig_infos);

        if (contigs.size() == 1)
            break;
    }

    int index = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() > min_contig_)
            contigs[index++].swap(contigs[i]);
    }
    contigs.resize(index);
}

