/**
 * @file idba_ud.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.6
 * @date 2011-08-06
 */

#include <cmath>
#include <cstdio>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include <sys/stat.h>
#include <unistd.h>

#include "assembly/assembly_utility.h"
#include "basic/bit_operation.h"
#include "basic/histgram.h"
#include "graph/connection_graph.h"
#include "graph/contig_graph.h"
#include "graph/hash_graph.h"
#include "misc/hash_aligner.h"
#include "misc/options_description.h"
#include "misc/sequence_io.h"
#include "misc/utils.h"
#include "misc/log.h"
#include "sequence/sequence.h"
#include "sequence/short_sequence.h"
#include "assembly/local_assembler.h"


using namespace std;

struct IDBAOption
{
    string directory;
    string read_file;
    string long_read_file;
    int mink;
    int maxk;
    int step;
    int inner_mink;
    int inner_step;
    int prefix_length;
    int min_count;
    int min_support;
    int min_contig;
    double similar;
    int max_mismatch;
    int seed_kmer_size;
    int num_threads;
    int min_pairs;
    bool is_no_local;
    bool is_use_all;
    bool is_no_coverage;
    bool is_no_correct;
    bool is_pre_correction;

    IDBAOption()
    {
        directory = "out";
        mink = 20;
        maxk = 100;
        step = 20;
        inner_mink = 10;
        inner_step = 5;
        prefix_length = 3;
        min_count = 2;
        min_support = 1;
        min_contig = 200;
        similar = 0.95;
        max_mismatch = 3;
        seed_kmer_size = 30;
        num_threads = 0;
        min_pairs = 3;
        is_no_local = false;
        is_use_all = false;
        is_no_coverage = false;
        is_no_correct = false;
        is_pre_correction = false;
    }

    string log_file()
    { return directory + "/log"; }

    string kmer_file()
    { return directory + "/kmer"; }

    string align_file(int kmer_size)
    { return directory + FormatString("/align-%d", kmer_size); }

    string graph_file(int kmer_size)
    { return directory + FormatString("/graph-%d.fa", kmer_size); }

    string contig_file(int kmer_size)
    { return directory + FormatString("/contig-%d.fa", kmer_size); }

    string contig_info_file(int kmer_size)
    { return directory + FormatString("/contig-info-%d.fa", kmer_size); }

    string local_contig_file(int kmer_size)
    { return directory + FormatString("/local-contig-%d.fa", kmer_size); }

    string contig_file()
    { return directory + "/contig.fa"; }

    string scaffold_file()
    { return directory + "/scaffold.fa"; }
};

AssemblyInfo assembly_info;
IDBAOption option;
int median = 0;
int sd = 0;
int read_length = 0;

void BuildHashGraph(int kmer_size);
void Assemble(HashGraph &hash_graph);
void AlignReads(int kmer_size, bool is_all);
void CorrectReads(int kmer_size);
void EstimateDistance(int kmer_size);
void LocalAssembly(int kmer_size, int new_kmer_size);
void Iterate(int kmer_size, int new_kmer_size);
void Scaffold(int kmer_size, int min_contig);
bool CompareLength(const Sequence &x, const Sequence &y);

int main(int argc, char *argv[])
{
    OptionsDescription desc;
    
    desc.AddOption("out", "o", option.directory, "output directory");
    desc.AddOption("read", "r", option.read_file, FormatString("fasta read file (<=%d)", ShortSequence::max_size()));
    desc.AddOption("long_read", "l", option.long_read_file, FormatString("fasta long read file (>%d)", ShortSequence::max_size()));
    desc.AddOption("mink", "", option.mink, FormatString("minimum k value (<=%d)", Kmer::max_size()));
    desc.AddOption("maxk", "", option.maxk, FormatString("maximum k value (<=%d)", Kmer::max_size()));
    desc.AddOption("step", "", option.step, "increment of k-mer of each iteration");
    desc.AddOption("inner_mink", "", option.inner_mink, "inner minimum k value");
    desc.AddOption("inner_step", "", option.inner_step, "inner increment of k-mer");
    desc.AddOption("prefix", "", option.prefix_length, "prefix length used to build sub k-mer table");
    desc.AddOption("min_count", "", option.min_count, "minimum multiplicity for filtering k-mer when building the graph");
    desc.AddOption("min_support", "", option.min_support, "minimum supoort in each iteration");
    desc.AddOption("num_threads", "", option.num_threads, "number of threads");
    desc.AddOption("seed_kmer", "", option.seed_kmer_size, "seed kmer size for alignment");
    desc.AddOption("min_contig", "", option.min_contig, "min size of contig");
    desc.AddOption("similar", "", option.similar, "similarity for alignment");
    desc.AddOption("max_mismatch", "", option.max_mismatch, "max mismatch of error correction");
    desc.AddOption("min_pairs", "", option.min_pairs, "minimum number of pairs");
    desc.AddOption("no_local", "", option.is_no_local, "do not use local assembly");
    desc.AddOption("no_coverage", "", option.is_no_coverage, "do not iterate on coverage");
    desc.AddOption("no_correct", "", option.is_no_correct, "do not do correction");
    desc.AddOption("pre_correction", "", option.is_pre_correction, "perform pre-correction before assembly");

    try
    {
        desc.Parse(argc, argv);

        if (option.read_file == "" && option.long_read_file == "")
            throw logic_error("not enough parameters");

        if (option.maxk < option.mink)
            throw invalid_argument("mink is larger than maxk");

        if (option.maxk > (int)Kmer::max_size())
            throw invalid_argument("maxk is too large");
    }
    catch (exception &e)
    {
        cerr << e.what() << endl;
        cerr << "IDBA-UD - Iterative de Bruijn Graph Assembler for sequencing data with highly uneven depth." << endl;
        cerr << "Usage: idba_ud -r read.fa -o output_dir" << endl;
        cerr << "Allowed Options: " << endl;
        cerr << desc << endl;
        exit(1);
    }


    MakeDir(option.directory);

    LogThread log_thread(option.log_file());

    string begin_file = option.directory + "/begin";
    fclose(OpenFile(begin_file, "wb"));

    if (option.num_threads == 0)
        option.num_threads = omp_get_max_threads();
    else
        omp_set_num_threads(option.num_threads);
    cout << "number of threads " << option.num_threads << endl;

    ReadInput(option.read_file, option.long_read_file, assembly_info);
    cout << "reads " << assembly_info.reads.size() << endl;
    cout << "long reads " << assembly_info.long_reads.size() << endl;

    read_length = assembly_info.read_length();
    cout << "read_length " << read_length << endl;

    if (option.maxk > read_length)
        option.maxk = read_length;

    if (option.is_pre_correction)
    {
        int kmer_size = (option.maxk + option.mink)/2;
        cout << "kmer " << kmer_size << endl;
        BuildHashGraph(kmer_size);
        AlignReads(kmer_size, true);
        CorrectReads(kmer_size);
        assembly_info.ClearStatus();
    }

    int old_kmer_size = 0;
    int kmer_size = option.mink;
    while (true)
    {
        cout << "kmer " << kmer_size << endl;
        if (kmer_size == option.mink)
            BuildHashGraph(kmer_size);
        else
            Iterate(old_kmer_size, kmer_size);

        AlignReads(kmer_size, true);

        if (kmer_size < option.maxk)
        {
            CorrectReads(kmer_size);
            assembly_info.ClearStatus();

            old_kmer_size = kmer_size;
            kmer_size = min(option.maxk, kmer_size + option.step);
            LocalAssembly(old_kmer_size, kmer_size);
        }
        else
            break;
    }

    kmer_size = option.maxk;

    deque<Sequence> contigs;
    deque<string> names;
    ReadSequence(option.contig_file(kmer_size), contigs, names);
    FastaWriter writer(option.contig_file());
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() >= option.min_contig)
            writer.Write(contigs[i], names[i]);
    }

    Scaffold(option.maxk, option.min_contig);

    string end_file = option.directory + "/end";
    fclose(OpenFile(end_file, "wb"));

    fflush(stdout);

    return 0;
}

void BuildHashGraph(int kmer_size)
{
    BuildKmerFile(assembly_info, kmer_size, option.min_count, option.prefix_length, option.kmer_file());

    HashGraph hash_graph(kmer_size);
    ReadKmerFile(option.kmer_file(), hash_graph);

    hash_graph.RefreshEdges();
    InsertInternalKmers(assembly_info, hash_graph, option.min_count);

    Assemble(hash_graph);
}

void Assemble(HashGraph &hash_graph)
{
    cout << "kmers " << hash_graph.num_vertices() << " "<< hash_graph.num_edges() << endl;

    int kmer_size = hash_graph.kmer_size();
    double min_cover = max(1, (kmer_size == option.mink ? option.min_count : option.min_support));

    deque<Sequence> contigs;
    deque<ContigInfo> contig_infos;
    hash_graph.Assemble(contigs, contig_infos);
    hash_graph.clear();

    {
        HashGraph tmp_hash_graph;
        tmp_hash_graph.swap(hash_graph);
    }

    ContigGraph contig_graph(kmer_size, contigs, contig_infos);
    contigs.clear();
    contig_infos.clear();

    contig_graph.RemoveDeadEnd(option.min_contig);
    int bubble = contig_graph.RemoveBubble();
    cout << "merge bubble " << bubble << endl;

    if (!option.is_no_coverage)
        contig_graph.RemoveLocalLowCoverage(min_cover, option.min_contig, 0.1);

    contig_graph.SortVertices();
    contig_graph.GetContigs(contigs, contig_infos);
    WriteSequence(option.graph_file(kmer_size), contigs);
    contigs.clear();
    contig_infos.clear();

    if (!option.is_no_coverage)
    {
        double ratio = (kmer_size < option.maxk) ? 0.5 : 0.1;
        contig_graph.IterateLocalCoverage(option.min_contig, ratio, min_cover, 1e100, 1.1);
    }

    deque<Sequence> multi_contigs;
    deque<ContigInfo> multi_contig_infos;
    contig_graph.SortVertices();
    contig_graph.GetContigs(multi_contigs, multi_contig_infos);
    PrintN50(multi_contigs);
    WriteSequence(option.contig_file(kmer_size), multi_contigs);
    //WriteContigInfo(option.contig_info_file(kmer_size), multi_contig_infos);
}

void AlignReads(int kmer_size, bool is_all)
{
    deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    HashAligner hash_aligner(min(kmer_size, option.seed_kmer_size), option.min_contig, 2);
    hash_aligner.Initialize(contigs);

    int64_t num_aligned_reads = AlignReads(assembly_info, hash_aligner, option.similar, option.align_file(kmer_size), is_all);
    cout << "aligned " << num_aligned_reads << " reads" << endl;
}

void CorrectReads(int kmer_size)
{
    if (option.is_no_correct)
        return;

    deque<Sequence> contigs;
    deque<string> names;
    deque<ContigInfo> contig_infos;
    ReadSequence(option.contig_file(kmer_size), contigs, names);
    CorrectReads(assembly_info, contigs, contig_infos, option.align_file(kmer_size), option.max_mismatch);
    //WriteSequence(option.contig_file(kmer_size), contigs);
    WriteContig(option.contig_file(kmer_size), contigs, contig_infos, FormatString("contig-%d", kmer_size));
}

void EstimateDistance(int kmer_size)
{
    deque<ShortSequence> &reads = assembly_info.reads;
    Histgram<int> histgram;
    FILE *falign = OpenFile(option.align_file(kmer_size), "rb");

    for (unsigned i = 0; i < reads.size(); i += 2)
    {
        deque<HashAlignerRecord> records1;
        deque<HashAlignerRecord> records2;

        ReadHashAlignerRecords(falign, records1);
        ReadHashAlignerRecords(falign, records2);

        if (records1.size() == 1 && records2.size() == 1)
        {
            HashAlignerRecord record1 = records1[0];
            HashAlignerRecord record2 = records2[0];
            record2.ReverseComplement();

            if (record1.ref_id == record2.ref_id && record1.is_reverse == record2.is_reverse)
            {
                histgram.insert((record2.ref_to + record2.query_length - record2.query_to)
                        - (record1.ref_from - record1.query_from));
            }
        }
    }

    histgram.Trim(0.01);
    cout << "distance mean " << histgram.mean() << " sd " << histgram.sd() << endl;
    median = histgram.median();
    sd = (int)histgram.sd();

    fclose(falign);
}

void LocalAssembly(int kmer_size, int new_kmer_size)
{
    if (median == 0)
        EstimateDistance(kmer_size);

    deque<ShortSequence> &reads = assembly_info.reads;

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    LocalAssembler local_assembler;
    local_assembler.Initialize(assembly_info, contigs);
    local_assembler.set_num_threads(option.num_threads);
    local_assembler.set_mink(option.inner_mink);
    local_assembler.set_maxk(new_kmer_size);
    local_assembler.set_step(option.inner_step);
    local_assembler.set_min_contig(option.min_contig);
    local_assembler.set_insert_distance(median, sd);

    FILE *falign = OpenFile(option.align_file(kmer_size), "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);

        ReadHashAlignerRecordBlock(falign, all_records);
#pragma omp parallel for
        for (int i = 0; i < size; ++i)
        {
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0)
                local_assembler.AddReadByHashAlignerRecord(record, offset + i);
        }
    }
    fclose(falign);

    deque<Sequence> local_contigs;

    if (!option.is_no_local)
        local_assembler.Assemble(local_contigs);
    
    int num_seed_contigs = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        if ((int)contigs[i].size() > option.min_contig)
            ++num_seed_contigs;
    }

    cout << "seed contigs " << num_seed_contigs <<  " local contigs " << local_contigs.size() << endl;
    WriteSequence(option.local_contig_file(kmer_size), local_contigs, FormatString("local_contig_%d", kmer_size));
}

void Iterate(int kmer_size, int new_kmer_size)
{
    deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    deque<Sequence> local_contigs;
    ReadSequence(option.local_contig_file(kmer_size), local_contigs);

    deque<Sequence> multi_contigs;
    ReadSequence(option.graph_file(kmer_size), multi_contigs);

    uint64_t sum = 0;
    for (unsigned i = 0; i < contigs.size(); ++i)
        sum += contigs[i].size();
    HashGraph hash_graph(kmer_size);
    hash_graph.reserve(sum);

    deque<Sequence> old_contigs;
    old_contigs.insert(old_contigs.end(), contigs.begin(), contigs.end());
    old_contigs.insert(old_contigs.end(), local_contigs.begin(), local_contigs.end());
    old_contigs.insert(old_contigs.end(), multi_contigs.begin(), multi_contigs.end());
    contigs.clear();
    local_contigs.clear();
    multi_contigs.clear();

    IterateHashGraph(assembly_info, new_kmer_size, option.min_support, hash_graph, old_contigs);
    kmer_size = new_kmer_size;
    old_contigs.clear();

    if (kmer_size < option.maxk)
        hash_graph.RefreshEdges();
    else
        hash_graph.AddAllEdges();

    Assemble(hash_graph);
}

void Scaffold(int kmer_size, int min_contig)
{
    if (median == 0)
        EstimateDistance(kmer_size);

    deque<Sequence> contigs;
    ReadSequence(option.contig_file(kmer_size), contigs);

    deque<ContigInfo> contig_infos(contigs.size());
    vector<int> num_aligned_reads(contigs.size(), 0);
    vector<double> coverage(contigs.size());

    deque<ShortSequence> &reads = assembly_info.reads;

    FILE *falign = OpenFile(option.align_file(kmer_size), "rb");
    int buffer_size = (1 << 20) * option.num_threads;
    for (int64_t offset = 0; offset < (int64_t)reads.size(); offset += buffer_size)
    {
        int64_t size = min((int64_t)buffer_size, (int64_t)(reads.size() - offset));
        vector<HashAlignerRecord> all_records(size);

        ReadHashAlignerRecordBlock(falign, all_records);
#pragma omp parallel for
        for (int i = 0; i < size; ++i)
        {
            HashAlignerRecord &record = all_records[i];

            if (record.match_length != 0)
            {
#pragma omp atomic
                ++num_aligned_reads[record.ref_id];
            }
        }
    }
    fclose(falign);
//    FILE *falign = OpenFile(option.align_file(kmer_size), "rb");
//    for (unsigned i = 0; i < reads.size(); ++i)
//    {
//        deque<HashAlignerRecord> records;
//        ReadHashAlignerRecords(falign, records);
//
//        for (unsigned j = 0; j < records.size(); ++j)
//            ++num_aligned_reads[records[j].ref_id];
//    }

    double sum_coverage = 0;
    double sum_length = 0;
#pragma omp parallel for reduction(+: sum_coverage, sum_length)
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        if ((int)contigs[i].size() > option.min_contig)
        {
            sum_coverage += num_aligned_reads[i];
            sum_length += contigs[i].size() - read_length + 1;
            coverage[i] = 1.0 * num_aligned_reads[i] / (contigs[i].size() - reads[0].size() + 1);
            contig_infos[i].set_kmer_count(num_aligned_reads[i]);
        }
    }
    double mean_coverage = sum_coverage / sum_length;
    cout << "expected coverage " << mean_coverage << endl;

    ConnectionGraph connection_graph(kmer_size);
    connection_graph.Initialize(contigs, contig_infos);

    int num_connections = 0;
    falign = OpenFile(option.align_file(kmer_size), "rb");
    for (unsigned i = 0; i < reads.size(); i += 2)
    {
        deque<HashAlignerRecord> records1;
        deque<HashAlignerRecord> records2;
        ReadHashAlignerRecords(falign, records1);
        ReadHashAlignerRecords(falign, records2);

        if (records1.size() != 1 || records2.size() != 1)
            continue;

        for (unsigned j = 0; j < records1.size(); ++j)
        {
            for (unsigned k = 0; k < records2.size(); ++k)
            {
                HashAlignerRecord &r1 = records1[j];
                HashAlignerRecord &r2 = records2[k];
                r2.ReverseComplement();

                if (r1.ref_length > min_contig && r2.ref_length > min_contig
                        && r1.ref_from - r1.query_from > r1.ref_length - median - 3*sd
                        && r2.ref_to + r2.query_length - r2.query_to < median + 3*sd
                        && r1.ref_id != r2.ref_id
                        )
                {
                    ContigGraphVertexAdaptor from(&connection_graph.vertices()[r1.ref_id], r1.is_reverse);
                    ContigGraphVertexAdaptor to(&connection_graph.vertices()[r2.ref_id], r2.is_reverse);
                    int d = median - (r1.ref_length - (r1.ref_from - r1.query_from)) - (r2.ref_to + r2.query_length - r2.query_to);
                    connection_graph.AddConnection(from, to, d);
                    ++num_connections;
                }
            }
        }
    }

    connection_graph.ParseConnections(option.min_pairs);
    cout << "num edges " << connection_graph.num_edges() << endl;

    deque<Sequence> scaffolds;
    deque<ContigInfo> scaffold_infos;
    connection_graph.Scaffold(mean_coverage, median, scaffolds, scaffold_infos);
    PrintN50(scaffolds);

    WriteContig(option.scaffold_file(), scaffolds, scaffold_infos, "scaffold", option.min_contig);
//    sort(scaffolds.begin(), scaffolds.end(), CompareLength);
//    FastaWriter writer(option.scaffold_file());
//
//    for (unsigned i = 0; i < scaffolds.size(); ++i)
//    {
//        if ((int)scaffolds[i].size() >= option.min_contig)
//            writer.Write(scaffolds[i], FormatString("scaffold_%d", i));
//    }
}

bool CompareLength(const Sequence &x, const Sequence &y)
{
    return x.size() > y.size();
}

