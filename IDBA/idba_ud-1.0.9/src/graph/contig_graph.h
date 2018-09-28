/**
 * @file contig_graph.h
 * @brief ContigGraph Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-16
 */

#ifndef __GRAPH_CONTIG_GRAPH_H_

#define __GRAPH_CONTIG_GRAPH_H_

#include <deque>
#include <algorithm>

#include "sequence/sequence.h"
#include "basic/kmer.h"
#include "container/hash_map.h"
#include "basic/bit_operation.h"

#include "graph/contig_graph_properties.h"
#include "graph/contig_graph_vertex.h"
#include "graph/contig_graph_path.h"
#include "graph/contig_info.h"

class ContigGraph
{
public:
    explicit ContigGraph(uint32_t kmer_size = 0)
        : num_edges_(0), kmer_size_(kmer_size)
    {}
    explicit ContigGraph(uint32_t kmer_size, const std::deque<Sequence> &contigs)
        : num_edges_(0), kmer_size_(kmer_size)
    { Initialize(contigs); }

    explicit ContigGraph(uint32_t kmer_size, const std::deque<Sequence> &contigs,
            const std::deque<ContigInfo> &contig_infos)
        : num_edges_(0), kmer_size_(kmer_size)
    { Initialize(contigs, contig_infos); }

    ~ContigGraph() { clear(); }

    void Initialize(const std::deque<Sequence> &contigs)
    {
        std::deque<ContigInfo> contig_infos(contigs.size());
        Initialize(contigs, contig_infos);
    }

    void Initialize(const std::deque<Sequence> &contigs, const std::deque<ContigInfo> &contig_infos);

    void Refresh();
    void RefreshVertices();
    void RefreshEdges();

    void AddEdge(ContigGraphVertexAdaptor from, ContigGraphVertexAdaptor to)
    {
        from.out_edges().Add(to.contig()[kmer_size_-1]);
        from.ReverseComplement();
        to.ReverseComplement();
        std::swap(from, to);
        from.out_edges().Add(to.contig()[kmer_size_-1]);
    }

    void AddAllEdges();
    void ClearStatus();

    void MergeSimplePaths();
    int64_t Trim(int min_length);
    int64_t RemoveDeadEnd(int min_length);
    int64_t RemoveBubble();

    double IterateCoverage(int min_length, double min_cover, double max_cover, double factor = 1.1);
    double IterateLocalCoverage(int min_length, double ratio, double min_cover, double max_cover, double factor = 1.1);

    bool RemoveLowCoverage(double min_cover, int min_length);
    bool RemoveLocalLowCoverage(double min_cover, int min_length, double ratio);
    double LocalCoverage(ContigGraphVertexAdaptor current, int region_length);
    double LocalCoverageSingle(ContigGraphVertexAdaptor current, int region_length, double &num_count, int &num_kmer);

    int64_t Assemble(std::deque<Sequence> &contigs, std::deque<ContigInfo> &contig_infos);

    ContigGraphVertexAdaptor GetNeighbor(const ContigGraphVertexAdaptor &current, int x)
    {
        Kmer kmer = current.end_kmer(kmer_size_);
        kmer.ShiftAppend(x);
        return FindVertexAdaptorByBeginKmer(kmer);
    }

    void GetComponents(std::deque<std::deque<ContigGraphVertexAdaptor> > &components, std::deque<std::string> &component_strings);

    void SortVertices()
    { std::sort(vertices_.begin(), vertices_.end(), CompareContigLength); Refresh(); }

    void GetContigs(std::deque<Sequence> &contigs, std::deque<ContigInfo> &contig_infos);

    std::deque<ContigGraphVertex> &vertices() { return vertices_; }
    const std::deque<ContigGraphVertex> &vertices() const { return vertices_; }

    void swap(ContigGraph &contig_graph)
    {
        begin_kmer_map_.swap(contig_graph.begin_kmer_map_);
        vertices_.swap(contig_graph.vertices_);
        std::swap(num_edges_, contig_graph.num_edges_);
        std::swap(kmer_size_, contig_graph.kmer_size_);
    }

    uint32_t kmer_size() const { return kmer_size_; }
    void set_kmer_size(uint32_t kmer_size) { kmer_size_ = kmer_size; }

    uint64_t num_vertices() const { return vertices_.size(); }
    uint64_t num_edges() const { return num_edges_; }

    void clear()
    {
        num_edges_ = 0;
        vertices_.clear();
        begin_kmer_map_.clear();
        in_kmer_count_table_.clear();
    }

private:
    ContigGraph(const ContigGraph &);
    const ContigGraph &operator =(const ContigGraph &);

    static bool CompareContigLength(const ContigGraphVertex &x, const ContigGraphVertex &y)
    { return x.contig_size() > y.contig_size(); }

    void BuildBeginKmerMap();

    bool GetNextVertexAdaptor(ContigGraphVertexAdaptor &current, ContigGraphVertexAdaptor &next)
    {
        if (current.out_edges().size() != 1)
            return false;

        next = GetNeighbor(current, bit_operation::BitToIndex(current.out_edges()));
        return next.in_edges().size() == 1 && !(next.contig_size() == kmer_size_ && next.contig().IsPalindrome());
    }

    bool IsLoop(const ContigGraphPath &path, const ContigGraphVertexAdaptor &next)
    { return path.front().id() == next.id() || path.back().id() == next.id(); }

    ContigGraphVertexAdaptor FindVertexAdaptorByBeginKmer(const Kmer &begin_kmer)
    {
        Kmer key = begin_kmer.unique_format();

        HashMap<Kmer, uint32_t>::iterator iter = begin_kmer_map_.find(key);
        if (iter != begin_kmer_map_.end())
        {
            ContigGraphVertexAdaptor current(&vertices_[iter->second]);
            if (current.begin_kmer(kmer_size_) == begin_kmer)
                return current;
            current.ReverseComplement();
            if (current.begin_kmer(kmer_size_) == begin_kmer)
                return current;
        }

        return ContigGraphVertexAdaptor();
    }

    HashMap<Kmer, uint32_t> begin_kmer_map_;
    std::deque<ContigGraphVertex> vertices_;
    uint64_t num_edges_;
    uint32_t kmer_size_;

    HashMap<Kmer, uint32_t> in_kmer_count_table_;
};

#endif

