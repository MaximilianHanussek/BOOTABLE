/**
 * @file connection_graph.h
 * @brief ConnectionGraph Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.1
 * @date 2011-09-02
 */

#ifndef __GRAPH_CONNECTION_GRAPH_H_

#define __GRAPH_CONNECTION_GRAPH_H_

#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

#include "graph/contig_graph.h"
#include "graph/contig_graph_properties.h"

class ContigConnection
{
public:
    ContigGraphVertexAdaptor to;
    int distance;
    std::vector<int> values;

    ContigConnection() { }

    ContigConnection(const ContigGraphVertexAdaptor &node, int value)
    { to = node; values.push_back(value); }

    static bool CompareVertex(const ContigConnection &x, const ContigConnection &y)
    { return x.to < y.to; }

    bool operator <(const ContigConnection &connection) const
    { return distance < connection.distance; }

    bool AddValue(const ContigGraphVertexAdaptor &node, int value)
    {
        if (node == to)
        {
            values.push_back(value);
            return true;
        }
        return false;
    }

    void Parse()
    {
        std::sort(values.begin(), values.end());
        if (values.empty())
            distance = 0;
        else
            distance = values[values.size()/2];
    }

    void Translate(int offset);
    bool IsConsistant(int delta);
    void RemoveInconsistantPairs(int delta);
};

class ConnectionGraph: public ContigGraph
{
public:
    typedef ContigGraphProperties<ContigConnection>::list_type connection_list_type;
    typedef connection_list_type::iterator connection_list_iterator;

    ConnectionGraph(uint32_t kmer_size)
        : ContigGraph(kmer_size) {}

    void Initialize(std::deque<Sequence> &contigs)
    { ContigGraph::Initialize(contigs); connectons_.resize(contigs.size()); }

    void Initialize(std::deque<Sequence> &contigs, std::deque<ContigInfo> &contig_infos)
    { ContigGraph::Initialize(contigs, contig_infos); connectons_.resize(contigs.size()); }

    void AddConnection(ContigGraphVertexAdaptor from, ContigGraphVertexAdaptor to, int distance)
    {
        AddConnection(connectons_[from], to, distance);
        from.ReverseComplement();
        to.ReverseComplement();
        std::swap(from, to);
        AddConnection(connectons_[from], to, distance);
    }

    void AddConnection(connection_list_type &connection_list, ContigGraphVertexAdaptor to, int distance)
    {
        for (connection_list_iterator p = connection_list.begin(); p != connection_list.end(); ++p)
        {
            if (p->AddValue(to, distance))
                return;
        }
        connection_list.push_front(ContigConnection(to, distance));
    }

    void RemoveConnection(ContigGraphVertexAdaptor from, ContigGraphVertexAdaptor to)
    {
        RemoveConnection(connectons_[from], to);
        from.ReverseComplement();
        to.ReverseComplement();
        std::swap(from, to);
        RemoveConnection(connectons_[from], to);
    }

    void RemoveConnection(connection_list_type &connection_list, ContigGraphVertexAdaptor to)
    {
        for (connection_list_iterator p = connection_list.begin(); p != connection_list.end(); ++p)
        {
            if (p->to == to)
            {
                connection_list.erase(p);
                break;
            }
        }
    }

    void RemoveAllConnections(ContigGraphVertexAdaptor current)
    {
        std::deque<ContigConnection> tmp;
        for (connection_list_iterator p = connections()[current].begin(); p != connections()[current].end(); ++p)
            tmp.push_back(*p);

        for (unsigned i = 0; i < tmp.size(); ++i)
            RemoveConnection(current, tmp[i].to);
    }

    void ParseConnections(int min_pairs = 0);
    void ParseConnections(int read_length, int insert_distance, double expected_coverage);

    bool IsConnected(ContigGraphVertexAdaptor from, ContigGraphVertexAdaptor to);
    int64_t RemoveTransitiveConnections();

    int64_t RemoveRepeat(double expected, int insert_distance);
    int64_t RemoveTips(int insert_distance);

    int Overlap(const Sequence &a, const Sequence &b);
    int64_t RemoveNonSupportConnections();

    bool IsConverged(ContigGraphVertexAdaptor from, ContigGraphVertexAdaptor to, std::deque<ContigGraphVertexAdaptor> &qu);
    int64_t MergeSimilarPath();

    int64_t LinearPaths(std::deque<ContigGraphPath> &paths);
    void Linearlize();

    int64_t GetContigs(std::deque<Sequence> &contigs, std::deque<ContigInfo> &contig_infos);
    int64_t Scaffold(double mean_coverage, double median, std::deque<Sequence> &scaffolds, std::deque<ContigInfo> &scaffold_infos);

    void GetComponents(std::deque<std::deque<ContigGraphVertexAdaptor> > &components, std::deque<std::string> &component_strings);

    ContigGraphProperties<ContigConnection> &connections() { return connectons_; }

    uint32_t num_edges() { return num_edges_; }

    void swap(ConnectionGraph &connection_graph)
    { ContigGraph::swap(connection_graph); connectons_.swap(connection_graph.connectons_); std::swap(num_edges_, connection_graph.num_edges_); }

    static const int kTimeLimit = 500;

private:
    uint32_t num_edges_;
    ContigGraphProperties<ContigConnection> connectons_;
};

#endif

