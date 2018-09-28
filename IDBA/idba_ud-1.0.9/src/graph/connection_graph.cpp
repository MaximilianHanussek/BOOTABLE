/**
 * @file connection_graph.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.1
 * @date 2011-09-02
 */

#include "graph/connection_graph.h"

#include <deque>
#include <map>
#include <iostream>
#include <sstream>

using namespace std;

void ConnectionGraph::ParseConnections(int min_pairs)
{
    deque<ContigGraphVertexAdaptor> removed_from;
    deque<ContigGraphVertexAdaptor> removed_to;
    num_edges_ = 0;
    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices()[i], strand);

            connection_list_type &connection_list = connectons_[current];
            connection_list_iterator p = connection_list.begin();
            int maximum = 0;
            while (p != connection_list.end())
            {
                if ((int)p->values.size() > maximum)
                    maximum = p->values.size();
                //total += p->values.size();
                ++p;
            }

            vector<ContigConnection> neighbors;

            p = connection_list.begin();
            while (p != connection_list.end())
            {
                p->Parse();

                neighbors.push_back(*p);

                if ((int)p->values.size() < min_pairs || (int)p->values.size() < maximum * 0.05)
                {
                    removed_from.push_back(current);
                    removed_to.push_back(p->to);
                }
                ++p;
            }

//            sort(neighbors.begin(), neighbors.end(), ContigConnection::CompareVertex);
//            for (unsigned i = 0; i+1 < neighbors.size(); ++i)
//            {
//                if (neighbors[i].to.id() == neighbors[i+1].to.id())
//                {
//                    cout << current.id() << " " << neighbors[i].to.id() << " " << neighbors[i].to.is_reverse()
//                        << " " << neighbors[i+1].to.id() <<  " " << neighbors[i+1].to.is_reverse() << endl;
//                    cout << neighbors[i].values.size() << " " << neighbors[i+1].values.size() << " " << neighbors[i].to.contig_size() << endl;
//                    if (neighbors[i].values.size() < neighbors[i+1].values.size() * 0.1)
//                    {
//                        removed_from.push_back(current);
//                        removed_to.push_back(neighbors[i].to);
//                    }
//
//                    if (neighbors[i+1].values.size() < neighbors[i].values.size() * 0.1)
//                    {
//                        removed_from.push_back(current);
//                        removed_to.push_back(neighbors[i+1].to);
//                    }
//                }
//            }
//
            num_edges_ += connection_list.size();
        }
    }

    for (unsigned i = 0; i < removed_from.size(); ++i)
        RemoveConnection(removed_from[i], removed_to[i]);

    num_edges_ /= 2;
}

void ConnectionGraph::ParseConnections(int read_length, int insert_distance, double expected_coverage)
{
    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices()[i], strand);

            connection_list_type &connection_list = connectons_[current];
            connection_list_iterator p = connection_list.begin();
            while (p != connection_list.end())
            {
                p->Parse();

                int from = std::max(-insert_distance, -(int)current.contig_size()) + insert_distance;
                int to = insert_distance;

                from = std::max(from, p->distance);
                to = std::min(to, p->distance + (int)current.contig_size());

                if (p->values.size() < 0.1 * expected_coverage * (to - from - read_length + 1)
                        || p->values.size() < 2)
                {
                    //std::cout << p->values.size() << " " << 0.1 * expected_coverage * (to - from - read_length + 1) << std::endl;
                    connection_list.erase(p++);
                }
                else
                    ++p;
            }
        }
    }
}

bool ConnectionGraph::IsConnected(ContigGraphVertexAdaptor from, ContigGraphVertexAdaptor to)
{
    deque<ContigGraphVertexAdaptor> qu;
    qu.push_back(from);
    from.status().SetUsedFlag();

    int time = 0;
    int index = 0;
    bool is_found = false;
    while (++time < kTimeLimit && index < (int)qu.size() && !is_found)
    {
        ContigGraphVertexAdaptor current = qu[index++];

        for (connection_list_iterator p = connections()[current].begin(); p != connections()[current].end(); ++p)
        {
            if (p->to == to)
            {
                is_found = true;
                break;
            }

            if (!p->to.status().IsUsed())
            {
                qu.push_back(p->to);
                p->to.status().SetUsedFlag();
            }
        }
    }

    for (unsigned i = 0; i < qu.size(); ++i)
        qu[i].status().clear();

    return is_found;
}

int64_t ConnectionGraph::RemoveTransitiveConnections()
{
    int removed = 0;
    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices()[i], strand);
            connection_list_type &connection_list = connections()[current];

            if (connection_list.size() < 2)
                continue;

            deque<ContigConnection> connections;
            for (connection_list_iterator p = connection_list.begin(); p != connection_list.end(); ++p)
                connections.push_back(*p);

            for (unsigned j = 0; j < connections.size(); ++j)
            {
                RemoveConnection(current, connections[j].to);
                if (!IsConnected(current, connections[j].to))
                {
                    for (unsigned k = 0; k < connections[j].values.size(); ++k)
                    {
                        AddConnection(current, connections[j].to, connections[j].values[k]);
                    }
                }
                else
                    ++removed;
            }
        }
    }

    return removed;
}

int64_t ConnectionGraph::RemoveRepeat(double expected, int insert_distance)
{
    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        if ((int)vertices()[i].contig().size() > insert_distance)
            continue;

        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices()[i], strand);
            connection_list_type &connection_list = connections()[current];

            if (connection_list.size() < 2)
                continue;

            deque<ContigConnection> connections;
            for (connection_list_iterator p = connection_list.begin(); p != connection_list.end(); ++p)
                connections.push_back(*p);

            for (unsigned j = 0; j < connections.size(); ++j)
            {
                RemoveConnection(current, connections[j].to);
                if (!IsConnected(current, connections[j].to))
                {
                    for (unsigned k = 0; k < connections[j].values.size(); ++k)
                    {
                        AddConnection(current, connections[j].to, connections[j].values[k]);
                    }
                }
            }
        }
    }

    int removed = 0;
    for (unsigned i = 0; i < vertices().size(); ++i)
    {
//        if (vertices()[i].kmer_count() / (vertices()[i].contig().size() - kmer_size() + 1) <= 1.5 * expected)
//            continue;

        if ((int)vertices()[i].contig().size() > insert_distance)
            continue;

        ContigGraphVertexAdaptor node(&vertices()[i], false);
        ContigGraphVertexAdaptor rev_comp(&vertices()[i], true);
        
        if (connections()[node].size() > 1 || connections()[rev_comp].size() > 1)
        {
            for (int strand = 0; strand < 2; ++strand)
            {
                ContigGraphVertexAdaptor current(&vertices()[i], strand);
                connection_list_type &connection_list = connections()[current];
    //
    //            if (connection_list.size() < 2)
    //                continue;
    //
                deque<ContigConnection> connections;
                for (connection_list_iterator p = connection_list.begin(); p != connection_list.end(); ++p)
                    connections.push_back(*p);

                for (unsigned j = 0; j < connections.size(); ++j)
                {
                    RemoveConnection(current, connections[j].to);
                    ++removed;
                }
            }
        }
    }

    return removed;
}

int64_t ConnectionGraph::RemoveTips(int insert_distance)
{
    int removed = 0;
    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices()[i], strand);
            connection_list_type &connection_list = connections()[current];

            if ((int)current.contig_size() >= insert_distance)
                continue;

            if (connection_list.size() > 1)
                continue;

            ContigGraphVertexAdaptor rev_comp = current;
            rev_comp.ReverseComplement();
            if (connections()[rev_comp].size() != 0)
                continue;

            ++removed;
            RemoveAllConnections(current);
        }
    }

    return removed;
}

int ConnectionGraph::Overlap(const Sequence &a, const Sequence &b)
{
    for (int len = 100; len > 0; --len)
    {
        int mismatch = 0;
        for (int j = 0; j < len; ++j)
        {
            if (a[a.size()-len+j] != b[j])
                ++mismatch;

            if (mismatch > 3)
                break;
        }

        if (mismatch <= 3)
            return len;
    }

    return 0;
}

int64_t ConnectionGraph::RemoveNonSupportConnections()
{
    int removed = 0;
    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices()[i], strand);
            connection_list_type &connection_list = connections()[current];
            deque<ContigConnection> tmp;
            for (connection_list_iterator p = connection_list.begin(); p != connection_list.end(); ++p)
                tmp.push_back(*p);

            for (unsigned j = 0; j < tmp.size(); ++j)
            {
                Sequence a = current.contig();
                Sequence b = tmp[j].to.contig();

                int overlap = Overlap(a, b);
                if (overlap < -tmp[j].distance - 50)
                {
                    ++removed;
                    RemoveConnection(current, tmp[j].to);
                }
            }
        }
    }

    return removed;
}

bool ConnectionGraph::IsConverged(ContigGraphVertexAdaptor from, ContigGraphVertexAdaptor to, deque<ContigGraphVertexAdaptor> &qu)
{
    qu.clear();

    qu.push_back(from);
    from.status().SetUsedFlag();

    int index = 0;
    int time = 0;
    bool is_deadend = false;
    while (index < (int)qu.size() && ++time < kTimeLimit)
    {
        ContigGraphVertexAdaptor current = qu[index++];

        if (current == to)
            continue;

        if (connections()[current].empty())
        {
            is_deadend = true;
            break;
        }

        for (connection_list_iterator p = connections()[current].begin(); p != connections()[current].end(); ++p)
        {
            if (!p->to.status().IsUsed())
            {
                qu.push_back(p->to);
                p->to.status().SetUsedFlag();
            }
        }
    }

    for (unsigned i = 0; i < qu.size(); ++i)
        qu[i].status().clear();

    if (time < kTimeLimit && !is_deadend && index == (int)qu.size())
        return true;
    else
        return false;
}

int64_t ConnectionGraph::MergeSimilarPath()
{
    int num_converged_paths = 0;

    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices()[i], strand);
            if (connections()[current].size() < 2)
                continue;

            deque<ContigGraphVertexAdaptor> qu;
            ContigGraphPath path;
            path.Append(current, 0);

            bool is_converged = false;
            for (int len = 2; len <= 5; ++len)
            {
                if (connections()[path.back()].size() == 0)
                    break;

                ContigConnection connection = connections()[path.back()].front();
                path.Append(connection.to, connection.distance);

                if (IsConverged(path.front(), path.back(), qu))
                {
                    path.ReverseComplement();
                    if (IsConverged(path.front(), path.back(), qu))
                    {
                        path.ReverseComplement();
                        is_converged = true;
                        break;
                    }
                    path.ReverseComplement();
                }
            }

            if (is_converged)
            {
                num_converged_paths++;

                for (unsigned j = 0; j < qu.size(); ++j)
                {
                    ContigGraphVertexAdaptor x = qu[j];

                    if (x != path.back())
                    {
                        RemoveAllConnections(x);
                    }

                    if (x != path.front())
                    {
                        x.ReverseComplement();
                        RemoveAllConnections(x);
                    }
                }

                for (unsigned j = 0; j+1 < path.num_nodes(); ++j)
                {
                    for (int k = 0; k < 10; ++k)
                        AddConnection(path[j], path[j+1], path.distances()[j]);
                }
            }
        }
    }

    return num_converged_paths;
}

int64_t ConnectionGraph::LinearPaths(deque<ContigGraphPath> &paths)
{
    paths.clear();

    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        if (vertices()[i].status().IsUsed())
            continue;

        ContigGraphVertexAdaptor current(&vertices()[i]);
        ContigGraphPath path;
        path.Append(current, 0);
        current.status().SetUsedFlag();

        for (int strand = 0; strand < 2; ++strand)
        {
            while (true)
            {
                ContigGraphVertexAdaptor current = path.back();

                if (connections()[current].size() != 1)
                    break;

                ContigGraphVertexAdaptor next = connections()[current].front().to;
                int d = connections()[current].front().distance;
                ContigGraphVertexAdaptor rev_next = next;
                rev_next.ReverseComplement();

                if (connections()[rev_next].size() != 1)
                    break;
                
                if (next.status().IsUsed())
                    break;

                path.Append(next, d);
                next.status().SetUsedFlag();
            }

            path.ReverseComplement();
        }

        paths.push_back(path);
    }

    ClearStatus();

    return paths.size();
}

void ConnectionGraph::Linearlize()
{
    deque<ContigGraphPath> paths;
    LinearPaths(paths);

    deque<Sequence> contigs(paths.size());
    deque<ContigInfo> contig_infos(paths.size());
    for (unsigned i = 0; i < paths.size(); ++i)
        paths[i].Assemble(contigs[i], contig_infos[i]);

    ConnectionGraph new_connection_graph(kmer_size());
    new_connection_graph.Initialize(contigs, contig_infos);

    map<ContigGraphVertexAdaptor, ContigConnection> super_contig_map;
    for (unsigned i = 0; i < contigs.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigConnection connection;
            connection.to = ContigGraphVertexAdaptor(&new_connection_graph.vertices()[i], strand);
            connection.distance = 0;
            for (unsigned j = 0; j < paths[i].num_nodes(); ++j)
            {
                super_contig_map[paths[i][j]] = connection;
                if (j+1 < paths[i].num_nodes())
                    connection.distance += paths[i][j].contig_size() + paths[i].distances()[j];
            }

            paths[i].ReverseComplement();
        }
    }

    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices()[i], strand);
            for (connection_list_iterator p = connections()[current].begin(); p != connections()[current].end(); ++p)
            {
                ContigGraphVertexAdaptor from = current;
                ContigGraphVertexAdaptor to = p->to;

                ContigGraphVertexAdaptor new_from = super_contig_map[from].to;
                ContigGraphVertexAdaptor new_to = super_contig_map[to].to;
                int distance = p->distance - (new_from.contig_size() - (super_contig_map[from].distance + from.contig_size()))
                    - (super_contig_map[to].distance);

                if (new_from.id() != new_to.id())
                {
                    for (unsigned k = 0; k < p->values.size(); ++k)
                        new_connection_graph.AddConnection(new_connection_graph.connections()[new_from], new_to, distance);
                }
            }
        }
    }

    swap(new_connection_graph);
}

int64_t ConnectionGraph::GetContigs(deque<Sequence> &contigs, deque<ContigInfo> &contig_infos)
{
    contigs.clear();
    contig_infos.clear();
    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        contigs.push_back(vertices()[i].contig());
        contig_infos.push_back(vertices()[i].contig_info());
    }
    return contigs.size();
}

int64_t ConnectionGraph::Scaffold(double mean_coverage, double median, deque<Sequence> &scaffolds, deque<ContigInfo> &scaffold_infos)
{
    scaffolds.clear();
    for (int round = 0; round < 3; ++round)
    {
        Linearlize();
        ParseConnections();

        RemoveRepeat(mean_coverage, median);
        ParseConnections();

        RemoveTransitiveConnections();
        ParseConnections();

        RemoveTips(median);
        ParseConnections();
    }

    GetContigs(scaffolds, scaffold_infos);
    return scaffolds.size();
}

void ConnectionGraph::GetComponents(deque<deque<ContigGraphVertexAdaptor> > &components, deque<string> &component_strings)
{
    components.clear();
    component_strings.clear();

    for (unsigned i = 0; i < vertices().size(); ++i)
    {
        if (vertices()[i].status().IsUsed())
            continue;

        deque<ContigGraphVertexAdaptor> qu;
        qu.push_back(ContigGraphVertexAdaptor(&vertices()[i], 0));
        vertices()[i].status().SetUsedFlag();

        stringstream ss;
        for (int index = 0; index < (int)qu.size(); ++index)
        {
            ContigGraphVertexAdaptor current = qu[index];

            for (int strand = 0; strand < 2; ++strand)
            {
                for (connection_list_iterator p = connections()[current].begin(); p != connections()[current].end(); ++p)
                {
                    ContigGraphVertexAdaptor next = p->to;
//                    if (next.status().IsUsed())
//                        continue;

                    if (strand == 0)
                    {
                        ss << current.id() << "_" << current.contig_size() << " " << next.id() << "_" << next.contig_size() <<  " " << current.contig_size() << " " << next.contig_size() << " " << p->distance << endl;

                        if (!next.status().IsUsed())
                            qu.push_back(next);
                    }
                    else
                    {
                        ss << next.id() << "_" << next.contig_size() << " " << current.id() << "_" << current.contig_size() << " " << next.contig_size() << " " << current.contig_size() << " " << p->distance <<  endl;

                        if (!next.status().IsUsed())
                            qu.push_back(next.ReverseComplement());
                    }

                    next.status().SetUsedFlag();
                }

                current.ReverseComplement();
            }
        }

        components.push_back(qu);
        component_strings.push_back(ss.str());
    }

    ClearStatus();
}

