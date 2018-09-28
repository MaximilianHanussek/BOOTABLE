/**
 * @file contig_graph.cpp
 * @brief 
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-26
 */

#include "graph/contig_graph.h"

#include <deque>

#include "sequence/sequence.h"
#include "graph/contig_graph_branch_group.h"

#include <sstream>
#include <iostream>
#include <map>

using namespace std;

void ContigGraph::Initialize(const deque<Sequence> &contigs, const deque<ContigInfo> &contig_infos)
{
    vertices_.clear();
    vertices_.resize(contigs.size());
#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)contigs.size(); ++i)
    {
        vertices_[i].clear();
        vertices_[i].set_contig(contigs[i]);
        vertices_[i].set_contig_info(contig_infos[i]);
        vertices_[i].set_id(i);
    }
    RefreshEdges();
}

void ContigGraph::Refresh()
{
    RefreshVertices();
    RefreshEdges();
}

void ContigGraph::RefreshVertices()
{
    uint64_t index = 0;
    for (unsigned i = 0; i < vertices_.size(); ++i)
    {
        if (!vertices_[i].status().IsDead())
        {
            vertices_[index].swap(vertices_[i]);
            vertices_[index].set_id(index);
            ++index;
        }
    }
    vertices_.resize(index);
}

void ContigGraph::RefreshEdges()
{
    BuildBeginKmerMap();

    uint64_t total_degree = 0;
#pragma omp parallel for reduction(+: total_degree)
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices_[i], strand);

            for (int x = 0; x < 4; ++x)
            {
                if (current.out_edges()[x])
                {
                    Kmer kmer = current.end_kmer(kmer_size_);
                    kmer.ShiftAppend(x);
                    if (FindVertexAdaptorByBeginKmer(kmer).is_null())
                        current.out_edges().Remove(x);
                }
            }

//#pragma omp atomic
            total_degree += current.out_edges().size();
        }

        if (vertices_[i].contig().size() == kmer_size_ 
                && vertices_[i].contig().IsPalindrome())
        {
            vertices_[i].in_edges() = vertices_[i].out_edges() | vertices_[i].out_edges(); 
            vertices_[i].out_edges() = vertices_[i].in_edges();
        }
    }

    num_edges_ = total_degree / 2;
}

void ContigGraph::AddAllEdges()
{
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        vertices_[i].in_edges() = 15;
        vertices_[i].out_edges() = 15;
    }
    RefreshEdges();
}

void ContigGraph::ClearStatus()
{
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
        vertices_[i].status().clear();
}

int64_t ContigGraph::Trim(int min_length)
{
    uint64_t old_num_vertices = vertices_.size();
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        if (vertices_[i].contig().size() == kmer_size_
                && vertices_[i].contig().IsPalindrome())
            continue;

        if ((vertices_[i].in_edges().empty() || vertices_[i].out_edges().empty())
                && vertices_[i].contig().size() < min_length + kmer_size_ - 1
                && (vertices_[i].in_edges().size() + vertices_[i].out_edges().size() <= 1)
           )
        {
            vertices_[i].status().SetDeadFlag();
        }
    }
    Refresh();
    MergeSimplePaths();

    return old_num_vertices - vertices_.size();
}

int64_t ContigGraph::RemoveDeadEnd(int min_length)
{
    uint64_t num_deadend = 0;
    int l = 1;
    while (true)
    {
        l = min(2*l, min_length);
        num_deadend += Trim(l);

        if (l == min_length)
            break;
    }
    num_deadend += Trim(min_length);
    return num_deadend;
}

int64_t ContigGraph::RemoveBubble()
{
    deque<ContigGraphVertexAdaptor> candidates;
    omp_lock_t bubble_lock;
    omp_init_lock(&bubble_lock);

#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices_[i], strand);

            if (current.out_edges().size() > 1 && current.contig_size() > kmer_size_)
            {
                ContigGraphBranchGroup branch_group(this, current, 4, kmer_size_ + 2);

                if (branch_group.Search())
                {
                    ContigGraphVertexAdaptor begin = branch_group.begin();
                    ContigGraphVertexAdaptor end = branch_group.end();

                    begin.ReverseComplement();
                    end.ReverseComplement();
                    std::swap(begin, end);
                    ContigGraphBranchGroup rev_branch_group(this, begin, 4, kmer_size_ + 2);

                    if (rev_branch_group.Search() && rev_branch_group.end() == end)
                    {
                        omp_set_lock(&bubble_lock);
                        candidates.push_back(current);
                        omp_unset_lock(&bubble_lock);
                    }
                }
            }
        }
    }

    int64_t bubble = 0;
    for (unsigned i = 0; i < candidates.size(); ++i)
    {
        ContigGraphVertexAdaptor current = candidates[i];

        if (current.out_edges().size() > 1)
        {
            ContigGraphBranchGroup branch_group(this, current, 4, kmer_size_ + 2);

            if (branch_group.Search())
            {
                ContigGraphVertexAdaptor begin = branch_group.begin();
                ContigGraphVertexAdaptor end = branch_group.end();

                begin.ReverseComplement();
                end.ReverseComplement();
                std::swap(begin, end);
                ContigGraphBranchGroup rev_branch_group(this, begin, 4, kmer_size_ + 2);

                if (rev_branch_group.Search() && rev_branch_group.end() == end)
                {
                    branch_group.Merge();
                    ++bubble;
                }
            }
        }
    }

    Refresh();
    MergeSimplePaths();

    return bubble;
}

double ContigGraph::IterateCoverage(int min_length, double min_cover, double max_cover, double factor)
{
    min_cover = min(min_cover, max_cover);
    while (true)
    {
        RemoveLowCoverage(min_cover, min_length);
        min_cover *= factor;
        if (min_cover >= max_cover)
            break;
    }
    return min_cover;
}

double ContigGraph::IterateLocalCoverage(int min_length, double ratio, double min_cover, double max_cover, double factor)
{
    in_kmer_count_table_.reserve(vertices_.size());

    min_cover = min(min_cover, max_cover);
    while (true)
    {
        bool is_changed = RemoveLocalLowCoverage(min_cover, min_length, ratio);

        if (!is_changed)
            break;

        if (min_cover >= max_cover)
            break;

        min_cover *= factor;
    }
    return min_cover;
}

bool ContigGraph::RemoveLowCoverage(double min_cover, int min_length)
{
    bool is_changed = false;

#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        ContigGraphVertexAdaptor current(&vertices_[i]);

        if (current.contig_size() < min_length + kmer_size_ - 1
//                && ((current.in_edges().size() <= 1 && current.out_edges().size() <= 1)
//                        || current.in_edges().size() == 0 || current.out_edges().size() == 0)
           )
        {
            if (current.coverage() < min_cover)
            {
                is_changed = true;
                current.status().SetDeadFlag();
            }
        }
    }

    Refresh();
    //Trim(min_length);
    MergeSimplePaths();

    return is_changed;
}

bool ContigGraph::RemoveLocalLowCoverage(double min_cover, int min_length, double ratio)
{
    int region_length = 1000;
    bool is_changed = false;
#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        ContigGraphVertexAdaptor current(&vertices_[i]);

        if (current.contig_size() < min_length + kmer_size_ - 1
                && ((current.in_edges().size() <= 1 && current.out_edges().size() <= 1)
                        || current.in_edges().size() == 0 || current.out_edges().size() == 0)
           )
        {
            if (is_changed && current.coverage() > min_cover)
                continue;

            double mean = LocalCoverage(current, region_length);
            double threshold = min_cover;
            if (min_cover < mean * ratio)
                is_changed = true;
            else
                threshold = mean * ratio;

            if (current.coverage() < threshold)
            {
                is_changed = true;
                current.status().SetDeadFlag();
            }
        }
    }

    Refresh();
    //Trim(min_length);
    MergeSimplePaths();

    return is_changed;
}

double ContigGraph::LocalCoverage(ContigGraphVertexAdaptor current, int region_length)
{
    double num_count = 0;
    int num_kmer = 0;
    double x = LocalCoverageSingle(current, region_length, num_count, num_kmer);
    double y = LocalCoverageSingle(current.ReverseComplement(), region_length, num_count, num_kmer);

    if (num_kmer == 0)
        return 0;
    else
        return num_count / num_kmer;
}

double ContigGraph::LocalCoverageSingle(ContigGraphVertexAdaptor current, int region_length, double &total_count, int &total_kmer)
{
    map<int, int> visited;
    deque<ContigGraphVertexAdaptor> qu;
    qu.push_back(current);
    visited[current.id()] = 0;

    int index = 0;
    int num_added = 0;
    int num_count = 0;
    int num_kmer = 0;
    while (index < (int)qu.size())
    {
        current = qu[index++];

        if (num_added >= 4 * region_length)
            break;

        if (visited.size() > 32)
            break;

        if (visited[current.id()] >= region_length)
            continue;

        int dist = visited[current.id()];

        for (int x = 0; x < 4; ++x)
        {
            if (current.out_edges()[x])
            {
                ContigGraphVertexAdaptor next = GetNeighbor(current, x);
                if (visited.find(next.id()) == visited.end())
                {
                    visited[next.id()] = dist + next.num_kmer();
                    qu.push_back(next);

                    if ((int)next.num_kmer() + dist > region_length)
                    {
                        if (next.num_kmer() < region_length)
                        {
                            num_count += (int64_t)next.kmer_count() * (region_length - dist) / next.num_kmer();
                            num_kmer += region_length - dist;
                            num_added += region_length - dist;
                        }
                        else
                        {
                            Kmer begin = next.begin_kmer(kmer_size_);
                            if (in_kmer_count_table_.find(begin) == in_kmer_count_table_.end())
                            {
                                int in_kmer_count = 0;
                                for (int i = 0; i < region_length; ++i)
                                    in_kmer_count += next.get_count(i);
                                in_kmer_count_table_[begin] = in_kmer_count;
                            }

                            num_count += (int64_t)in_kmer_count_table_[begin] * (region_length - dist) / region_length;
                            num_kmer += region_length - dist;
                            num_added += region_length - dist;

//                            for (int i = 0; i < region_length - dist; ++i)
//                            {
//                                num_count += next.get_count(i);
//                                num_kmer++;
//                                num_added++;
//                            }
                        }
                    }
                    else
                    {
                        num_count += next.kmer_count();
                        num_kmer += next.num_kmer();
                        num_added += next.num_kmer();
                    }
                }
            }
        }
    }

    total_count += num_count;
    total_kmer += num_kmer;

    if (num_kmer == 0)
        return 0;
    else
        return num_count * 1.0 / num_kmer;
}

void ContigGraph::MergeSimplePaths()
{
    deque<Sequence> contigs;
    deque<ContigInfo> contig_infos;
    Assemble(contigs, contig_infos);
    Initialize(contigs, contig_infos);
}


int64_t ContigGraph::Assemble(deque<Sequence> &contigs, deque<ContigInfo> &contig_infos)
{
    contigs.clear();
    contig_infos.clear();

    omp_lock_t contig_lock;
    omp_init_lock(&contig_lock);

#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        if (vertices_[i].contig().size() == kmer_size_
                && vertices_[i].contig().IsPalindrome())
        {
            vertices_[i].status().Lock(omp_get_max_threads());

            Sequence contig = vertices_[i].contig();
            //ContigInfo contig_info(vertices_[i].kmer_count(), vertices_[i].in_edges(), vertices_[i].out_edges());
            ContigInfo contig_info;
            contig_info.set_kmer_count(vertices_[i].kmer_count());
            contig_info.in_edges() = vertices_[i].in_edges();
            contig_info.out_edges() = vertices_[i].out_edges();

            omp_set_lock(&contig_lock);
            contigs.push_back(contig);
            contig_infos.push_back(contig_info);
            omp_unset_lock(&contig_lock);
        }
    }

    //cout << "palindrome " << contigs.size() << endl;

#pragma omp parallel for schedule(static, 1)
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        if (!vertices_[i].status().Lock(omp_get_thread_num()))
            continue;

        ContigGraphPath path;
        path.Append(ContigGraphVertexAdaptor(&vertices_[i]), 0);

        Sequence contig;
        ContigInfo contig_info;
        for (int strand = 0; strand < 2; ++strand)
        {
            while (true)
            {
                ContigGraphVertexAdaptor current = path.back();
                ContigGraphVertexAdaptor next;

                if (!GetNextVertexAdaptor(current, next))
                    break;

                if (IsLoop(path, next))
                    break;

                if (!next.status().LockPreempt(omp_get_thread_num()))
                    goto FAIL;

                path.Append(next, -kmer_size_ + 1);
            }

            path.ReverseComplement();
        }

        path.Assemble(contig, contig_info);
        omp_set_lock(&contig_lock);
        contigs.push_back(contig);
        contig_infos.push_back(contig_info);
        omp_unset_lock(&contig_lock);
FAIL:
        ;
    }

    omp_destroy_lock(&contig_lock);
    ClearStatus();

    return contigs.size();
}

void ContigGraph::GetComponents(deque<deque<ContigGraphVertexAdaptor> > &components, deque<string> &component_strings)
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
                //for (connection_list_iterator p = connections()[current].begin(); p != connections()[current].end(); ++p)
                for (int x = 0; x < 4; ++x)
                {
                    if (current.out_edges()[x])
                    {
                        ContigGraphVertexAdaptor next = GetNeighbor(current, x);

                        if (strand == 0)
                        {
                            ss << current.id() << "_" << current.contig_size() << " " << next.id() << "_" << next.contig_size() << endl;

                            if (!next.status().IsUsed())
                                qu.push_back(next);
                        }
                        else
                        {
                            ss << next.id() << "_" << next.contig_size() << " " << current.id() << "_" << current.contig_size() << endl;

                            if (!next.status().IsUsed())
                                qu.push_back(next.ReverseComplement());
                        }

                        next.status().SetUsedFlag();
                    }
                }

                current.ReverseComplement();
            }
        }

        components.push_back(qu);
        component_strings.push_back(ss.str());
    }

    ClearStatus();
}

void ContigGraph::GetContigs(deque<Sequence> &contigs, deque<ContigInfo> &contig_infos)
{
    contigs.resize(vertices_.size());
    contig_infos.resize(vertices_.size());

#pragma omp parallel for 
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        contigs[i] = vertices_[i].contig();
        contig_infos[i] = vertices_[i].contig_info();
    }
}

void ContigGraph::BuildBeginKmerMap()
{
    begin_kmer_map_.clear();
    begin_kmer_map_.reserve(vertices_.size()*2);
#pragma omp parallel for
    for (int64_t i = 0; i < (int64_t)vertices_.size(); ++i)
    {
        for (int strand = 0; strand < 2; ++strand)
        {
            ContigGraphVertexAdaptor current(&vertices_[i], strand);
            Kmer kmer = current.begin_kmer(kmer_size_);

            Kmer key = kmer.unique_format();
            begin_kmer_map_[key] = i;
        }
    }
}


