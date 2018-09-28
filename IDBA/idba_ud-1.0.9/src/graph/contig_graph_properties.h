/**
 * @file contig_graph_properties.h
 * @brief ContigGraphProperties Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-25
 */

#ifndef __GRAPH_CONTIG_GRAPH_PROPERTIES_H_

#define __GRAPH_CONTIG_GRAPH_PROPERTIES_H_

#include "graph/contig_graph_vertex.h"
#include "container/managed_list.h"
#include "basic/pool.h"

#include <deque>
#include <cstddef>

template <class T>
class ContigGraphProperties
{
public:
    typedef T value_type;
    typedef size_t size_type;
    typedef std::ptrdiff_t difference_type;
    typedef ManagedList<value_type> list_type;
    typedef typename list_type::iterator list_iterator;
    typedef typename list_type::node_pool_type pool_type;

    explicit ContigGraphProperties(size_type size = 0)
    { resize(size); pool_ = new pool_type(); }

    ~ContigGraphProperties()
    { in_properties_.clear(); out_properties_.clear(); delete pool_; }

    list_type &operator [](const ContigGraphVertexAdaptor &x)
    { return (!x.is_reverse() ? out_properties_[x.id()] : in_properties_[x.id()]); }

    void swap(ContigGraphProperties &properties)
    { 
        if (this != &properties)
        {
            in_properties_.swap(properties.in_properties_);
            out_properties_.swap(properties.out_properties_);
            std::swap(pool_, properties.pool_);
            //pool_.swap(properties.pool_);

//            for (unsigned i = 0; i < in_properties_.size(); ++i)
//            {
//                in_properties_[i].set_pool(pool_);
//                out_properties_[i].set_pool(pool_);
//            }
        }
    }

    size_type size() const
    { return in_properties_.size(); }

    void resize(size_type size)
    {
        in_properties_.resize(size, list_type(*pool_));
        out_properties_.resize(size, list_type(*pool_));
    }

    void clear()
    {
        in_properties_.clear();
        out_properties_.clear();
        pool_->clear();
    }

private:
    ContigGraphProperties(const ContigGraphProperties &);
    const ContigGraphProperties &operator =(const ContigGraphProperties &);
    
    std::deque<list_type> in_properties_;
    std::deque<list_type> out_properties_;
    pool_type *pool_;
};

namespace std
{
template <typename T> inline void swap(ContigGraphProperties<T> &x, ContigGraphProperties<T> &y)
{ x.swap(y); }
}

#endif

