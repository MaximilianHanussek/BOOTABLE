/**
 * @file bit_edges.h
 * @brief BitEdges Class.
 * @author Yu Peng (ypeng@cs.hku.hk)
 * @version 1.0.0
 * @date 2011-08-26
 */

#ifndef __GRAPH_BIT_EDGES_H_

#define __GRAPH_BIT_EDGES_H_

#include "basic/bit_operation.h"
#include "basic/atomic_integer.h"

class BitEdges
{
public:
    BitEdges() {}
    BitEdges(const BitEdges &bit_edges)
        : edges_(bit_edges.edges_) {}
    explicit BitEdges(uint8_t edges)
        : edges_(edges) {}

    ~BitEdges() {}

    const BitEdges &operator =(const BitEdges &bit_edges)
    { edges_ = bit_edges.edges_; return *this; }

    const BitEdges &operator =(uint8_t edges)
    { edges_ = edges; return *this; }

    operator uint8_t () const { return edges_; }

    void Add(int x) { edges_ |= uint8_t(1 << x); }
    void Remove(int x) { edges_ &= ~uint8_t(1 << x); }

    void swap(BitEdges &bit_edges) 
    { if (this != &bit_edges) edges_.swap(bit_edges.edges_); }

    bool operator [] (int index) const { return edges_ & (1 << index); }
    int size() const { return bit_operation::BitCount(edges_); }
    bool empty() const { return edges_ == 0; }

    void clear() { edges_ = 0; }

private:
    AtomicInteger<uint8_t> edges_;
};

namespace std
{
template <> inline void swap(BitEdges &x, BitEdges &y)
{ x.swap(y); }
}

#endif

