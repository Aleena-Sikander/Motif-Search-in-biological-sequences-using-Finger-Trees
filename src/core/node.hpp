#pragma once
#include <vector>
#include <cassert>
#include "bio_types.hpp"

// A Node groups 2 or 3 elements together for storage in the spine.
// It caches their combined measure so we never recompute it.
template <typename T, typename M>
struct Node {
    std::vector<T> elems;       // exactly 2 or 3 elements
    M cached = M::zero();       // pre-computed combined measure

    Node() = default;

    explicit Node(std::vector<T> e) : elems(e), cached(M::zero()) {
        assert(elems.size() == 2 || elems.size() == 3);
        for (const T& x : elems)
            cached = cached + Measured<T, M>::measure(x);
    }

    const M& measure() const { return cached; }
};

// lets the tree treat a Node as a measurable element (used in the spine)
template <typename T, typename M>
struct Measured<Node<T, M>, M> {
    static M measure(const Node<T, M>& n) { return n.measure(); }
};
