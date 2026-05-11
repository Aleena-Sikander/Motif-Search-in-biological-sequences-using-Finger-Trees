#pragma once

#include "core/finger_tree.hpp"
#include "core/bio_types.hpp"
#include "core/utils.hpp"

// ============================================================
// GCResult: result of a range GC content query
// ============================================================
struct GCResult {
    std::size_t length;
    std::size_t gc_count;
    double      gc_fraction;   // 0.0 if length == 0
};

// ============================================================
// gc_content(tree, start, end)
//   Compute GC content for the half-open range [start, end).
//
//   Uses two O(log n) splits to isolate the sub-tree, then
//   reads the cached Measure from the sub-tree root — O(1).
//   Total: O(log n).  No character-by-character scan.
// ============================================================
inline GCResult
gc_content(const FingerTree<char, Measure>& tree,
           std::size_t start, std::size_t end)
{
    const std::size_t n = tree.measure().length;
    if (start >= end || start >= n)
        return { 0, 0, 0.0 };
    end = std::min(end, n);

    auto [left, rest]  = split_at(tree, start);
    (void)left;
    auto [mid,  right] = split_at(rest, end - start);
    (void)right;

    const Measure& m = mid.measure();
    double frac = (m.length == 0) ? 0.0
                                  : static_cast<double>(m.gc_count) / m.length;
    return { m.length, m.gc_count, frac };
}
