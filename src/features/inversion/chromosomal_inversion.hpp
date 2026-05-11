#pragma once

#include <algorithm>
#include <string>
#include "core/finger_tree.hpp"
#include "core/bio_types.hpp"
#include "core/utils.hpp"

// ============================================================
// chromosomal_inversion(tree, start, end)
//   Reverse the subsequence in the half-open range [start, end).
//
//   Algorithm (all O(log n)):
//     1. split_at(start)  -> left | rest
//     2. split_at(end-start) on rest -> segment | right
//     3. flatten segment to string, std::reverse it
//     4. rebuild segment tree, concat left + rev_seg + right
//
//   The string reverse step is O(end - start), which is fine
//   for genomic-scale inversions (typically << full genome).
// ============================================================
inline FingerTree<char, Measure>
chromosomal_inversion(const FingerTree<char, Measure>& tree,
                      std::size_t start, std::size_t end)
{
    const std::size_t n = tree.measure().length;
    if (start >= end || start >= n) return tree;
    end = std::min(end, n);

    auto [left,  rest]    = split_at(tree, start);
    auto [segment, right] = split_at(rest, end - start);

    std::string seg = to_sequence(segment);
    std::reverse(seg.begin(), seg.end());

    FingerTree<char, Measure> rev_seg;
    for (char c : seg) rev_seg = rev_seg.push_back(c);

    return FingerTree<char, Measure>::concat(
               FingerTree<char, Measure>::concat(left, rev_seg),
               right);
}
