#pragma once
#include <string>
#include "finger_tree.hpp"
#include "bio_types.hpp"

// shorthand so we don't write FingerTree<char, Measure> everywhere
using DNATree = FingerTree<char, Measure>;

// build a finger tree from a plain string, one character at a time
inline DNATree build_tree(const std::string& seq) {
    DNATree t;
    for (char c : seq) t = t.push_back(c);
    return t;
}

// flatten a finger tree back into a plain string
inline std::string to_string(const DNATree& tree) {
    std::string result;
    result.reserve(tree.measure().length);
    tree.for_each([&](char c){ result += c; });
    return result;
}

// alias used by chromosomal_inversion.hpp
inline std::string to_sequence(const DNATree& tree) { return to_string(tree); }

// ---------------------------------------------------------------
// split_at: cut the tree into two halves at a given position
//   left  = characters [0, pos)
//   right = characters [pos, end)
//
// Uses the tree's O(log n) split operation.
// The predicate fires when we have seen >= pos characters,
// which puts the character at index pos-1 in the pivot slot.
// We move the pivot into the left side to get exactly pos chars.
// ---------------------------------------------------------------
struct SplitResult { DNATree left; DNATree right; };

inline SplitResult split_at(const DNATree& tree, size_t pos) {
    size_t len = tree.measure().length;
    if (pos == 0)   return { DNATree{}, tree };
    if (pos >= len) return { tree, DNATree{} };

    auto [left, pivot, right] = tree.split([pos](const Measure& m){
        return m.length >= pos;
    });

    // pivot is the character exactly at index pos-1; put it in the left half
    return { left.push_back(pivot), right };
}

// return the character at 0-based index i  O(log n)
inline char lookup_at(const DNATree& tree, size_t i) {
    size_t len = tree.measure().length;
    if (i >= len) throw std::out_of_range("lookup_at: index out of range");
    auto [left, right] = split_at(tree, i + 1);
    (void)right;
    auto vr = left.view_right();
    return vr.last;
}
