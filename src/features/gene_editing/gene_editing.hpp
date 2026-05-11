#pragma once
#include <string>
#include <algorithm>
#include "core/finger_tree.hpp"
#include "core/bio_types.hpp"
#include "core/utils.hpp"

// ---------------------------------------------------------------
// All three gene editing operations follow the same pattern:
//
//   1. split the tree at the target position(s)  — O(log n)
//   2. do something to the pieces
//   3. concat the pieces back together            — O(log n)
//
// The finger tree makes this efficient because split and concat
// are both O(log n), unlike a string where inserting in the
// middle means copying everything.
// ---------------------------------------------------------------

// insertAt: insert a new sequence starting at position pos
//
//   before:  [  left part  |  right part  ]
//   after:   [  left part  | new_seq | right part  ]
//
inline DNATree insertAt(const DNATree& tree, size_t pos, const std::string& new_seq) {
    auto [left, right] = split_at(tree, pos);

    // build a small tree from the new sequence
    DNATree inserted;
    for (char c : new_seq) inserted = inserted.push_back(c);

    return DNATree::concat(DNATree::concat(left, inserted), right);
}

// deleteRange: remove characters in [start, end)
//
//   before:  [  left  |  deleted  |  right  ]
//   after:   [  left  |  right  ]
//
inline DNATree deleteRange(const DNATree& tree, size_t start, size_t end) {
    size_t len = tree.measure().length;
    if (start >= end || start >= len) return tree;  // nothing to delete
    end = std::min(end, len);

    auto [left,  rest]    = split_at(tree, start);
    auto [deleted, right] = split_at(rest, end - start);
    (void)deleted;  // we just discard the middle part

    return DNATree::concat(left, right);
}

// updateRange: replace characters starting at pos with replacement
//
// This is just a delete followed by an insert.
//
inline DNATree updateRange(const DNATree& tree, size_t pos, const std::string& replacement) {
    DNATree after_delete = deleteRange(tree, pos, pos + replacement.size());
    return insertAt(after_delete, pos, replacement);
}
