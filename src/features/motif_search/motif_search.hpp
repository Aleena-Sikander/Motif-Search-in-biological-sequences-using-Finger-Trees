#pragma once
#include <string>
#include <vector>
#include "core/finger_tree.hpp"
#include "core/bio_types.hpp"
#include "core/utils.hpp"

// ---------------------------------------------------------------
// findMotif: search for a pattern inside a DNA finger tree
//
// Uses the KMP (Knuth-Morris-Pratt) algorithm, which runs in
// O(n + k) time where n = sequence length, k = motif length.
//
// The key point: we use the tree's for_each traversal to stream
// characters one by one through the KMP state machine. We never
// copy the whole sequence into a string first.
//
// Returns a list of all starting positions (0-based) where the
// motif was found.
// ---------------------------------------------------------------
inline std::vector<size_t> findMotif(const DNATree& tree, const std::string& motif) {
    if (motif.empty()) return {};

    size_t k = motif.size();

    // --- build the KMP failure table ---
    // fail[i] = length of the longest proper prefix of motif[0..i]
    // that is also a suffix. This lets us skip re-checking characters
    // we already matched when a mismatch happens.
    std::vector<int> fail(k, 0);
    for (size_t i = 1; i < k; i++) {
        int j = fail[i - 1];
        while (j > 0 && motif[i] != motif[j]) j = fail[j - 1];
        if (motif[i] == motif[j]) j++;
        fail[i] = j;
    }

    std::vector<size_t> hits;
    int matched = 0;   // how many characters of the motif we have matched so far
    size_t pos  = 0;   // current position in the sequence

    // stream characters from the tree one by one using for_each
    tree.for_each([&](char c) {
        // if mismatch, fall back using the failure table (don't restart from 0)
        while (matched > 0 && c != motif[matched])
            matched = fail[matched - 1];

        if (c == motif[matched]) matched++;

        if ((size_t)matched == k) {
            hits.push_back(pos - k + 1);  // found a complete match
            matched = fail[matched - 1];   // look for overlapping matches
        }
        pos++;
    });

    return hits;
}
