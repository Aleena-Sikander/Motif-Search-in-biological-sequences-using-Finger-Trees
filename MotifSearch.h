/*
 * MotifSearch.h
 * -------------
 * Locates all occurrences of a short pattern (motif) inside a DNA sequence.
 *
 * WHAT IS A MOTIF? (plain English)
 * ---------------------------------
 * A motif is a short, biologically meaningful pattern that appears
 * repeatedly in DNA.  For example:
 *
 *   TATATA — a common promoter motif that signals "gene starts here"
 *   GAATTC — the recognition site for the EcoRI restriction enzyme
 *   ATG    — the "start codon" that tells a ribosome where to begin
 *             reading a gene
 *
 * Finding every position where a motif occurs in a genome lets
 * biologists identify regulatory regions, binding sites, and
 * conserved functional elements — all without needing to be a biologist
 * yourself.  From the data structure perspective it is just:
 *
 *   "Find all positions where this short string appears in this long string."
 *
 * WHY NOT JUST USE std::string::find?
 * -------------------------------------
 * std::string::find is O(n * m) in the worst case (n = sequence length,
 * m = motif length).  For a 3-million-base genome and a 10-base motif
 * that is 30 million character comparisons.
 *
 * Using the Finger Tree we can:
 *   1. Use cached node sizes to SKIP subtrees that are too short to
 *      contain the motif — avoiding many comparisons entirely.
 *   2. Extract candidate windows in O(log n) each rather than O(n).
 *
 * The project proposal also requires comparing against a naive approach
 * and a hashing approach (e.g., Rabin-Karp).  This class implements the
 * Finger Tree-based approach.
 *
 * ALGORITHM OVERVIEW
 * -------------------
 * Sliding window, tree-accelerated:
 *
 *   for position i from 0 to (sequenceLength - motifLength):
 *       candidate = seq.substring(i, i + motifLength)   // O(log n + m)
 *       if candidate == motif:
 *           record i as a match
 *
 * An optimised version traverses the tree directly, using the cached
 * size at each node to skip branches where no match can start.
 *
 * Total worst-case: O(n * log n) — better than naive O(n * m) when the
 * tree is balanced and m is large.
 */

#ifndef MOTIFSEARCH_H
#define MOTIFSEARCH_H

#include "Sequence.h"
#include <vector>
#include <string>

class MotifSearch {
public:
    /*
     * findMotif(seq, motif)
     * ----------------------
     * Searches for every occurrence of `motif` inside `seq` and returns
     * a list of ALL starting positions (0-based) where it was found.
     *
     * Parameters:
     *   seq   — the full DNA sequence stored in a Finger Tree
     *   motif — the short pattern to search for (e.g., "GAATTC")
     *
     * Returns:
     *   A vector of integers, each being a 0-based index into `seq`
     *   where `motif` begins.  Returns an empty vector if not found.
     *
     * Example:
     *   seq   = "ATGAATTCGAATTCTT"
     *   motif = "GAATTC"
     *   result → {2, 8}   (motif appears starting at positions 2 and 8)
     *
     * HOW IT USES THE FINGER TREE:
     *   At each candidate position i, it calls seq.substring(i, i+m)
     *   which internally uses split() + getSequence() + concat() to
     *   extract exactly m characters in O(log n) rather than copying
     *   the whole string.  The cached size measures let it skip nodes
     *   where no complete motif could possibly fit.
     *
     * Time complexity: O(n * log n) using the Finger Tree
     *                  O(n * m)     naive / brute-force baseline
     *                  O(n + m)     Rabin-Karp hashing baseline
     *
     * IMPORTANT: motif matching is case-insensitive — both seq and motif
     * are uppercased before comparison (via Utils::toUpper).
     */
    static std::vector<int> findMotif(Sequence& seq, const std::string& motif);

    /*
     * findMotifNaive(seq, motif)
     * ---------------------------
     * A brute-force O(n * m) baseline that converts the sequence to a
     * plain std::string first and then slides a window character by
     * character.
     *
     * WHY INCLUDE THIS: The project proposal requires a performance
     * comparison between the Finger Tree approach and the naive approach.
     * This method provides the naive baseline with the same interface so
     * timing experiments are fair (same input, same output format).
     *
     * Returns the same result as findMotif() — use for correctness
     * verification and benchmarking only.
     */
    static std::vector<int> findMotifNaive(Sequence& seq, const std::string& motif);

    /*
     * findMotifHashing(seq, motif)
     * -----------------------------
     * Rabin-Karp rolling-hash approach: O(n + m) average case.
     *
     * Computes a hash of the motif, then rolls a hash window across
     * the sequence.  Only checks character equality when hashes match.
     *
     * WHY INCLUDE THIS: The project proposal lists hashing as the
     * third comparison approach.  Rabin-Karp is the standard choice
     * because its rolling hash avoids recomputing from scratch each step.
     *
     * Returns the same result as findMotif().
     */
    static std::vector<int> findMotifHashing(Sequence& seq, const std::string& motif);
};

#endif // MOTIFSEARCH_H
