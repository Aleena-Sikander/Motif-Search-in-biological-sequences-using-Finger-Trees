/*
 * Sequence.h
 * ----------
 * A DNA sequence stored inside a Finger Tree.
 *
 * WHAT IS THIS FILE?
 * -------------------
 * Think of this as the "friendly wrapper" around the Finger Tree.
 * The Finger Tree is powerful but low-level (it deals with raw Measured
 * pointers and tree states).  Sequence hides all that complexity and
 * gives the rest of the program a clean, human-readable interface:
 *
 *   seq.insert(3, 'A')   — insert base 'A' at position 3
 *   seq.remove(5)        — delete the base at position 5
 *   seq.toString()       — get the full sequence as a plain string
 *
 * WHY NOT JUST USE std::string?
 * ------------------------------
 * std::string insert/erase are O(n) — they shift every character after
 * the edit point.  On a 1,000,000-base genome that is a million copies.
 * Sequence backed by a Finger Tree does the same edit in O(log n) steps.
 *
 * RELATIONSHIP TO OTHER FILES
 * ----------------------------
 *   FingerTree.h  ← Sequence.h ← MotifSearch.h
 *                              ← GCContent.h
 *                              ← Inversion.h
 *
 * Sequence is the central object that every algorithm operates on.
 */

#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "FingerTree.h"
#include <string>

class Sequence {
private:
    /*
     * tree
     * -----
     * The underlying Finger Tree that physically stores the DNA characters.
     * Each leaf of the tree is a Base object holding one character (A/T/G/C).
     *
     * WHY shared_ptr<FingerTree>: The Finger Tree is persistent/immutable —
     * every operation returns a NEW tree rather than modifying the old one.
     * Storing it in a shared_ptr lets us cheaply re-point to the new tree
     * after each edit without copying the data.
     */
    std::shared_ptr<FingerTree> tree;

public:
    /*
     * Sequence()
     * -----------
     * Creates an empty DNA sequence (no bases).
     * Useful as a starting point before loading data, or when building
     * a sequence programmatically base by base.
     */
    Sequence();

    /*
     * Sequence(const std::string& str)
     * ---------------------------------
     * Constructs a sequence directly from a plain string such as "ATGCGT".
     * Each character is validated against {A, T, G, C} and pushed into
     * the Finger Tree one by one using pushBack().
     *
     * Time: O(n) — one pushBack per character, each O(1) amortised.
     *
     * Used by: main(), test cases, loading sequences from FASTA files.
     */
    explicit Sequence(const std::string& str);

    /* -------------------------------------------------------
     *  GENE EDITING OPERATIONS
     *  (these correspond to Objective 2 in the project proposal)
     * ------------------------------------------------------- */

    /*
     * insert(index, value)
     * ---------------------
     * Inserts a DNA base at the given position, shifting everything after
     * it one position to the right.
     *
     * HOW IT WORKS UNDER THE HOOD:
     *   1. split(tree, index)  → leftPart + rightPart         O(log n)
     *   2. pushFront(rightPart, Base(value))                   O(1)
     *   3. concat(leftPart, newRight)                          O(log n)
     *   Total: O(log n)
     *
     * Analogy: snipping the necklace at bead #index, adding the new bead,
     * then tying the two halves back together.
     *
     * Example: "ATGC".insert(2, 'A') → "ATAGC"
     */
    void insert(int index, char value);

    /*
     * remove(index)
     * --------------
     * Deletes the DNA base at the given position.
     *
     * HOW IT WORKS:
     *   1. split(tree, index)       → leftPart + rightPart    O(log n)
     *   2. split(rightPart, 1)      → discarded + rest        O(log n)
     *   3. concat(leftPart, rest)                             O(log n)
     *   Total: O(log n)
     *
     * Example: "ATGC".remove(1) → "AGC"
     */
    void remove(int index);

    /*
     * update(index, value)
     * ---------------------
     * Replaces the base at `index` with a new character.
     * Equivalent to remove(index) + insert(index, value), but
     * can be done in a single pass for efficiency.
     *
     * Example: "ATGC".update(0, 'G') → "GTGC"
     */
    void update(int index, char value);

    /* -------------------------------------------------------
     *  QUERY OPERATIONS
     * ------------------------------------------------------- */

    /*
     * length()
     * ---------
     * Returns the total number of DNA bases in the sequence.
     * Reads the cached Measure from the tree root — O(1), no traversal.
     *
     * Used by: MotifSearch (loop bound), GCContent (range validation),
     *          Inversion (range validation).
     */
    int length() const;

    /*
     * get(index)
     * -----------
     * Returns the DNA character at position `index` (0-based).
     * Uses the tree's cached sizes to navigate without visiting every
     * node — O(log n).
     *
     * Used by: Sequence::toString(), MotifSearch sliding-window comparison.
     */
    char get(int index) const;

    /*
     * substring(left, right)
     * -----------------------
     * Extracts the subsequence from position `left` to `right` (exclusive)
     * as a plain std::string.
     *
     * HOW IT WORKS:
     *   1. split at `left`  → discard left part, keep right
     *   2. split at (right - left) → keep that many chars, discard rest
     *   3. convert the middle tree to string
     *   4. reconstruct the original tree (concat all three parts back)
     *   Total: O(log n) for splits + O(k) for string conversion where k
     *          is the length of the extracted substring.
     *
     * Used by: MotifSearch to extract each candidate window and compare
     *          it against the query motif.
     */
    std::string substring(int left, int right) const;

    /*
     * toString()
     * -----------
     * Converts the entire sequence to a plain std::string.
     * Traverses every leaf of the tree — O(n).
     *
     * Used by: printing results to the console, correctness testing,
     *          writing output files.
     */
    std::string toString() const;

    /*
     * getTree()
     * ----------
     * Returns a reference to the raw Finger Tree.
     * Used by GCContent and Inversion which need to call split/concat
     * directly on the tree for efficiency.
     *
     * This exposes an implementation detail, so use sparingly.
     * Prefer the higher-level methods above when possible.
     */
    std::shared_ptr<FingerTree>& getTree();
    const std::shared_ptr<FingerTree>& getTree() const;
};

#endif // SEQUENCE_H