/*
 * FingerTree.h
 * ------------
 * The heart of the entire project.
 *
 * WHAT IS A FINGER TREE? (plain English)
 * ----------------------------------------
 * Imagine a very long necklace of DNA beads (A, T, G, C).
 * A regular array is like laying the necklace flat on a table —
 * finding the middle bead means counting from one end every time.
 *
 * A Finger Tree is like holding the necklace with both hands ("fingers"):
 *   - Your LEFT hand grips a few beads at the front  (prefix)
 *   - Your RIGHT hand grips a few beads at the back   (suffix)
 *   - Everything in the middle is bundled neatly into groups (nodes)
 *
 * Because you're already touching both ends, adding or removing a bead
 * from either end is nearly instant (O(1) amortised).
 * Splitting the necklace in half, or joining two necklaces, only takes
 * O(log n) steps — much faster than rebuilding from scratch.
 *
 * WHY DO WE NEED THIS FOR THE PROJECT?
 * --------------------------------------
 * DNA sequences can be millions of characters long. The project needs to:
 *   1. Search for a short pattern (motif) inside that long sequence
 *   2. Insert/delete/update individual bases (gene editing)
 *   3. Cut out a chunk and flip it (chromosomal inversion)
 *   4. Measure the G/C content of any range
 *
 * A plain std::string or std::vector makes cuts and joins O(n) — slow.
 * The Finger Tree makes them O(log n) — fast even on huge genomes.
 *
 * NOTE ON YOUR IMPLEMENTATION vs YOUR FRIEND'S
 * ----------------------------------------------
 * Your friend's implementation (the single-file main.cpp) is a TRUE
 * Finger Tree: it has prefix/suffix fingers, internal 2-3 nodes, cached
 * size measurements (the Monoid), O(log n) split, and O(log n) concat
 * via app3().  It is architecturally correct.
 *
 * YOUR implementation in FingerTrees.cpp is actually a LINKED LIST
 * disguised as a Finger Tree.  The Node struct has only `left` and
 * `right` pointers (making it a doubly-linked list), and the critical
 * operations split() and concat() are left as placeholders that do
 * nothing.  The project cannot currently do motif search, GC content,
 * or inversion correctly because those all depend on split/concat.
 *
 * RECOMMENDATION:
 *   Adopt your friend's structural design (Empty / Single / Deep with
 *   real prefix+suffix+middle) but keep YOUR modular file layout
 *   (Sequence.h, MotifSearch.h, GCContent.h, Inversion.h, Utils.h).
 *   That gives you the best of both worlds.
 */

#ifndef FINGERTREE_H
#define FINGERTREE_H

#include <memory>
#include <vector>
#include <string>
#include <functional>

/*
 * Measure (the Monoid)
 * ---------------------
 * Every node in the tree caches a small summary of everything beneath it.
 * Right now the summary is just the total number of characters (size).
 *
 * WHY: When searching for a motif at position i, we can look at a node's
 * cached size and immediately know "the match can't start here — this
 * whole subtree has only 3 chars but the motif is 5 chars long".
 * We skip entire branches without touching individual characters.
 * That is what makes O(log n) search possible.
 *
 * You can extend Measure later to cache prefix/suffix substrings so
 * motif candidates that straddle node boundaries are also caught quickly.
 */
struct Measure {
    size_t size;   // total number of DNA characters in this subtree

    // Combining two measures: just add their sizes.
    // This "addition" is the monoid operation.
    Measure operator+(const Measure& other) const {
        return { size + other.size };
    }
};

/*
 * Measured (interface)
 * ---------------------
 * Anything that can live inside the tree (a single DNA base, or an
 * internal node grouping several bases) must be able to:
 *   1. Report its measure (size)
 *   2. Write its characters into an output string
 */
struct Measured {
    virtual ~Measured() = default;

    // Returns the cached size of this element.
    virtual Measure getMeasure() const = 0;

    // Appends all DNA characters under this element to `out`.
    virtual void gatherBases(std::string& out) const = 0;
};

/*
 * Base
 * -----
 * A single DNA character: 'A', 'T', 'G', or 'C'.
 * This is a LEAF of the finger tree — the actual data being stored.
 *
 * WHY: The tree stores characters as leaf objects so every node can hold
 * a uniform Measured* pointer regardless of whether it points to a
 * single base or a grouped node.
 */
struct Base : public Measured {
    char type;   // one of: A, T, G, C (or amino acid letters for proteins)

    explicit Base(char t);

    // A single base always has measure {1} — it contributes exactly one
    // character to the total size.
    Measure getMeasure() const override;

    // Appends this one character to the output string.
    void gatherBases(std::string& out) const override;
};

/*
 * Node (2-3 internal node)
 * -------------------------
 * Groups 2 or 3 child Measured objects together into a single unit.
 * Internal nodes are what make the tree "balanced" — the middle spine
 * of a Deep tree holds Nodes of Nodes of Nodes, each level grouping
 * elements by 2 or 3, keeping height at O(log n).
 *
 * WHY: Without grouping, adding one element to the prefix when it
 * already has 4 items would just shove an element into the middle tree.
 * Grouping 3 of the prefix items into a Node before pushing keeps the
 * prefix/suffix sizes bounded (1-4 items each), preserving amortised O(1)
 * prepend/append.
 */
struct Node : public Measured {
    std::vector<std::shared_ptr<Measured>> children;  // 2 or 3 children
    Measure m;   // cached combined measure of all children

    explicit Node(std::vector<std::shared_ptr<Measured>> c);

    // Returns the pre-cached measure (no recomputation needed).
    Measure getMeasure() const override;

    // Recursively gather all base characters from all children.
    void gatherBases(std::string& out) const override;
};

/*
 * FingerTree (abstract base)
 * ---------------------------
 * The three possible states of a Finger Tree:
 *   Empty  — no elements at all
 *   Single — exactly one element
 *   Deep   — two or more elements (has a prefix, middle, and suffix)
 *
 * All three inherit from FingerTree so they can be stored in the same
 * shared_ptr<FingerTree> without the caller knowing which state it is.
 */
class FingerTree {
public:
    virtual ~FingerTree() = default;

    // Returns the total number of DNA characters stored in this tree.
    virtual Measure getMeasure() const = 0;

    // Returns true only for the Empty state.
    virtual bool isEmpty() const = 0;
};

/*
 * Empty
 * ------
 * Represents an empty sequence — no DNA bases at all.
 * Acts as the identity element for concat (concat(t, Empty) == t).
 */
class Empty : public FingerTree {
public:
    Measure getMeasure() const override;   // always {0}
    bool isEmpty() const override;         // always true
};

/*
 * Single
 * -------
 * Represents a sequence with exactly one element.
 * Sits between Empty and Deep to avoid needing empty prefix/suffix.
 * Promotes to Deep when a second element is pushed in.
 */
class Single : public FingerTree {
public:
    std::shared_ptr<Measured> val;   // the one element

    explicit Single(std::shared_ptr<Measured> v);

    Measure getMeasure() const override;   // = val->getMeasure()
    bool isEmpty() const override;         // always false
};

/*
 * Deep
 * -----
 * The normal working state once a sequence has >= 2 elements.
 *
 *   prefix  : 1-4 elements accessible from the FRONT  (the left "finger")
 *   middle  : a recursive FingerTree<Node> holding the bulk of the data
 *   suffix  : 1-4 elements accessible from the BACK   (the right "finger")
 *
 * WHY TWO FINGERS: Biological operations frequently touch both ends.
 * For example, reading a sequence left-to-right (motif scan) accesses
 * the prefix constantly; a reverse-complement operation touches the
 * suffix.  Having O(1) access to both ends avoids unnecessary traversal.
 *
 * The cached measure `m` = prefix sizes + middle size + suffix sizes,
 * computed once at construction.  No traversal needed to know total length.
 */
class Deep : public FingerTree {
public:
    std::vector<std::shared_ptr<Measured>> prefix;   // left finger  (1-4 items)
    std::shared_ptr<FingerTree>            middle;   // recursive spine
    std::vector<std::shared_ptr<Measured>> suffix;   // right finger (1-4 items)
    Measure m;   // cached total measure

    Deep(std::vector<std::shared_ptr<Measured>> p,
         std::shared_ptr<FingerTree>            mid,
         std::vector<std::shared_ptr<Measured>> s);

    Measure getMeasure() const override;   // returns cached m
    bool isEmpty() const override;         // always false
};

/* =========================================================
 *  CORE TREE OPERATIONS
 * =========================================================
 *
 * All operations return a NEW tree; the original is unchanged.
 * This is the "functional / persistent" style — safe for concurrent
 * reads and easy to reason about correctness.
 */

/*
 * pushFront
 * ----------
 * Adds one element to the FRONT of the sequence in O(1) amortised time.
 *
 * Analogy: threading one new bead onto the left end of the necklace.
 * When the left hand (prefix) already holds 4 beads, it bundles 3 of
 * them into a Node and pushes that Node down into the middle tree before
 * accepting the new bead.  This keeps the prefix size bounded at 1-4.
 *
 * Used by: Sequence::insert(0, c), motif search (building a window),
 *          chromosomal inversion (rebuilding a reversed segment).
 */
std::shared_ptr<FingerTree> pushFront(std::shared_ptr<FingerTree> t,
                                      std::shared_ptr<Measured>   v);

/*
 * pushBack
 * ---------
 * Adds one element to the BACK of the sequence in O(1) amortised time.
 * Mirror image of pushFront.
 *
 * Used by: Sequence construction from a string, appending new bases,
 *          concatenating sequences (concat calls pushBack internally).
 */
std::shared_ptr<FingerTree> pushBack(std::shared_ptr<FingerTree> t,
                                     std::shared_ptr<Measured>   v);

/*
 * toNodes  (helper for concat)
 * -----------------------------
 * Takes a flat list of elements and groups them into 2-3 Nodes.
 * This is needed when the middle of two trees being concatenated
 * contains leftover suffix+buffer+prefix elements that must be
 * re-packed into proper nodes for the new middle spine.
 *
 * Rule: prefer groups of 3; use a group of 2 only when needed to
 * avoid a remainder of 1 (Nodes must have at least 2 children).
 */
std::vector<std::shared_ptr<Measured>> toNodes(
    std::vector<std::shared_ptr<Measured>> b);

/*
 * app3  (the real concatenation engine)
 * ---------------------------------------
 * Merges two trees with a list of "middle" elements `ts` in between.
 * Called recursively, carrying the suffix of t1 and prefix of t2
 * as it descends through the middle spines of both trees.
 *
 * Public-facing wrapper: concat(t1, t2) = app3(t1, {}, t2)
 *
 * Time: O(log n) — only the spine depth (O(log n)) levels of recursion.
 *
 * WHY IT MATTERS: Chromosomal inversion = split + reverse + concat.
 * GC content of a range = split + measure + concat back.
 * Without O(log n) concat this would be O(n) — too slow for large genomes.
 */
std::shared_ptr<FingerTree> app3(std::shared_ptr<FingerTree>             t1,
                                 std::vector<std::shared_ptr<Measured>>  ts,
                                 std::shared_ptr<FingerTree>             t2);

/*
 * concat
 * -------
 * Joins two sequences end-to-end in O(log n) time.
 *
 * Example: concat("ATGC", "GCTA") → "ATGCGCTA"
 *
 * Used by: gene editing (insert = split + pushFront + concat),
 *          inversion (split + reverse + concat),
 *          any operation that temporarily dismantles and re-joins the tree.
 */
std::shared_ptr<FingerTree> concat(std::shared_ptr<FingerTree> t1,
                                   std::shared_ptr<FingerTree> t2);

/*
 * split
 * ------
 * Divides the sequence at position `i` into two trees:
 *   left  — characters [0 .. i-1]
 *   right — characters [i .. n-1]
 *
 * Time: O(log n) — uses cached measures to navigate without visiting
 * every element.
 *
 * Used by:
 *   - Sequence::insert(i, c)  → split at i, pushFront c onto right, concat
 *   - Sequence::remove(i)     → split at i, split right at 1, concat
 *   - GCContent::calculate()  → split at left, split at right, measure middle
 *   - Inversion::invert()     → split at left, split at right, reverse middle, concat
 *   - MotifSearch::findMotif()→ extract a window of length |motif| to compare
 */
std::pair<std::shared_ptr<FingerTree>, std::shared_ptr<FingerTree>>
split(std::shared_ptr<FingerTree> t, size_t i);

/*
 * getSequence
 * ------------
 * Traverses the entire tree in order and appends every DNA character
 * to the output string `out`.
 *
 * Time: O(n) — must visit all n leaves.
 *
 * Used by: Sequence::toString(), printing results, building substrings
 *          during motif search window extraction.
 */
void getSequence(std::shared_ptr<FingerTree> t, std::string& out);

#endif // FINGERTREE_H