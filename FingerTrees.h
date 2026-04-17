//generic fingertree stucture
/*
 * FingerTree.h
 * ------------
 * WHY DO WE NEED THIS FILE?
 *   Think of a Finger Tree like a very smart bookshelf. A regular array is like a
 *   bookshelf where you can only grab books from the middle (slow!). A Finger Tree
 *   lets you grab from BOTH ends instantly, and split or join shelves in O(log n)
 *   time. For a DNA sequence of millions of characters, this matters enormously.
 *
 *   This is the HEART of the entire project. Every other file (Sequence, MotifSearch,
 *   GCContent, Inversion) depends on what is defined here.
 *
 * REAL-WORLD ANALOGY:
 *   Imagine a train where you always need to add/remove carriages from either end,
 *   and occasionally split the train in half or merge two trains. A linked list is
 *   slow in the middle. A Finger Tree puts "fingers" (fast access points) on both
 *   ends so those operations are cheap.
 *
 * KEY OPERATIONS AND THEIR TIME COMPLEXITIES:
 *   - pushFront / pushBack   → O(1) amortised  (add a DNA base to either end)
 *   - split(i)               → O(log n)         (cut the sequence at position i)
 *   - concat(t1, t2)         → O(log n)         (glue two sequences together)
 *   - lookup(i)              → O(log n)         (read base at position i)
 *
 * STRUCTURE OVERVIEW:
 *   A FingerTree is one of three shapes:
 *     Empty  — the tree has no elements.
 *     Single — exactly one element.
 *     Deep   — has a prefix (up to 4 elements), a recursive middle spine,
 *              and a suffix (up to 4 elements). The middle holds Node objects
 *              that group children in groups of 2 or 3.
 *
 *   Cached "measurements" (just the size count for this project) are stored at
 *   every internal node so we can navigate without scanning the whole tree.
 */

#ifndef FINGERTREE_H
#define FINGERTREE_H

#include <memory>
#include <vector>
#include <string>

// ---------------------------------------------------------------------------
// Measure (Monoid)
// ---------------------------------------------------------------------------
/*
 * WHY A Measure?
 *   Every node in the tree caches a small piece of summary information called
 *   a "measure". Here the measure is simply the count of DNA bases under that
 *   node. Because it is a monoid (it has a zero = {0} and can be combined with
 *   +), the tree can propagate it automatically.
 *
 *   When you later want to split at position 50, the tree checks the cached
 *   size of each subtree and decides which branch to descend into — exactly
 *   like binary search. Without the cached measure you would have to count
 *   every leaf manually, turning an O(log n) split into O(n).
 */
struct Measure {
    size_t size;   // number of DNA bases in this subtree

    /* Combine two measures (used when building internal nodes) */
    Measure operator+(const Measure& other) const {
        return { size + other.size };
    }
};

// ---------------------------------------------------------------------------
// Measured (interface)
// ---------------------------------------------------------------------------
/*
 * WHY AN INTERFACE?
 *   Both individual DNA bases (Base) and internal grouping nodes (Node) need
 *   to report a measure. Making them share this interface lets the tree treat
 *   them uniformly without caring which one it currently holds.
 */
struct Measured {
    virtual ~Measured() = default;

    /* Return the cached measure for this element or subtree. */
    virtual Measure getMeasure() const = 0;

    /* Append every leaf character to 'out'. Used to reconstruct a string. */
    virtual void gatherBases(std::string& out) const = 0;
};

// ---------------------------------------------------------------------------
// Base (a single DNA character: A, T, G, or C)
// ---------------------------------------------------------------------------
/*
 * WHY A SEPARATE CLASS FOR A SINGLE CHARACTER?
 *   In a plain std::string a character is just a byte. Here we wrap it so it
 *   can participate in the Measured interface — it reports a size of 1 and
 *   knows how to append itself to an output string. These are the LEAVES of
 *   the tree.
 */
struct Base : public Measured {
    char nucleotide;   // one of 'A', 'T', 'G', 'C'

    explicit Base(char c) : nucleotide(c) {}

    /* A single base always has size 1. */
    Measure getMeasure() const override { return { 1 }; }

    /* Append this nucleotide to the accumulator string. */
    void gatherBases(std::string& out) const override { out += nucleotide; }
};

// ---------------------------------------------------------------------------
// Node (internal 2-3 grouping node)
// ---------------------------------------------------------------------------
/*
 * WHY GROUP ELEMENTS INTO NODES?
 *   The Finger Tree stores bases individually only at the outer "digits"
 *   (prefix/suffix). Elements that get pushed into the inner spine are
 *   first bundled into Node objects containing 2 or 3 children. This
 *   bundling keeps the tree height logarithmic even for very long sequences.
 *
 *   Think of it as packing books into boxes before putting them on the
 *   high shelf — you handle fewer items at each level.
 */
struct Node : public Measured {
    std::vector<std::shared_ptr<Measured>> children;  // 2 or 3 children
    Measure cached;   // pre-computed total size of all children

    explicit Node(std::vector<std::shared_ptr<Measured>> c);

    /* Return the pre-computed measure (O(1) lookup). */
    Measure getMeasure() const override { return cached; }

    /* Recurse into each child and collect their characters. */
    void gatherBases(std::string& out) const override;
};

// ---------------------------------------------------------------------------
// FingerTree variants
// ---------------------------------------------------------------------------

/* Abstract base — lets us store any of the three shapes behind a pointer. */
struct FingerTree {
    virtual ~FingerTree() = default;
    virtual Measure getMeasure() const = 0;
    virtual bool isEmpty() const = 0;
};

/*
 * Empty — the base case; a sequence with no elements.
 * Measure is {0}.
 */
struct Empty : public FingerTree {
    Measure getMeasure() const override { return { 0 }; }
    bool isEmpty() const override { return true; }
};

/*
 * Single — exactly one element (a Base or a Node from a deeper level).
 * Avoids allocating prefix/suffix vectors when the tree contains just one thing.
 */
struct Single : public FingerTree {
    std::shared_ptr<Measured> value;

    explicit Single(std::shared_ptr<Measured> v) : value(v) {}
    Measure getMeasure() const override { return value->getMeasure(); }
    bool isEmpty() const override { return false; }
};

/*
 * Deep — the general case.
 *
 *   prefix  — 1-4 elements at the left "finger" (fast prepend/pop-front)
 *   middle  — a recursive FingerTree of Nodes (deeper elements)
 *   suffix  — 1-4 elements at the right "finger" (fast append/pop-back)
 *   cached  — total measure of prefix + middle + suffix
 *
 * WHY PREFIX AND SUFFIX HAVE AT MOST 4 ELEMENTS?
 *   When either digit grows beyond 4, the tree "overflows" three of those
 *   elements into the middle spine as a bundled Node. This keeps the digits
 *   small (O(1) access) and the height logarithmic.
 */
struct Deep : public FingerTree {
    std::vector<std::shared_ptr<Measured>> prefix;
    std::shared_ptr<FingerTree>            middle;
    std::vector<std::shared_ptr<Measured>> suffix;
    Measure cached;

    Deep(std::vector<std::shared_ptr<Measured>> p,
         std::shared_ptr<FingerTree>            mid,
         std::vector<std::shared_ptr<Measured>> s);

    Measure getMeasure() const override { return cached; }
    bool isEmpty() const override { return false; }
};

// ---------------------------------------------------------------------------
// Core tree operations (free functions)
// ---------------------------------------------------------------------------

/*
 * pushFront(tree, element)
 *   Add one element to the left end of the tree.
 *   Cost: O(1) amortised.
 *
 *   Used in: building a sequence left-to-right, or after a split when we
 *   need to re-attach a head element.
 */
std::shared_ptr<FingerTree> pushFront(std::shared_ptr<FingerTree> t,
                                      std::shared_ptr<Measured>   v);

/*
 * pushBack(tree, element)
 *   Add one element to the right end of the tree.
 *   Cost: O(1) amortised.
 *
 *   Used in: building a sequence from a DNA string character by character,
 *   and in the concatenation helper (app3).
 */
std::shared_ptr<FingerTree> pushBack(std::shared_ptr<FingerTree> t,
                                     std::shared_ptr<Measured>   v);

/*
 * concat(left, right)
 *   Merge two trees into one, preserving order.
 *   Cost: O(log n).
 *
 *   Used in: Inversion (split → reverse segment → concat back),
 *   GCContent (split to isolate range → measure → concat to restore).
 *
 *   Internally delegates to app3() which threads the boundary digits
 *   through the recursive spine as grouped Nodes.
 */
std::shared_ptr<FingerTree> concat(std::shared_ptr<FingerTree> left,
                                   std::shared_ptr<FingerTree> right);

/*
 * getSequence(tree, out)
 *   Perform an in-order traversal and append every leaf character to 'out'.
 *   Cost: O(n).
 *
 *   Used whenever you need to convert the tree back to a plain string,
 *   for example to display results or run a naive substring compare during
 *   motif search verification.
 */
void getSequence(std::shared_ptr<FingerTree> t, std::string& out);

// ---------------------------------------------------------------------------
// Internal helper — not part of public API
// ---------------------------------------------------------------------------

/*
 * app3(left, middle_elements, right)
 *   The recursive engine behind concat(). It carries a list of "in-between"
 *   elements that need to be absorbed into the spine as Nodes. End-users
 *   should call concat() instead.
 */
std::shared_ptr<FingerTree> app3(std::shared_ptr<FingerTree>              t1,
                                 std::vector<std::shared_ptr<Measured>>   ts,
                                 std::shared_ptr<FingerTree>              t2);

/*
 * toNodes(elements)
 *   Bundle a flat vector of elements into 2-child or 3-child Nodes.
 *   Called by app3 to re-group the boundary digits before pushing them
 *   into the middle spine.
 */
std::vector<std::shared_ptr<Measured>> toNodes(
        std::vector<std::shared_ptr<Measured>> elements);

#endif // FINGERTREE_H