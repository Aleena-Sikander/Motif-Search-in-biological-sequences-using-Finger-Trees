/*
 * Inversion.h
 * ------------
 * Performs a chromosomal inversion: cuts out a segment of DNA,
 * reverses it, and puts it back.
 *
 * WHAT IS A CHROMOSOMAL INVERSION? (plain English)
 * --------------------------------------------------
 * Sometimes a chunk of DNA gets flipped around during cell division.
 * For example, if you had:
 *
 *   Original:  A T G [C G A T T] C C
 *                      ↑ this segment gets inverted ↑
 *   Inverted:  A T G [T T A G C] C C
 *
 * The segment "CGATT" became "TTAGC" — the same characters, backwards.
 *
 * In real biology, chromosomal inversions:
 *   - Can disrupt genes if the inversion breaks in the middle of one
 *   - Are used as markers to trace evolutionary history
 *   - Sometimes cause genetic disorders (e.g., some cases of haemophilia)
 *
 * For this project, inversion is a SEQUENCE EDITING OPERATION:
 *   "Reverse the characters between position L and R in place."
 *
 * HOW THE FINGER TREE MAKES THIS EFFICIENT
 * -----------------------------------------
 * Doing this on a plain std::string requires:
 *   std::reverse(str.begin()+L, str.begin()+R)
 * That is O(R - L) — must touch every character in the range.
 *
 * With the Finger Tree, the operation is:
 *   1. split at L                    O(log n)
 *   2. split rightPart at R-L        O(log n)
 *   3. reverse the middle segment    O(k) where k = R-L
 *      (rebuild the segment by prepending each element of the original
 *       middle onto a new tree — like reading it backwards)
 *   4. concat leftPart + reversedMiddle + restPart   O(log n)
 *
 * Steps 1, 2, 4 are O(log n).  Step 3 is still O(k) because we must
 * physically reorder k elements.  However, because the Finger Tree
 * supports fast split and concat, the STRUCTURAL work is O(log n) and
 * only the actual reversal work is O(k).  For small inversions on large
 * genomes this is significantly faster than in-place string reversal
 * which also has to shift the surrounding data.
 *
 * Advanced enhancement (not required but possible):
 *   Add a "reverse" lazy tag to Deep nodes, flipping prefix/suffix and
 *   propagating on demand.  Then inversion becomes O(log n) total
 *   including step 3 — the reversal is deferred until individual
 *   characters are actually accessed.
 */

#ifndef INVERSION_H
#define INVERSION_H

#include "Sequence.h"

class Inversion {
public:
    /*
     * invert(seq, left, right)
     * -------------------------
     * Reverses the segment seq[left .. right-1] in place.
     * The sequence outside [left, right) is left unchanged.
     *
     * Parameters:
     *   seq   — the DNA sequence (modified in place)
     *   left  — start index (0-based, inclusive)
     *   right — end index   (0-based, exclusive)
     *
     * Example:
     *   seq = "ATGCGATT"
     *   invert(seq, 3, 7)
     *   → seq becomes "ATGTTAG" ... wait, let's be precise:
     *     segment is seq[3..6] = "CGAT"  reversed = "TAGC"
     *     result: "ATG" + "TAGC" + "T" = "ATGTAGCT"
     *
     * HOW IT WORKS (step by step):
     *   1. auto [left_tree, right_tree] = split(seq.getTree(), left)
     *      left_tree  = tree for seq[0..left-1]
     *      right_tree = tree for seq[left..n-1]
     *
     *   2. auto [mid_tree, rest_tree] = split(right_tree, right - left)
     *      mid_tree   = tree for seq[left..right-1]  ← this gets reversed
     *      rest_tree  = tree for seq[right..n-1]
     *
     *   3. reversed_mid = reverseTree(mid_tree)
     *      Builds a new tree by gathering all bases from mid_tree into
     *      a string, reversing the string, then pushing each character
     *      back with pushBack().
     *
     *   4. seq.getTree() = concat(concat(left_tree, reversed_mid), rest_tree)
     *      Reassembles the full sequence.
     *
     * Time: O(log n) for splits/concats + O(k) for the reversal,
     *       where k = right - left.
     *
     * Preconditions:
     *   0 <= left <= right <= seq.length()
     *
     * Used by: genomic rearrangement simulations, testing structural
     *          integrity of the Finger Tree after complex edits.
     */
    static void invert(Sequence& seq, int left, int right);

private:
    /*
     * reverseTree(t)
     * ---------------
     * Converts a Finger Tree to a string, reverses it, and rebuilds
     * a new Finger Tree from the reversed string.
     *
     * This is an internal helper — callers should use invert() above.
     *
     * Time: O(k) where k is the number of elements in t.
     */
    static std::shared_ptr<FingerTree> reverseTree(std::shared_ptr<FingerTree> t);
};

#endif // INVERSION_H
