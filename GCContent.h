//since GC bonds are stronger than AT bonds so GC content is important for stability of DNA and for pattern recognition
/*
 * GCContent.h
 * ------------
 * Calculates the GC content (%) of any range within a DNA sequence.
 *
 * WHAT IS GC CONTENT? (plain English)
 * -------------------------------------
 * DNA is made of four bases: Adenine (A), Thymine (T), Guanine (G),
 * and Cytosine (C).  A always pairs with T; G always pairs with C.
 *
 * G-C pairs are held together by THREE hydrogen bonds.
 * A-T pairs are held together by only TWO hydrogen bonds.
 *
 * So a region of DNA with more G/C bases is physically STRONGER —
 * it needs more heat to "melt" the two strands apart.
 *
 * GC content = (number of G + C bases) / (total bases) × 100%
 *
 * WHY IS THIS USEFUL?
 * --------------------
 * In bioinformatics, GC content of a region tells you:
 *   - How thermally stable that region is (high GC = harder to open)
 *   - Whether a region might be a CpG island (gene promoter indicator)
 *   - Whether a motif you found is likely to be functionally real
 *     (some organisms have very high or low genome-wide GC content)
 *
 * For this project, GC content is a RANGE QUERY:
 *   "What percentage of bases between position L and R are G or C?"
 *
 * HOW THE FINGER TREE HELPS
 * --------------------------
 * A naive approach converts the whole sequence to a string and counts
 * in the range — O(n) every time.
 *
 * With the Finger Tree we use split/concat:
 *   1. split at L              → leftPart + rightPart
 *   2. split rightPart at R-L  → middle (the range) + restPart
 *   3. traverse middle to count G/C
 *   4. concat leftPart + middle + restPart to restore the tree
 *
 * Steps 1, 2, 4 are each O(log n).  Step 3 is O(k) where k = R - L.
 * For small ranges this is much faster than always processing all n bases.
 *
 * You can also extend the Measure struct to CACHE the GC count at each
 * node — then step 3 becomes O(log n) too, making the entire query O(log n).
 * (This is an advanced enhancement noted in the project proposal under
 *  "metadata in nodes".)
 */

#ifndef GCCONTENT_H
#define GCCONTENT_H

#include "Sequence.h"

class GCContent {
public:
    /*
     * calculate(seq, left, right)
     * ----------------------------
     * Returns the GC content percentage for the subsequence
     * seq[left .. right-1] (right is exclusive, like std::string::substr).
     *
     * Parameters:
     *   seq   — the full DNA sequence
     *   left  — start index (0-based, inclusive)
     *   right — end index   (0-based, exclusive)
     *
     * Returns:
     *   A double in [0.0, 100.0] representing the GC percentage.
     *   Returns 0.0 if the range is empty or invalid.
     *
     * Example:
     *   seq = "ATGCGCTA"
     *   calculate(seq, 2, 6) operates on "GCGC" → 4/4 × 100 = 100.0%
     *   calculate(seq, 0, 4) operates on "ATGC" → 2/4 × 100 = 50.0%
     *
     * HOW IT USES THE FINGER TREE:
     *   Calls seq.substring(left, right) to extract just the target
     *   range using split+concat (O(log n)), then iterates over the
     *   extracted string to count G and C (O(k) where k = right - left).
     *
     * Preconditions:
     *   0 <= left <= right <= seq.length()
     */
    static double calculate(const Sequence& seq, int left, int right);

    /*
     * countGC(str)
     * -------------
     * Helper: counts the number of G and C characters in a plain string.
     * Called internally by calculate() after extracting the substring.
     *
     * Case-insensitive — 'g' and 'c' count the same as 'G' and 'C'.
     *
     * Example: countGC("ATGCGC") → 4
     */
    static int countGC(const std::string& str);
};

#endif // GCCONTENT_H
