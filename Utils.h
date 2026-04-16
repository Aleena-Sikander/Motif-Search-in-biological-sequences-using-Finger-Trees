/*
 * Utils.h
 * -------
 * Shared utility functions used across the entire project.
 *
 * WHY HAVE A SEPARATE UTILS FILE?
 * --------------------------------
 * Several modules need the same small helper operations:
 *   - "Is this character a valid DNA base?"
 *   - "Convert this string to uppercase"
 *
 * Rather than copy-pasting these into MotifSearch.cpp, GCContent.cpp,
 * Sequence.cpp etc., we define them once here.  This avoids bugs where
 * one copy gets updated and another doesn't.
 *
 * It also keeps the other headers clean — they describe WHAT they do
 * (search for motifs, calculate GC content) without being cluttered by
 * low-level character validation details.
 */

#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <cctype>

class Utils {
public:
    /*
     * isValidDNA(c)
     * --------------
     * Returns true if `c` is a legal DNA base character: A, T, G, or C
     * (case-insensitive).  Returns false for anything else (digits,
     * spaces, 'N', 'U', etc.).
     *
     * WHY WE NEED THIS:
     * When a user types a sequence at the command line, they might
     * accidentally include spaces, newlines, or invalid letters.
     * Calling isValidDNA() on each character during Sequence construction
     * catches these errors early, before they corrupt the tree or produce
     * nonsensical motif search results.
     *
     * It also prevents the Finger Tree from storing junk data — every
     * Base leaf in the tree should be a real DNA nucleotide.
     *
     * Example:
     *   isValidDNA('A')  → true
     *   isValidDNA('g')  → true   (lowercase accepted)
     *   isValidDNA('N')  → false  (N means "unknown base" in FASTA format;
     *                              not supported in this project)
     *   isValidDNA('5')  → false
     *   isValidDNA(' ')  → false
     */
    static bool isValidDNA(char c);

    /*
     * toUpper(str)
     * -------------
     * Returns a new string with every character converted to uppercase.
     *
     * WHY WE NEED THIS:
     * Motif search should be case-insensitive: searching for "atg" and
     * "ATG" should find the same matches.  Rather than handling both
     * cases in the matching logic, we normalise everything to uppercase
     * once at the entry point.
     *
     * Both the input sequence and the query motif are uppercased before
     * comparison, so the rest of the code only ever deals with uppercase.
     *
     * Example:
     *   toUpper("atgc") → "ATGC"
     *   toUpper("AtGc") → "ATGC"
     *   toUpper("ATGC") → "ATGC"   (no change, already uppercase)
     */
    static std::string toUpper(const std::string& str);

    /*
     * isValidMotif(motif)
     * --------------------
     * Returns true if every character in `motif` is a valid DNA base.
     * Validates the entire query string before starting the search.
     *
     * WHY WE NEED THIS:
     * A motif with an invalid character (e.g., "AT5C") can never match
     * a real DNA sequence.  Detecting this upfront avoids wasting time
     * running the search and confuses the user less ("invalid motif"
     * is a clearer error than "0 matches found").
     *
     * Example:
     *   isValidMotif("GAATTC") → true
     *   isValidMotif("GN4TC")  → false
     *   isValidMotif("")       → false  (empty motif is meaningless)
     */
    static bool isValidMotif(const std::string& motif);
};

#endif // UTILS_H
