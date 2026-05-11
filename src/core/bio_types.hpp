#pragma once
#include <cstddef>

// Stores two running totals for any subtree: how many chars, how many are G/C
struct Measure {
    size_t length   = 0;
    size_t gc_count = 0;

    // combining two measures just adds the numbers
    Measure operator+(const Measure& other) const {
        return { length + other.length, gc_count + other.gc_count };
    }

    // identity value — used as a starting point before any elements are added
    static Measure zero() { return {0, 0}; }
};

// tells the tree how to measure one element of type T
template <typename T, typename M>
struct Measured;

// specialisation for a single DNA character
template <>
struct Measured<char, Measure> {
    static Measure measure(char c) {
        size_t gc = (c == 'G' || c == 'C' || c == 'g' || c == 'c') ? 1 : 0;
        return {1, gc};  // every char adds 1 to length; G/C also adds 1 to gc_count
    }
};
