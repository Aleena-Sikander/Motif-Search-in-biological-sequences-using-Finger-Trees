#ifndef GCCONTENT_H
#define GCCONTENT_H

#include "Sequence.h"

class GCContent {
public:
    // Calculate GC % in range [l, r]
    static double calculate(Sequence& seq, int left, int right);
};

#endif
