//since GC bonds are stronger than AT bonds so GC content is important for stability of DNA and for pattern recognition
#ifndef GCCONTENT_H
#define GCCONTENT_H

#include "Sequence.h"

class GCContent {
public:
    //calculating GC % 
    static double calculate(Sequence& seq, int left, int right);
};

#endif
