#ifndef INVERSION_H
#define INVERSION_H

#include "Sequence.h"

class Inversion {
public:
    // Reverse segment [l, r]
    static void invert(Sequence& seq, int left, int right);
};

#endif
