#ifndef INVERSION_H
#define INVERSION_H

#include "Sequence.h"

class Inversion {
public:
    //reversing a part of DNA
    static void invert(Sequence& seq, int left, int right);
};

#endif
