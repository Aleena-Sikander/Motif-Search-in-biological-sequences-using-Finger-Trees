#ifndef MOTIFSEARCH_H
#define MOTIFSEARCH_H

#include "Sequence.h"
#include <vector>
#include <string>

class MotifSearch {
public:
    // Returns all starting indices of motif
    static std::vector<int> findMotif(Sequence& seq, const std::string& motif);
};

#endif
