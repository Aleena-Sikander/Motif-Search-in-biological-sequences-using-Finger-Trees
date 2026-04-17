#include <iostream>
#include "Sequence.h"

using namespace std;

int main() {
    std::shared_ptr<FingerTree> s1 = std::make_shared<Empty>(), s2 = std::make_shared<Empty>();
    for (char c : std::string("ATGC")) s1 = pushBack(s1, std::make_shared<Base>(c));
    for (char c : std::string("GCTA")) s2 = pushBack(s2, std::make_shared<Base>(c));

    std::string out1, out2, outComb;
    getSequence(s1, out1); getSequence(s2, out2);
    auto combined = concat(s1, s2);
    getSequence(combined, outComb);

    std::cout << "Seq 1: " << out1 << " (Size: " << s1->getMeasure().size << ")\n";
    std::cout << "Seq 2: " << out2 << " (Size: " << s2->getMeasure().size << ")\n";
    std::cout << "Combined: " << outComb << " (Size: " << combined->getMeasure().size << ")\n";
    return 0;
}