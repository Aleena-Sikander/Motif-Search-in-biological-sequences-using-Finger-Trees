#ifndef SEQUENCE_H
#define SEQUENCE_H

#include "FingerTree.h"
#include <string>

class Sequence {
private:
    FingerTree<char> tree;

public:
    Sequence();
    Sequence(const std::string& str);

    // Gene editing
    void insert(int index, char value);
    void remove(int index);
    void update(int index, char value);

    // Utilities
    int length() const;
    std::string toString() const;

    FingerTree<char>& getTree();
};

#endif
