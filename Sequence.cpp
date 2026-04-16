#include "Sequence.h"

Sequence::Sequence() : tree() {}

Sequence::Sequence(const std::string& str) : tree() {
    for (char c : str) {
        tree.insert(tree.length(), c);
    }
}

void Sequence::insert(int index, char value) {
    tree.insert(index, value);
}

void Sequence::remove(int index) {
    tree.remove(index);
}

void Sequence::update(int index, char value) {
    tree.update(index, value);
}

int Sequence::length() const {
    return tree.length();
}

std::string Sequence::toString() const {
    std::string result;
    for (int i = 0; i < tree.length(); ++i) {
        result += tree.get(i);
    }
    return result;
}

FingerTree<char>& Sequence::getTree() {
    return tree;
}
