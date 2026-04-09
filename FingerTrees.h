#ifndef FINGERTREE_H
#define FINGERTREE_H

#include <memory>
#include <vector>
#include <string>

template <typename T>
class FingerTree {
private:
    struct Node {
        T value;
        int size;

        std::shared_ptr<Node> left;
        std::shared_ptr<Node> right;

        Node(T val);
    };

    std::shared_ptr<Node> root;

public:
    FingerTree();

    // Basic operations
    void append(const T& value);
    void prepend(const T& value);

    // Split & concatenate
    FingerTree<T> split(int index);
    void concatenate(FingerTree<T>& other);

    // Access
    T lookup(int index) const;
    int getSize() const;

    // Convert to list/string
    std::vector<T> toVector() const;
};

#endif
