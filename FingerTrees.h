//generic fingertree stucture
#ifndef FINGERTREE_H
#define FINGERTREE_H

#include <memory>
#include <vector>
#include <string>

template <typename T>
class FingerTree {
private:
    struct Node {
        T value; //DNA character
        int size; //subtree size

        std::shared_ptr<Node> left;
        std::shared_ptr<Node> right;

        Node(T val);
    };

    std::shared_ptr<Node> root;

public:
    FingerTree();

    void append(const T& value);
    void prepend(const T& value);

    FingerTree<T> split(int index);
    void concat(FingerTree<T>& other);

    T lookup(int index) const;
    int getSize() const;

    //converting to list
    std::vector<T> toVector() const;
};

#endif
