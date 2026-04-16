#include "FingerTrees.h"

template <typename T>
FingerTree<T>::Node::Node(T val) : value(val), size(1), left(nullptr), right(nullptr) {};

template <typename T>
FingerTree<T>::FingerTree() : root(nullptr) {}

template <typename T>
void FingerTree<T>::append(const T& value) {
    auto newNode = std::make_shared<Node>(value);
    if (!root) {
        root = newNode;
    } else {
        auto current = root;
        while (current->right) {
            current = current->right;
        }
        current->right = newNode;
        current->size += 1;
    }
}

template <typename T>
void FingerTree<T>::prepend(const T& value) {
    auto newNode = std::make_shared<Node>(value);
    if (!root) {
        root = newNode;
    } else {
        newNode->right = root;
        root->size += 1;
        root = newNode;
    }
}

template <typename T>
FingerTree<T> FingerTree<T>::split(int index) {
    // Implement split logic here
    // This is a placeholder implementation
    return FingerTree<T>();
}

template <typename T>
void FingerTree<T>::concat(FingerTree<T>& other) {
    // Implement concat logic here
    // This is a placeholder implementation
}

template <typename T>
T FingerTree<T>::lookup(int index) const {
    // Implement lookup logic here
    // This is a placeholder implementation
    return T();
}

template <typename T>
int FingerTree<T>::getSize() const {
    if(root) {
        return root->size;
    } else {
        return 0;
    }
}   

template <typename T>
std::vector<T> FingerTree<T>::toVector() const {
    std::vector<T> result;
    auto current = root;
    while (current) {
        result.push_back(current->value);
        current = current->right;
    }
    return result;
}