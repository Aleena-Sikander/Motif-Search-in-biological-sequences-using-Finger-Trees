#ifndef FINGERTREE_H
#define FINGERTREE_H

#include <memory>
#include <vector>
#include <string>
#include <functional>

//number of DNA bases in a subtree; used for O(log n) navigation and splits
struct Measure {
    size_t size; // total number of DNA characters
    Measure operator+(const Measure& other) const { return {size + other.size}; }
};

struct Measured {
    virtual ~Measured() = default;

    //size of the element
    virtual Measure getMeasure() const = 0;

    //gathering DNA bases into a singgle string
    virtual void gatherBases(std::string& out) const = 0;
};

struct Base : public Measured {
    char type; // one of: A, T, G, C 

    Base(char t);
    Measure getMeasure() const override; //always one bc base
    void gatherBases(std::string& out) const override; //gathering single  char to the output str
};

struct Node : public Measured {
    std::vector<std::shared_ptr<Measured>> children;
    Measure m; //grouping into one node to keep tree balanced

    Node(std::vector<std::shared_ptr<Measured>> c);
    Measure getMeasure() const override;
    void gatherBases(std::string& out) const override;
};

class FingerTree {
public:
    virtual ~FingerTree() = default;
    virtual Measure getMeasure() const = 0; //total number of DNA characters stored in this tree
    virtual bool isEmpty() const = 0;
};

class Empty : public FingerTree {
public:
    Measure getMeasure() const override; //always 0 bc empty
    bool isEmpty() const override; //always true bc empty
};

class Single : public FingerTree {
public:
    std::shared_ptr<Measured> val; 
    Single(std::shared_ptr<Measured> v);
    Measure getMeasure() const override;
    bool isEmpty() const override; // always false
};

class Deep : public FingerTree {
public:
    std::vector<std::shared_ptr<Measured>> prefix;
    std::shared_ptr<FingerTree>middle;
    std::vector<std::shared_ptr<Measured>> suffix;
    Measure m; //total measure

    Deep(std::vector<std::shared_ptr<Measured>> p,std::shared_ptr<FingerTree> mid,std::vector<std::shared_ptr<Measured>> s);
    Measure getMeasure() const override;
    bool isEmpty() const override;
};

//core tree ops: return a new tree instead of changing the initial one

std::shared_ptr<FingerTree> pushFront(std::shared_ptr<FingerTree> t, std::shared_ptr<Measured> v);

std::shared_ptr<FingerTree> pushBack(std::shared_ptr<FingerTree> t, std::shared_ptr<Measured> v);

std::vector<std::shared_ptr<Measured>> toNodes(std::vector<std::shared_ptr<Measured>> b); //joining into nodes for the spine

std::shared_ptr<FingerTree> app3(std::shared_ptr<FingerTree> t1, std::vector<std::shared_ptr<Measured>> ts, std::shared_ptr<FingerTree> t2); //joining trees

std::shared_ptr<FingerTree> concat(std::shared_ptr<FingerTree> t1, std::shared_ptr<FingerTree> t2);//joining sequences

std::pair<std::shared_ptr<FingerTree>, std::shared_ptr<FingerTree>>split(std::shared_ptr<FingerTree> t, size_t i);

void getSequence(std::shared_ptr<FingerTree> t, std::string& out);

#endif 