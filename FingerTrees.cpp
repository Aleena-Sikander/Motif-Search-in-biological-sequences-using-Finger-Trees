#include "FingerTree.h"
#include <stdexcept>

Base::Base(char t) : type(t) {}

Measure Base::getMeasure() const {
    return { 1 };
}

void Base::gatherBases(std::string& out) const {
    out += type;
}

Node::Node(std::vector<std::shared_ptr<Measured>> c) : children(std::move(c)) {
    m = { 0 };
    for (auto& child : children)
        m = m + child->getMeasure();
}

Measure Node::getMeasure() const {return m;}

void Node::gatherBases(std::string& out) const {
    for (auto& child : children)
        child->gatherBases(out);
}

Measure Empty::getMeasure() const { return { 0 }; }
bool Empty::isEmpty() const { return true;  }

Single::Single(std::shared_ptr<Measured> v) : val(std::move(v)) {}

Measure Single::getMeasure() const { return val->getMeasure(); }
bool Single::isEmpty() const { return false; }

Deep::Deep(std::vector<std::shared_ptr<Measured>> p, std::shared_ptr<FingerTree> mid, std::vector<std::shared_ptr<Measured>> s): prefix(std::move(p)), middle(std::move(mid)), suffix(std::move(s))
{
    m = { 0 };
    for (auto& x : prefix){ m = m + x->getMeasure()};

    m = m + middle->getMeasure();
    for (auto& x : suffix){ m = m + x->getMeasure()};
}

Measure Deep::getMeasure() const { return m;}
bool Deep::isEmpty() const { return false; }

std::shared_ptr<FingerTree> pushFront(std::shared_ptr<FingerTree> t, std::shared_ptr<Measured> v)
{
    if (t->isEmpty())
        return std::make_shared<Single>(v);

    //if tree is single then it becomes deep
    if (auto s = std::dynamic_pointer_cast<Single>(t)){
        return std::make_shared<Deep>(
            std::vector<std::shared_ptr<Measured>>{ v },
            std::make_shared<Empty>(),
            std::vector<std::shared_ptr<Measured>>{ s->val }
        );
    }

    auto d = std::dynamic_pointer_cast<Deep>(t);

    //if prefix has space then prepend v to prefix.
    if (d->prefix.size() < 4) {
        auto p = d->prefix;
        p.insert(p.begin(), v);
        return std::make_shared<Deep>(p, d->middle, d->suffix);
    }

    //rebalancing when prefix gets full, convert into node and push it into the middle
    auto overflow = std::make_shared<Node>(std::vector<std::shared_ptr<Measured>>{ d->prefix[1], d->prefix[2], d->prefix[3] });
    return std::make_shared<Deep>(
        std::vector<std::shared_ptr<Measured>>{ v, d->prefix[0] },
        pushFront(d->middle, overflow),
        d->suffix
    );
}

std::shared_ptr<FingerTree> pushBack(std::shared_ptr<FingerTree> t, std::shared_ptr<Measured> v){

    if (t->isEmpty()){
        return std::make_shared<Single>(v);}

    if (auto s = std::dynamic_pointer_cast<Single>(t)){
        return std::make_shared<Deep>(
            std::vector<std::shared_ptr<Measured>>{ s->val },
            std::make_shared<Empty>(),
            std::vector<std::shared_ptr<Measured>>{ v }
        );
    }

    auto d = std::dynamic_pointer_cast<Deep>(t);

    if (d->suffix.size() < 4) {
        auto s = d->suffix;
        s.push_back(v);
        return std::make_shared<Deep>(d->prefix, d->middle, s);
    }

    //rebalancing when sufix gets full by converting into node and pushing it in the middle
    auto overflow = std::make_shared<Node>(std::vector<std::shared_ptr<Measured>>{ d->suffix[0], d->suffix[1], d->suffix[2] });
    return std::make_shared<Deep>(
        d->prefix,
        pushBack(d->middle, overflow),
        std::vector<std::shared_ptr<Measured>>{ d->suffix[3], v }
    );
}

std::vector<std::shared_ptr<Measured>> toNodes(std::vector<std::shared_ptr<Measured>> b){.
    //grouping into 3, 2 only when 2/4 are left to avoid 1-node
    std::vector<std::shared_ptr<Measured>> result;

    for (size_t i = 0; i < b.size(); ) {
        size_t rem = b.size() - i;

        if (rem == 2 || rem == 3) {
            // Take all remaining elements as one node (2 or 3 children).
            result.push_back(std::make_shared<Node>(
                std::vector<std::shared_ptr<Measured>>(b.begin() + i, b.end())
            ));
            break;
        }

        if (rem == 4) {
            //spliting 4 into two 2-nodes 
            result.push_back(std::make_shared<Node>(std::vector<std::shared_ptr<Measured>>{ b[i], b[i + 1] }));
            result.push_back(std::make_shared<Node>(std::vector<std::shared_ptr<Measured>>{ b[i + 2], b[i + 3] }));
            break;
        }

        //takeing 3 elements as one node
        result.push_back(std::make_shared<Node>(std::vector<std::shared_ptr<Measured>>{ b[i], b[i + 1], b[i + 2] }));
        i += 3;
    }
    return result;
}

std::shared_ptr<FingerTree> app3(std::shared_ptr<FingerTree> t1,std::vector<std::shared_ptr<Measured>> ts, std::shared_ptr<FingerTree>t2)
{
    if (t1->isEmpty()) {
        // Push all ts elements and then t2's content onto an empty left.
        std::shared_ptr<FingerTree> r = t2;
        for (auto it = ts.rbegin(); it != ts.rend(); ++it){
            r = pushFront(r, *it);
        return r;}
    }

    if (t2->isEmpty()) {
        std::shared_ptr<FingerTree> r = t1;
        for (auto& x : ts){
            r = pushBack(r, x);
        return r;}
    }

    //pushinf single and recursing with an empty placeholder.
    if (auto s1 = std::dynamic_pointer_cast<Single>(t1))
        return pushFront(app3(std::make_shared<Empty>(), ts, t2), s1->val);

    if (auto s2 = std::dynamic_pointer_cast<Single>(t2))
        return pushBack(app3(t1, ts, std::make_shared<Empty>()), s2->val);

    //deep+deep: combining suffix,middle and prefix then converting into nodes
    auto d1 = std::dynamic_pointer_cast<Deep>(t1);
    auto d2 = std::dynamic_pointer_cast<Deep>(t2);

    // Build the combined middle buffer:  d1.suffix + ts + d2.prefix
    std::vector<std::shared_ptr<Measured>> b = d1->suffix;
    b.insert(b.end(), ts.begin(), ts.end());
    b.insert(b.end(), d2->prefix.begin(), d2->prefix.end());

    return std::make_shared<Deep>(
        d1->prefix,
        app3(d1->middle, toNodes(b), d2->middle),
        d2->suffix
    );
}

/* =========================================================
 *  concat — public O(log n) join of two trees
 * ========================================================= */

std::shared_ptr<FingerTree> concat(std::shared_ptr<FingerTree> t1,
                                   std::shared_ptr<FingerTree> t2)
{
    // app3 with an empty middle buffer is exactly concat.
    return app3(t1, {}, t2);
}


static std::pair<std::shared_ptr<FingerTree>, std::shared_ptr<FingerTree>>
splitMeasured(std::shared_ptr<Measured> elem, size_t i);


static std::pair<std::shared_ptr<FingerTree>, std::shared_ptr<FingerTree>>
splitDigit(const std::vector<std::shared_ptr<Measured>>& digit, size_t i);


static std::pair<std::shared_ptr<FingerTree>, std::shared_ptr<FingerTree>>
splitMeasured(std::shared_ptr<Measured> elem, size_t i)
{
    //If Base the split will happen at the leaf level
    if (auto b = std::dynamic_pointer_cast<Base>(elem)) {
        if (i == 0) //base going to right tree
            return { std::make_shared<Empty>(), std::make_shared<Single>(b) };
        else
            return { std::make_shared<Single>(b), std::make_shared<Empty>() };
    }

    //recursing into its children if it's a node
    auto nd = std::dynamic_pointer_cast<Node>(elem);

    size_t offset = 0;
    for (size_t ci = 0; ci < nd->children.size(); ++ci) {
        size_t childSize = nd->children[ci]->getMeasure().size;

        if (i < offset + childSize) {
            auto [childLeft, childRight] = splitMeasured(nd->children[ci], i - offset);

            std::shared_ptr<FingerTree> left = childLeft;
            for (size_t j = 0; j < ci; ++j)
                left = pushBack(left, nd->children[j]);

            std::shared_ptr<FingerTree> right = childRight;
            for (size_t j = ci + 1; j < nd->children.size(); ++j)
                right = pushBack(right, nd->children[j]);

            return { left, right };
        }
        offset += childSize;
    }

    //if i >= total size of this node then everything goes left
    std::shared_ptr<FingerTree> left = std::make_shared<Empty>();
    for (auto& c : nd->children) left = pushBack(left, c);
    return { left, std::make_shared<Empty>() };
}

static std::pair<std::shared_ptr<FingerTree>, std::shared_ptr<FingerTree>>
splitDigit(const std::vector<std::shared_ptr<Measured>>& digit, size_t i)
{
    size_t offset = 0;
    for (size_t di = 0; di < digit.size(); ++di) {
        size_t elemSize = digit[di]->getMeasure().size;

        if (i < offset + elemSize) {
            // Split falls inside digit[di].
            auto [elemLeft, elemRight] = splitMeasured(digit[di], i - offset);

            // All digits before di go into left; all after di go into right.
            std::shared_ptr<FingerTree> left = elemLeft;
            for (size_t j = 0; j < di; ++j)
                left = pushBack(left, digit[j]);

            std::shared_ptr<FingerTree> right = elemRight;
            for (size_t j = di + 1; j < digit.size(); ++j)
                right = pushBack(right, digit[j]);

            return { left, right };
        }
        offset += elemSize;
    }

    // i >= total size of digi so everything goes left.
    std::shared_ptr<FingerTree> left = std::make_shared<Empty>();
    for (auto& d : digit) left = pushBack(left, d);
    return { left, std::make_shared<Empty>() };
}

std::pair<std::shared_ptr<FingerTree>, std::shared_ptr<FingerTree>>
split(std::shared_ptr<FingerTree> t, size_t i)
{
    if (t->isEmpty())
        return { std::make_shared<Empty>(), std::make_shared<Empty>() };

    if (i >= t->getMeasure().size)
        return { t, std::make_shared<Empty>() };

    if (i == 0)
        return { std::make_shared<Empty>(), t };

    if (auto s = std::dynamic_pointer_cast<Single>(t)) {
        if (i >= s->val->getMeasure().size)
            return { t, std::make_shared<Empty>() };
        return { std::make_shared<Empty>(), t };
    }

    auto d = std::dynamic_pointer_cast<Deep>(t);

    size_t prefixSize = 0;
    for (auto& x : d->prefix) prefixSize += x->getMeasure().size;

    if (i < prefixSize) {
        // Split is entirely within the prefix.
        auto [prefLeft, prefRight] = splitDigit(d->prefix, i);
        // Reassemble: right part = prefRight concatenated with middle + suffix.
        std::shared_ptr<FingerTree> rightTree = prefRight;
        // Push middle elements onto right
        std::string midStr;
        getSequence(d->middle, midStr); // flatten middle
        // Actually we need to concat properly:
        // right = concat(prefRight, concat(middleAsTree, suffixTree))
        std::shared_ptr<FingerTree> suffTree = std::make_shared<Empty>();
        for (auto& x : d->suffix) suffTree = pushBack(suffTree, x);
        rightTree = concat(prefRight, concat(d->middle, suffTree));
        return { prefLeft, rightTree };
    }

    size_t afterPrefix = prefixSize;

    // --- Check if the split falls in the MIDDLE ---
    size_t middleSize = d->middle->getMeasure().size;

    if (i < afterPrefix + middleSize) {
        // Split is inside the middle spine.
        size_t middleOffset = i - afterPrefix;
        auto [midLeft, midRight] = split(d->middle, middleOffset);

        // Left part: prefix + midLeft
        std::shared_ptr<FingerTree> leftTree = std::make_shared<Empty>();
        for (auto& x : d->prefix) leftTree = pushBack(leftTree, x);
        leftTree = concat(leftTree, midLeft);

        // Right part: midRight + suffix
        std::shared_ptr<FingerTree> suffTree = std::make_shared<Empty>();
        for (auto& x : d->suffix) suffTree = pushBack(suffTree, x);
        std::shared_ptr<FingerTree> rightTree = concat(midRight, suffTree);

        return { leftTree, rightTree };
    }

    // --- The split falls in the SUFFIX ---
    size_t suffixOffset = i - afterPrefix - middleSize;
    auto [sufLeft, sufRight] = splitDigit(d->suffix, suffixOffset);

    // Left part: prefix + middle + sufLeft
    std::shared_ptr<FingerTree> leftTree = std::make_shared<Empty>();
    for (auto& x : d->prefix) leftTree = pushBack(leftTree, x);
    leftTree = concat(leftTree, concat(d->middle, sufLeft));

    return { leftTree, sufRight };
}

void getSequence(std::shared_ptr<FingerTree> t, std::string& out)
{
    if (t->isEmpty()) return;

    if (auto s = std::dynamic_pointer_cast<Single>(t)) {
        // Ask the single element (Base or Node) to write itself into out.
        s->val->gatherBases(out);
        return;
    }
    auto d = std::dynamic_pointer_cast<Deep>(t);

    //traversing
    for (auto& x : d->prefix) x->gatherBases(out);
    getSequence(d->middle, out);
    for (auto& x : d->suffix) x->gatherBases(out);
}