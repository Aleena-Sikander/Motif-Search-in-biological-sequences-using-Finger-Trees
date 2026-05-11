#pragma once

// 2-3 Finger Tree  —  Hinze & Paterson (2006)
//
// A finger tree stores a sequence with:
//   push_front / push_back  in amortised O(1)
//   split / concat          in O(log n)
//   measure()               in O(1)  — cached at every node
//
// Shape of a Deep node:
//   [ prefix (1-4 elems) | spine | suffix (1-4 elems) ]
//
// The spine is another finger tree that stores Nodes (groups of
// 2 or 3 elements) instead of raw elements. That nesting is what
// makes it a "2-3 finger tree" and gives the O(log n) depth.
//
// --- Why the depth number D? ---
// Writing  FingerTree<Node<T>>  inside  FingerTree<T>  would make
// the compiler instantiate an infinite chain of types.  We break
// that by adding D: FT<T,M,D> has a spine of type FT<Node<T,M>,M,D-1>.
// At D=0 there is no spine. Depth 12 supports ~89 million elements.

#include "bio_types.hpp"
#include "node.hpp"
#include <memory>
#include <vector>
#include <cassert>
#include <stdexcept>
#include <functional>

template <typename T, typename M, int D> class FT;

// ---------------------------------------------------------------
// Base case (D = 0): holds at most one element, no spine.
// ---------------------------------------------------------------
template <typename T, typename M>
class FT<T, M, 0> {
    bool has = false;
    M    anno = M::zero();
    T    val{};
public:
    FT() = default;
    explicit FT(const T& x) : has(true), val(x) { anno = Measured<T,M>::measure(x); }

    bool     empty()   const { return !has; }
    const M& measure() const { return anno; }

    template <typename F> void for_each(F&& f) const { if (has) f(val); }

    FT push_back (const T& x) const { if (!has) return FT(x); throw std::overflow_error("tree depth exhausted"); }
    FT push_front(const T& x) const { if (!has) return FT(x); throw std::overflow_error("tree depth exhausted"); }

    struct ViewLeft  { T head; FT tail; };
    struct ViewRight { FT init; T last; };
    ViewLeft  view_left()  const { return {val, FT{}}; }
    ViewRight view_right() const { return {FT{}, val}; }

    struct Split { FT left; T pivot; FT right; };
    // only one element — it always becomes the pivot
    Split split(const std::function<bool(const M&)>&, M = M::zero()) const {
        return {FT{}, val, FT{}};
    }

    static FT concat(const FT& a, const FT& b) {
        FT r = a;
        b.for_each([&](const T& x){ r = r.push_back(x); });
        return r;
    }
    bool check_invariants() const { return true; }
};

// ---------------------------------------------------------------
// Recursive case (D > 0)
// ---------------------------------------------------------------
template <typename T, typename M, int D>
class FT {
    using NodeT  = Node<T, M>;        // a group of 2-3 elements
    using SpineT = FT<NodeT, M, D-1>; // the spine stores nodes, not raw elements

    enum class Shape { Empty, Single, Deep };

    Shape  shape  = Shape::Empty;
    M      anno   = M::zero();   // cached total measure of this whole tree
    T      solo{};               // the element when shape == Single

    // shape == Deep:
    std::vector<T>             pre;    // prefix:  1-4 elements
    std::shared_ptr<SpineT>    mid;    // spine:   a deeper finger tree of nodes
    std::vector<T>             suf;    // suffix:  1-4 elements

    // recompute the total measure from parts
    void recalc() {
        anno = M::zero();
        if (shape == Shape::Single) { anno = Measured<T,M>::measure(solo); return; }
        if (shape == Shape::Deep) {
            for (const T& x : pre) anno = anno + Measured<T,M>::measure(x);
            if (mid)               anno = anno + mid->measure();
            for (const T& x : suf) anno = anno + Measured<T,M>::measure(x);
        }
    }

    static FT deep(std::vector<T> p, std::shared_ptr<SpineT> m, std::vector<T> s) {
        assert(p.size() >= 1 && p.size() <= 4);
        assert(s.size() >= 1 && s.size() <= 4);
        FT t; t.shape = Shape::Deep; t.pre = p; t.mid = m; t.suf = s; t.recalc();
        return t;
    }

    // group a flat list of elements into nodes of size 2 or 3
    static std::vector<NodeT> to_nodes(const std::vector<T>& xs) {
        std::vector<NodeT> out;
        size_t i = 0, n = xs.size();
        while (i < n) {
            size_t left = n - i;
            if (left == 2 || left == 4) { out.emplace_back(std::vector<T>{xs[i], xs[i+1]}); i += 2; }
            else                        { out.emplace_back(std::vector<T>{xs[i], xs[i+1], xs[i+2]}); i += 3; }
        }
        return out;
    }

    // build a tree by pushing a small list of elements one by one
    static FT from_list(const std::vector<T>& xs) {
        FT t;
        for (const T& x : xs) t = t.push_back(x);
        return t;
    }

    // reassemble a tree from prefix + spine + suffix (any part may be empty)
    static FT reassemble(const std::vector<T>& p,
                         const std::shared_ptr<SpineT>& m,
                         const std::vector<T>& s)
    {
        FT t;
        for (const T& x : p) t = t.push_back(x);
        if (m) m->for_each([&](const NodeT& nd){ for (const T& e : nd.elems) t = t.push_back(e); });
        for (const T& x : s) t = t.push_back(x);
        return t;
    }

    // concat helper: joins a, some bridging elements, and b
    // (this is "app3" from the Hinze & Paterson paper)
    static FT join3(const FT& a, std::vector<T> bridge, const FT& b) {
        if (a.empty()) {
            FT r = b;
            for (auto it = bridge.rbegin(); it != bridge.rend(); ++it) r = r.push_front(*it);
            return r;
        }
        if (b.empty()) { FT r = a; for (const T& x : bridge) r = r.push_back(x); return r; }
        if (a.shape == Shape::Single) {
            FT r = b;
            for (auto it = bridge.rbegin(); it != bridge.rend(); ++it) r = r.push_front(*it);
            return r.push_front(a.solo);
        }
        if (b.shape == Shape::Single) {
            FT r = a; for (const T& x : bridge) r = r.push_back(x); return r.push_back(b.solo);
        }
        // both Deep: pack a.suffix + bridge + b.prefix into nodes,
        // then recursively join the two spines with those nodes in between
        std::vector<T> merged;
        for (const T& x : a.suf)    merged.push_back(x);
        for (const T& x : bridge)   merged.push_back(x);
        for (const T& x : b.pre)    merged.push_back(x);

        auto nodes = to_nodes(merged);
        SpineT new_spine = a.mid ? *a.mid : SpineT{};
        for (const NodeT& nd : nodes) new_spine = new_spine.push_back(nd);
        if (b.mid && !b.mid->empty()) {
            std::vector<NodeT> bn;
            b.mid->for_each([&](const NodeT& n){ bn.push_back(n); });
            for (const NodeT& nd : bn) new_spine = new_spine.push_back(nd);
        }
        return deep(a.pre, std::make_shared<SpineT>(new_spine), b.suf);
    }

public:
    FT() = default;
    explicit FT(const T& x) : shape(Shape::Single), solo(x) { anno = Measured<T,M>::measure(x); }

    bool     empty()   const { return shape == Shape::Empty; }
    const M& measure() const { return anno; }

    // ---- push_back: add to the right end   O(1) amortised ----
    // Normal case: suffix has room, just append.
    // Overflow:    suffix is full (4 elements), pack the first 3 into a
    //              Node and push that Node into the spine.
    FT push_back(const T& x) const {
        if (shape == Shape::Empty)  return FT(x);
        if (shape == Shape::Single) return deep({solo}, std::make_shared<SpineT>(), {x});
        if (suf.size() < 4) { auto s = suf; s.push_back(x); return deep(pre, mid, s); }

        // suffix full — overflow 3 elements into the spine as a Node
        NodeT overflow({suf[0], suf[1], suf[2]});
        auto new_mid = mid ? std::make_shared<SpineT>(mid->push_back(overflow))
                           : std::make_shared<SpineT>(SpineT().push_back(overflow));
        return deep(pre, new_mid, {suf[3], x});
    }

    // ---- push_front: add to the left end   O(1) amortised ----
    FT push_front(const T& x) const {
        if (shape == Shape::Empty)  return FT(x);
        if (shape == Shape::Single) return deep({x}, std::make_shared<SpineT>(), {solo});
        if (pre.size() < 4) {
            std::vector<T> p = {x};
            p.insert(p.end(), pre.begin(), pre.end());
            return deep(p, mid, suf);
        }
        NodeT overflow({pre[0], pre[1], pre[2]});
        auto new_mid = mid ? std::make_shared<SpineT>(mid->push_front(overflow))
                           : std::make_shared<SpineT>(SpineT().push_front(overflow));
        return deep({x, pre[3]}, new_mid, suf);
    }

    // ---- view_left / view_right ----
    struct ViewLeft  { T head; FT tail; };
    struct ViewRight { FT init; T last; };

    ViewLeft view_left() const {
        if (shape == Shape::Empty)  throw std::runtime_error("view_left: empty tree");
        if (shape == Shape::Single) return {solo, FT{}};
        T h = pre.front();
        if (pre.size() > 1) { return {h, deep({pre.begin()+1, pre.end()}, mid, suf)}; }
        if (mid && !mid->empty()) {
            auto vl = mid->view_left();  // borrow a node from the spine
            return {h, deep(vl.head.elems, std::make_shared<SpineT>(vl.tail), suf)};
        }
        if (suf.size() == 1) return {h, FT(suf[0])};
        return {h, deep({suf.front()}, std::make_shared<SpineT>(), {suf.begin()+1, suf.end()})};
    }

    ViewRight view_right() const {
        if (shape == Shape::Empty)  throw std::runtime_error("view_right: empty tree");
        if (shape == Shape::Single) return {FT{}, solo};
        T l = suf.back();
        if (suf.size() > 1) { return {deep(pre, mid, {suf.begin(), suf.end()-1}), l}; }
        if (mid && !mid->empty()) {
            auto vr = mid->view_right();
            return {deep(pre, std::make_shared<SpineT>(vr.init), vr.last.elems), l};
        }
        if (pre.size() == 1) return {FT(pre[0]), l};
        return {deep({pre.begin(), pre.end()-1}, std::make_shared<SpineT>(), {pre.back()}), l};
    }

    // ---- for_each: visit every element left to right   O(n) ----
    template <typename F>
    void for_each(F&& f) const {
        if (shape == Shape::Empty)  return;
        if (shape == Shape::Single) { f(solo); return; }
        for (const T& x : pre) f(x);
        if (mid) mid->for_each([&](const NodeT& nd){ for (const T& e : nd.elems) f(e); });
        for (const T& x : suf) f(x);
    }

    // ---- split: find the position where pred first becomes true  O(log n) ----
    //
    // pred is a function on the accumulated measure.
    // Returns { left, pivot, right } where:
    //   - left  contains everything before the split point
    //   - pivot is the element exactly at the split point
    //   - right contains everything after
    //
    // This works by checking prefix, then recursing into the spine
    // (which is itself a finger tree), then checking suffix.
    struct Split { FT left; T pivot; FT right; };

    Split split(const std::function<bool(const M&)>& pred, M acc = M::zero()) const {
        if (shape == Shape::Single) return {FT{}, solo, FT{}};
        if (shape == Shape::Empty)  throw std::runtime_error("split: empty tree");

        // check prefix left to right
        {
            M run = acc;
            for (size_t i = 0; i < pre.size(); i++) {
                run = run + Measured<T,M>::measure(pre[i]);
                if (pred(run)) {
                    // split point is pre[i]
                    std::vector<T> left_part(pre.begin(), pre.begin() + i);
                    std::vector<T> right_part(pre.begin() + i + 1, pre.end());
                    return { from_list(left_part), pre[i], reassemble(right_part, mid, suf) };
                }
            }
        }

        M pre_total = acc;
        for (const T& x : pre) pre_total = pre_total + Measured<T,M>::measure(x);

        // check spine — recurse into it (O(log n) descent)
        if (mid && !mid->empty() && pred(pre_total + mid->measure())) {
            // ask the spine tree to split itself at the right node
            auto node_split = mid->split([&](const M& m){ return pred(pre_total + m); });

            // the pivot from the spine split is a Node — find which element inside it
            const NodeT& pivot_node = node_split.pivot;
            M node_acc = pre_total + node_split.left.measure();
            size_t ei = 0;
            for (; ei < pivot_node.elems.size(); ei++) {
                node_acc = node_acc + Measured<T,M>::measure(pivot_node.elems[ei]);
                if (pred(node_acc)) break;
            }

            // build left: prefix + spine-left + node elements before ei
            FT left_tree = reassemble(pre, std::make_shared<SpineT>(node_split.left), {});
            for (size_t j = 0; j < ei; j++) left_tree = left_tree.push_back(pivot_node.elems[j]);

            // build right: node elements after ei + spine-right + suffix
            FT right_tree;
            for (size_t j = ei + 1; j < pivot_node.elems.size(); j++)
                right_tree = right_tree.push_back(pivot_node.elems[j]);
            right_tree = FT::concat(right_tree, reassemble({}, std::make_shared<SpineT>(node_split.right), suf));

            return { left_tree, pivot_node.elems[ei], right_tree };
        }

        // check suffix left to right
        {
            M run = pre_total;
            if (mid) run = run + mid->measure();
            for (size_t i = 0; i < suf.size(); i++) {
                run = run + Measured<T,M>::measure(suf[i]);
                if (pred(run)) {
                    std::vector<T> left_part(suf.begin(), suf.begin() + i);
                    std::vector<T> right_part(suf.begin() + i + 1, suf.end());
                    return { reassemble(pre, mid, left_part), suf[i], from_list(right_part) };
                }
            }
        }
        throw std::runtime_error("split: predicate never became true");
    }

    // ---- concat: join two trees   O(log n) ----
    static FT concat(const FT& a, const FT& b) { return join3(a, {}, b); }

    bool check_invariants() const {
        if (shape == Shape::Deep) {
            if (pre.size() < 1 || pre.size() > 4) return false;
            if (suf.size() < 1 || suf.size() > 4) return false;
            if (mid && !mid->check_invariants())   return false;
        }
        return true;
    }
};

// Public alias — depth 12 is enough for ~89 million characters
template <typename T, typename M>
using FingerTree = FT<T, M, 12>;
