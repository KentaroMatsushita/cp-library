#pragma once
#include "HLD.cpp"
// Edge type for utilities requiring endpoint conversion and explicit edge ids.
template <typename T> struct Edge {
    int from = -1, to = -1;
    T cost{};
    int id = -1;
    operator int() const { return to; }
};
// Build HLD once and expose the traditional utility getter names.
template <typename G> struct TreeUtil : HLDecomposition<G> {
    using Base = HLDecomposition<G>;
    TreeUtil(G& g, std::vector<int> roots = {}) : Base(g, roots) { this->build(); }
    int subtree_size(int v) const { return Base::sz[v]; }
    int depth(int v) const { return Base::d[v]; }
    int lid(int v) const { return Base::in[v]; }
    int rid(int v) const { return Base::out[v]; }
    int rev(int i) const { return Base::rev[i]; }
};