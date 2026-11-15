// Standalone includes must precede using-directives and any personal template.
#include "tree/HLD.cpp"
#include "tree/tree_compressing.cpp"
#include "tree/tree_util.hpp"
#include "data_structure/Wavelet_Matrix.cpp"
#include <cassert>
#include <numeric>
#include <random>
#include <string>

struct AuditEdge {
    int to;
    long long cost;
    operator int() const { return to; }
};

static void audit_hld() {
    std::mt19937 rng(71237);
    for(int iteration = 0; iteration < 160; ++iteration) {
        int n = 1 + rng() % 70;
        std::vector<std::vector<int>> graph(n);
        for(int v = 1; v < n; ++v) {
            int p = rng() % v;
            graph[v].push_back(p);
            graph[p].push_back(v);
        }
        int root = rng() % n;
        HLDecomposition hld(graph, {root});
        hld.build();
        std::vector<int> parent(n, -1), depth(n), order{root};
        for(std::size_t i = 0; i < order.size(); ++i) {
            int v = order[i];
            for(int to : graph[v]) if(to != parent[v]) {
                parent[to] = v;
                depth[to] = depth[v] + 1;
                order.push_back(to);
            }
        }
        std::vector<int> values(n);
        std::iota(values.begin(), values.end(), 7);
        auto flattened = hld.build(values);
        for(int v = 0; v < n; ++v) {
            assert(hld.d[v] == depth[v]);
            assert(hld.par[v] == parent[v]);
            assert(hld.rev[hld.in[v]] == v);
            assert(flattened[hld.in[v]] == values[v]);
            assert(hld.la(v, depth[v] + 1) == -1);
            int u = v;
            for(int k = 0; k <= depth[v]; ++k, u = parent[u]) assert(hld.la(v, k) == u);
            for(int w = 0; w < n; ++w) {
                int ancestor = w;
                while(ancestor != -1 && ancestor != v) ancestor = parent[ancestor];
                assert((hld.in[v] <= hld.in[w] && hld.in[w] < hld.out[v]) == (ancestor == v));
            }
        }
        for(int query = 0; query < 150; ++query) {
            int a = rng() % n, b = rng() % n;
            int u = a, v = b;
            std::vector<int> left, right;
            while(depth[u] > depth[v]) { left.push_back(u); u = parent[u]; }
            while(depth[v] > depth[u]) { right.push_back(v); v = parent[v]; }
            while(u != v) { left.push_back(u); right.push_back(v); u = parent[u]; v = parent[v]; }
            left.push_back(u);
            std::reverse(right.begin(), right.end());
            left.insert(left.end(), right.begin(), right.end());
            assert(hld.lca(a, b) == u);
            assert(hld.road(a, b) == left);
            assert(hld.dist(a, b) == int(left.size()) - 1);
            assert(hld.jump(a, b, -1) == -1);
            for(int k = 0; k < int(left.size()); ++k) assert(hld.jump(a, b, k) == left[k]);
            assert(hld.jump(a, b, left.size()) == -1);
            std::vector<int> decoded;
            for(auto interval : hld.compressed_road(a, b)) {
                int l = interval.first, r = interval.second;
                for(int i = l;; i += l <= r ? 1 : -1) {
                    decoded.push_back(hld.rev[i]);
                    if(i == r) break;
                }
            }
            assert(decoded == left);
            auto query_forward = [&](int l, int r) {
                std::string answer;
                for(int i = l; i < r; ++i) answer.push_back(char(hld.rev[i]));
                return answer;
            };
            auto query_reverse = [&](int l, int r) {
                std::string answer;
                for(int i = r; i-- > l;) answer.push_back(char(hld.rev[i]));
                return answer;
            };
            std::string expected;
            for(int x : left) expected.push_back(char(x));
            assert(hld.query(a, b, std::string(), query_forward, query_reverse,
                [](const auto &x, const auto &y) { return x + y; }) == expected);
        }
    }
    std::vector<std::vector<AuditEdge>> weighted(5);
    for(int i = 1; i < 5; ++i) {
        weighted[i - 1].push_back({i, i * 13});
        weighted[i].push_back({i - 1, i * 13});
    }
    HLDecomposition weighted_hld(weighted, {2});
    auto weights = weighted_hld.build<long long>();
    for(int i = 0; i < 5; ++i) {
        long long expected = i == 2 ? 0 : (i < 2 ? (i + 1) * 13 : i * 13);
        assert(weights[weighted_hld.in[i]] == expected);
    }
    std::vector<std::vector<int>> forest{{1}, {0}, {}, {4}, {3}};
    HLDecomposition forest_hld(forest);
    forest_hld.build();
    auto permutation = forest_hld.in;
    std::sort(permutation.begin(), permutation.end());
    for(int i = 0; i < 5; ++i) assert(permutation[i] == i);
    std::vector<std::vector<int>> empty;
    HLDecomposition(empty).build();
    const int length = 300000;
    std::vector<std::vector<int>> chain(length);
    for(int i = 1; i < length; ++i) { chain[i - 1].push_back(i); chain[i].push_back(i - 1); }
    HLDecomposition chain_hld(chain);
    chain_hld.build();
    assert(chain_hld.lca(length - 1, length / 2) == length / 2);
    assert(chain_hld.la(length - 1, length - 1) == 0);
}

static void audit_wavelet_matrix() {
    std::mt19937 rng(82131);
    for(int iteration = 0; iteration < 120; ++iteration) {
        const int n = rng() % 85;
        std::vector<int> values(n);
        for(int &v : values) v = rng() % 128;
        WaveletMatrix<int, 7> wm(values);
        WaveletMatrix_sum<int, long long, 7> sums(n);
        for(int i = 0; i < n; ++i) sums.set(i, values[i]);
        sums.build();
        for(int repeat = 0; repeat < 2; ++repeat) {
            for(int i = 0; i < n; ++i) assert(wm.access(i) == values[i]);
            for(int query = 0; query < 200; ++query) {
                int l = rng() % (n + 1), r = rng() % (n + 1);
                if(l > r) std::swap(l, r);
                std::vector<int> ordered(values.begin() + l, values.begin() + r);
                std::sort(ordered.begin(), ordered.end());
                long long prefix = 0;
                for(int k = 0; k <= r - l; ++k) {
                    assert(sums.bottom_k_sum(l, r, k) == prefix);
                    assert(sums.top_k_sum(l, r, r - l - k) == std::accumulate(ordered.begin() + k, ordered.end(), 0LL));
                    if(k < r - l) { assert(wm.kth_smallest(l, r, k) == ordered[k]); prefix += ordered[k]; }
                }
                for(int threshold : {-5, 0, 1, 62, 127, 128, 1000}) {
                    assert(wm.range_freq(l, r, threshold) == std::lower_bound(ordered.begin(), ordered.end(), threshold) - ordered.begin());
                }
            }
            for(int &v : values) v = rng() % 128;
            for(int i = 0; i < n; ++i) { wm.set(i, values[i]); sums.set(i, values[i]); }
            wm.build();
            sums.build2();
        }
    }
    WaveletMatrix<int, 5> empty(0);
    empty.build();
    assert(empty.n == 0 && empty.range_freq(0, 0, 30) == 0);
    WaveletMatrix_sum<int, long long, 5> empty_sum;
    empty_sum.build();
    assert(empty_sum.top_k_sum(0, 0, 0) == 0);
    WaveletMatrixBitVector bits(128);
    bits.set(0); bits.set(63); bits.set(64); bits.set(127); bits.build();
    assert(bits.rank1(128) == 4 && bits.rank1(64) == 2 && bits.rank1(0) == 0);
    bits.init(128); bits.build();
    assert(bits.rank1(128) == 0);
    using U64 = std::uint64_t;
    std::vector<U64> extremes{0, U64(1) << 63, ~U64(0), 4};
    WaveletMatrix<U64, 64> wide(extremes);
    assert(wide.kth_smallest(0, 4, 3) == ~U64(0));
    assert(wide.range_freq(0, 4, ~U64(0)) == 3);
    assert(wide.access(1) == (U64(1) << 63));
    bool threw = false;
    try { WaveletMatrix<int, 4> invalid(std::vector<int>{16}); } catch(const std::out_of_range &) { threw = true; }
    assert(threw);
    threw = false;
    try { WaveletMatrix<int, 4> invalid(std::vector<int>{-1}); } catch(const std::invalid_argument &) { threw = true; }
    assert(threw);
}

static void audit_virtual_tree() {
    std::mt19937 rng(172196);
    for(int iteration = 0; iteration < 180; ++iteration) {
        const int n = 1 + rng() % 70;
        std::vector<std::vector<int>> graph(n);
        for(int i = 1; i < n; ++i) if(rng() % 4 != 0) {
            int p = rng() % i;
            graph[i].push_back(p);
            graph[p].push_back(i);
        }
        LCA_Tree tree(graph);
        std::vector<int> selected;
        for(int i = 0; i < n; ++i) if(rng() % 3 == 0) selected.push_back(i);
        std::vector<int> required(n);
        for(int v : selected) { required[v] = 1; required[tree.hld.root[v]] = 1; }
        for(int u : selected) for(int v : selected) {
            if(tree.hld.root[u] == tree.hld.root[v]) required[tree.lca(u, v)] = 1;
        }
        if(!selected.empty()) selected.push_back(selected.front());
        auto [vertices, compressed] = tree.make(selected);
        assert(vertices == selected);
        std::vector<int> position(n, -1);
        for(int i = 0; i < int(vertices.size()); ++i) {
            position[vertices[i]] = i;
            assert(required[vertices[i]]);
        }
        for(int i = 0; i < n; ++i) assert((position[i] >= 0) == bool(required[i]));
        std::size_t expected_edges = 0;
        for(int i = 0; i < int(vertices.size()); ++i) {
            int p = tree.hld.par[vertices[i]];
            while(p != -1 && position[p] == -1) p = tree.hld.par[p];
            if(p != -1) {
                assert(std::count(compressed[i].begin(), compressed[i].end(), position[p]) == 1);
                assert(std::count(compressed[position[p]].begin(), compressed[position[p]].end(), i) == 1);
                expected_edges += 2;
            }
        }
        std::size_t actual_edges = 0;
        for(auto &neighbors : compressed) actual_edges += neighbors.size();
        assert(actual_edges == expected_edges);
        TreeUtil utilities(graph);
        for(int i = 0; i < n; ++i) {
            assert(utilities.subtree_size(i) == tree.hld.sz[i]);
            assert(utilities.depth(i) == tree.hld.d[i]);
            assert(utilities.rev(utilities.lid(i)) == i);
        }
    }
}

int main() {
    audit_virtual_tree();
    audit_hld();
    audit_wavelet_matrix();
}
