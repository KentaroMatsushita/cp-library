#pragma once
#include <algorithm>
#include <queue>
#include <utility>
#include <vector>
template <typename G> struct HLDecomposition {
    G &g;
    std::vector<int> sz, in, out, head, rev, par, d, vis, roots, root;

    HLDecomposition(G &g, std::vector<int> r = std::vector<int>())
        : g(g), sz(g.size()), in(g.size()), out(g.size()), head(g.size()), rev(g.size()), par(g.size()), d(g.size()), root(g.size()) {
        if(std::empty(r)) {
            std::vector<bool> vis(g.size());
            for(int i = 0; i < g.size(); i++) {
                if(vis[i]) continue;
                roots.emplace_back(i);
                std::queue<int> q;
                q.emplace(i);
                while(!std::empty(q)) {
                    int x = q.front();
                    vis[x] = true;
                    q.pop();
                    for(auto e : g[x]) {
                        if(!vis[e]) q.emplace(e);
                    }
                }
            }
        } else
            roots = r;
    }

    // Iterative traversals keep construction safe for path-shaped trees.
    void dfs_sz(int idx, int p) {
        std::vector<int> order{idx};
        root[idx] = p < 0 ? idx : root[p];
        par[idx] = p;
        d[idx] = p < 0 ? 0 : d[p] + 1;
        for(std::size_t i = 0; i < order.size(); ++i) {
            int v = order[i];
            for(auto &to : g[v]) {
                if(to == par[v]) continue;
                par[to] = v;
                root[to] = root[v];
                d[to] = d[v] + 1;
                order.push_back(to);
            }
        }
        for(auto it = order.rbegin(); it != order.rend(); ++it) {
            int v = *it;
            sz[v] = 1;
            int heavy = -1;
            for(int i = 0; i < int(g[v].size()); ++i) {
                int to = g[v][i];
                if(to == par[v]) continue;
                sz[v] += sz[to];
                if(heavy < 0 || sz[to] > sz[int(g[v][heavy])]) heavy = i;
            }
            if(heavy >= 0) std::swap(g[v][0], g[v][heavy]);
        }
    }

    template <typename Visit> void traverse_hld(int idx, int &times, const Visit &visit) {
        std::vector<int> stack{idx};
        while(!stack.empty()) {
            int v = stack.back();
            stack.pop_back();
            in[v] = times++;
            rev[in[v]] = v;
            out[v] = in[v] + sz[v];
            visit(v);
            for(int i = int(g[v].size()); i-- > 0;) {
                int to = g[v][i];
                if(to == par[v]) continue;
                head[to] = i == 0 ? head[v] : to;
                stack.push_back(to);
            }
        }
    }

    void dfs_hld(int idx, int, int &times) {
        traverse_hld(idx, times, [](int) {});
    }

    template <typename T> void dfs_hld(int idx, int, int &times, std::vector<T> &v) {
        traverse_hld(idx, times, [&](int x) {
            for(auto &to : g[x]) if(to == par[x]) { v[in[x]] = to.cost; break; }
        });
    }

    template <typename T> void dfs_hld(int idx, int, int &times, std::vector<T> &v, std::vector<T> &a) {
        traverse_hld(idx, times, [&](int x) { v[in[x]] = a[x]; });
    }

    void build() {
        int t = 0;
        for(auto root : roots) {
            head[root] = root;
            dfs_sz(root, -1);
            dfs_hld(root, -1, t);
        }
    }

    template <typename T> std::vector<T> build(int /* root is selected in the constructor */ = 0) {
        std::vector<T> res(g.size());
        int t = 0;
        for(auto root : roots) {
            head[root] = root;
            dfs_sz(root, -1);
            dfs_hld(root, -1, t, res);
        }
        return res;
    }

    template <typename T> std::vector<T> build(std::vector<T> &a, int /* root is selected in the constructor */ = 0) {
        std::vector<T> res(g.size());
        int t = 0;
        for(auto root : roots) {
            head[root] = root;
            dfs_sz(root, -1);
            dfs_hld(root, -1, t, res, a);
        }
        return res;
    }

    int la(int v, int k) {
        if(k < 0 || k > d[v]) return -1;
        while(1) {
            int u = head[v];
            if(in[v] - k >= in[u]) return rev[in[v] - k];
            k -= in[v] - in[u] + 1;
            v = par[u];
        }
    }

    int lca(int u, int v) {
        for(;; v = par[head[v]]) {
            if(in[u] > in[v]) std::swap(u, v);
            if(head[u] == head[v]) return u;
        }
    }

    template <typename T, typename Q, typename F> T query(int u, int v, const T &e, const Q &q, const F &f, bool edge = false) {
        T l = e, r = e;
        for(;; v = par[head[v]]) {
            if(in[u] > in[v]) std::swap(u, v), std::swap(l, r);
            if(head[u] == head[v]) break;
            l = f(q(in[head[v]], in[v] + 1), l);
        }
        return f(f(q(in[u] + edge, in[v] + 1), l), r);
    }

    template <typename T, typename Q, typename Q2, typename F> T query(int u, int v, const T &e, const Q &q1, const Q2 &q2, const F &f, bool edge = false) {
        T l = e, r = e;
        for(;;) {
            if(head[u] == head[v]) break;
            if(in[u] > in[v]) {
                l = f(l, q2(in[head[u]], in[u] + 1));
                u = par[head[u]];
            } else {
                r = f(q1(in[head[v]], in[v] + 1), r);
                v = par[head[v]];
            }
        }
        if(in[u] > in[v]) return f(f(l, q2(in[v] + edge, in[u] + 1)), r);
        return f(f(l, q1(in[u] + edge, in[v] + 1)), r);
    }

    template <typename Q> void add(int u, int v, const Q &q, bool edge = false) {
        for(;; v = par[head[v]]) {
            if(in[u] > in[v]) std::swap(u, v);
            if(head[u] == head[v]) break;
            q(in[head[v]], in[v] + 1);
        }
        q(in[u] + edge, in[v] + 1);
    }

    template <typename Q1, typename Q2> void add(int u, int v, const Q1 &q1, const Q2 &q2, bool edge = false) {
        for(;;) {
            if(head[u] == head[v]) break;
            if(in[u] > in[v]) {
                q2(in[head[u]], in[u] + 1);
                u = par[head[u]];
            } else {
                q1(in[head[v]], in[v] + 1);
                v = par[head[v]];
            }
        }
        if(in[u] > in[v])
            q2(in[v] + edge, in[u] + 1);
        else
            q1(in[u] + edge, in[v] + 1);
    }

    constexpr int operator[](int k) { return in[k]; }

    constexpr int dist(int u, int v) { return d[u] + d[v] - 2 * d[lca(u, v)]; }

    // u -> v の unique path
    std::vector<int> road(int u, int v) {
        int l = lca(u, v);
        std::vector<int> a, b;
        for(; v != l; v = par[v]) b.emplace_back(v);
        for(; u != l; u = par[u]) a.emplace_back(u);
        a.emplace_back(l);
        for(auto it = b.rbegin(); it != b.rend(); ++it) a.emplace_back(*it);
        return a;
    }

    // u -> v の unique path の in[x] を連続部分ごとに固めたもの
    // l > r なら降順
    std::vector<std::pair<int, int>> compressed_road(int u, int v) {
        std::vector<std::pair<int, int>> l, r, res;
        add(u, v, [&](int x, int y) { return r.emplace_back(x, y - 1); }, [&](int x, int y) { l.emplace_back(x, y - 1); });
        for(auto [x, y] : l) res.emplace_back(y, x);
        std::reverse(std::begin(r), std::end(r));
        for(auto [x, y] : r) res.emplace_back(x, y);
        return res;
    }
    int jump(int s, int t, int i) {
        if(i < 0) return -1;
        if(!i) return s;
        int l = lca(s, t);
        int dst = d[s] + d[t] - d[l] * 2;
        if(dst < i) return -1;
        if(d[s] - d[l] >= i) return la(s, i);
        i -= d[s] - d[l];
        return la(t, d[t] - d[l] - i);
    }
};
