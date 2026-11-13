
template <typename G> struct HLDecomposition {
    G g;
    vector<int> sz, in, out, head, rev, par, d, vis, roots;

    HLDecomposition() = default;
    HLDecomposition(G &g, vector<int> r = vector<int>())
        : g(g), d(g.size()), sz(g.size()), in(g.size()), out(g.size()), head(g.size()), rev(g.size()), par(g.size()) {
        if(empty(r)) {
            vector<bool> vis(g.size());
            for(int i = 0; i < g.size(); i++) {
                if(vis[i]) continue;
                roots.emplace_back(i);
                queue<int> q;
                q.emplace(i);
                while(!empty(q)) {
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

    void dfs_sz(int idx, int p) {
        par[idx] = p;
        sz[idx] = 1;
        if(g[idx].size() and g[idx][0] == p) swap(g[idx][0], g[idx].back());
        for(auto &to : g[idx]) {
            if(to == p) continue;
            d[to] = d[idx] + 1;
            dfs_sz(to, idx);
            sz[idx] += sz[to];
            if(sz[g[idx][0]] < sz[to]) swap(g[idx][0], to);
        }
    }

    void dfs_hld(int idx, int par, int &times) {
        in[idx] = times++;
        rev[in[idx]] = idx;
        for(auto &to : g[idx]) {
            if(to == par) continue;
            head[to] = (g[idx][0] == to ? head[idx] : to);
            dfs_hld(to, idx, times);
        }
        out[idx] = times;
    }

    template <typename T> void dfs_hld(int idx, int par, int &times, vector<T> &v) {
        in[idx] = times++;
        rev[in[idx]] = idx;
        for(auto &to : g[idx]) {
            if(to == par) {
                v[in[idx]] = to.cost;
                continue;
            }
            head[to] = (g[idx][0] == to ? head[idx] : to);
            dfs_hld(to, idx, times, v);
        }
        out[idx] = times;
    }

    template <typename T> void dfs_hld(int idx, int par, int &times, vector<T> &v, vector<T> &a) {
        in[idx] = times++;
        rev[in[idx]] = idx;
        v[in[idx]] = a[idx];
        for(auto &to : g[idx]) {
            if(to == par) continue;
            head[to] = (g[idx][0] == to ? head[idx] : to);
            dfs_hld(to, idx, times, v, a);
        }
        out[idx] = times;
    }

    void build() {
        int t = 0;
        for(auto root : roots) {
            head[root] = root;
            dfs_sz(root, -1);
            dfs_hld(root, -1, t);
        }
    }

    template <typename T> vector<T> build(int root = 0) {
        vector<T> res(g.size());
        int t = 0;
        for(auto root : roots) {
            head[root] = root;
            dfs_sz(root, -1);
            dfs_hld(root, -1, t, res);
        }
        return res;
    }

    template <typename T> vector<T> build(vector<T> &a, int root = 0) {
        vector<T> res(g.size());
        for(auto root : roots) {
            head[root] = root;
            dfs_sz(root, -1);
            int t = 0;
            dfs_hld(root, -1, t, res, a);
        }
        return res;
    }

    int la(int v, int k) {
        while(1) {
            int u = head[v];
            if(in[v] - k >= in[u]) return rev[in[v] - k];
            k -= in[v] - in[u] + 1;
            v = par[u];
        }
    }

    int lca(int u, int v) {
        for(;; v = par[head[v]]) {
            if(in[u] > in[v]) swap(u, v);
            if(head[u] == head[v]) return u;
        }
    }

    template <typename T, typename Q, typename F> T query(int u, int v, const T &e, const Q &q, const F &f, bool edge = false) {
        T l = e, r = e;
        for(;; v = par[head[v]]) {
            if(in[u] > in[v]) swap(u, v), swap(l, r);
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
            if(in[u] > in[v]) swap(u, v);
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
    vector<int> road(int u, int v) {
        int l = lca(u, v);
        vector<int> a, b;
        for(; v != l; v = par[v]) b.eb(v);
        for(; u != l; u = par[u]) a.eb(u);
        a.eb(l);
        per(i, si(b), 0) a.eb(b[i]);
        return a;
    }

    // u -> v の unique path の in[x] を連続部分ごとに固めたもの
    // l > r なら降順
    vector<pair<int, int>> compressed_road(int u, int v) {
        vector<pair<int, int>> l, r, res;
        add(u, v, [&](int x, int y) { return r.emplace_back(x, y - 1); }, [&](int x, int y) { l.emplace_back(x, y - 1); });
        for(auto [x, y] : l) res.emplace_back(y, x);
        reverse(begin(r), end(r));
        for(auto [x, y] : r) res.emplace_back(x, y);
        return res;
    }
    int jump(int s, int t, int i) {
        if(!i) return s;
        int l = lca(s, t);
        int dst = d[s] + d[t] - d[l] * 2;
        if(dst < i) return -1;
        if(d[s] - d[l] >= i) return la(s, i);
        i -= d[s] - d[l];
        return la(t, d[t] - d[l] - i);
    }
};

template <typename G> struct LowLink {
    const G &g;
    vector<int> used, ord, low;
    vector<int> articulation;
    vector<pair<int, int>> bridge;

    LowLink(const G &g) : g(g) {}

    int dfs(int idx, int k, int par) {
        used[idx] = true;
        ord[idx] = k++;
        low[idx] = ord[idx];
        bool is_articulation = false;
        int cnt = 0;
        for(auto &to : g[idx]) {
            if(!used[to]) {
                ++cnt;
                k = dfs(to, k, idx);
                low[idx] = min(low[idx], low[to]);
                is_articulation |= ~par && low[to] >= ord[idx];
                if(ord[idx] < low[to]) bridge.emplace_back(minmax(idx, (int)to));
            } else if(to != par) {
                low[idx] = min(low[idx], ord[to]);
            }
        }
        is_articulation |= par == -1 && cnt > 1;
        if(is_articulation) articulation.push_back(idx);
        return k;
    }

    virtual void build() {
        articulation.clear(); bridge.clear();
        used.assign(g.size(), 0);
        ord.assign(g.size(), 0);
        low.assign(g.size(), 0);
        int k = 0;
        for(int i = 0; i < g.size(); i++) {
            if(!used[i]) k = dfs(i, k, -1);
        }
    }
};
template <typename G> struct BiConnectedComponents : LowLink<G> {
    using LoL = LowLink<G>;

    vector<int> used;
    vector<vector<pair<int, int>>> bc;
    vector<pair<int, int>> tmp;

    BiConnectedComponents(const G &g) : LoL(g) {}

    void dfs(int idx, int par) {
        used[idx] = true;
        for(auto &to : this->g[idx]) {
            if(to == par) continue;
            if(!used[to] || this->ord[to] < this->ord[idx]) { tmp.emplace_back(minmax(idx, to)); }
            if(!used[to]) {
                dfs(to, idx);
                if(this->low[to] >= this->ord[idx]) {
                    bc.emplace_back();
                    for(;;) {
                        auto e = tmp.back();
                        bc.back().emplace_back(e);
                        tmp.pop_back();
                        if(e.first == min(idx, to) && e.second == max(idx, to)) { break; }
                    }
                }
            }
        }
    }

    void build() override {
        bc.clear(); tmp.clear();
        LoL::build();
        used.assign(this->g.size(), 0);
        for(int i = 0; i < used.size(); i++) {
            if(!used[i]) dfs(i, -1);
        }
    }
};

// aux : block cut tree
// is_arti(V) : return if V is articulation
// 0, ..., N - 1 : 元の頂点に対応する点(間接点)
// N, ..., N + si(bcc) - 1 : bcc
template <typename G> struct BlockCutTree {
    const G &g;
    BiConnectedComponents<G> bcc;
    vector<vector<int>> aux;
    vector<int> idar, idcc, arti;

    HLDecomposition<vector<vector<int>>> hld;

    BlockCutTree(const G &_g) : g(_g), bcc(_g) {
        bcc.build();
        build();
    }

    void build() {
        auto ar = bcc.articulation;
        arti.assign(g.size(),0);
        for(int i = 0; i < (int)ar.size(); i++) arti[ar[i]] = true;
        aux.assign(g.size() + bcc.bc.size(),{});

        vector<int> used(g.size());
        for(int i = 0; i < (int)bcc.bc.size(); i++) {
            for(auto &[u, v] : bcc.bc[i]) {
                if(!exchange(used[u], 1)) { add(u, g.size() + i); }
                if(!exchange(used[v], 1)) { add(v, g.size() + i); }
            }
            for(auto &[u, v] : bcc.bc[i]) used[u] = used[v] = 0;
        }
        hld = HLDecomposition<decltype(aux)>(aux);
        hld.build();
    }

    vector<int> &operator[](int i) { return aux[i]; }

    int size() const { return aux.size(); }

    bool is_arti(int i) { return arti[i]; }

    int arti_num() const { return bcc.articulation.size(); }

    // bcc の番号から aux での番号
    int bcc_id(int x) { return x + g.size(); }

    // 辺 -> その辺が所属する BCC
    int edge_to_bcc(int u, int v) { return hld.jump(u, v, 1) - (int)g.size(); }

    // a - b - c の順に現れる単純パスが存在するか
    int exist_path(int a, int b, int c) { return hld.dist(a, b) + hld.dist(b, c) - 2 <= hld.dist(a, c); }

  private:
    void add(int i, int j) {
        if(i == -1 or j == -1) return;
        aux[i].push_back(j);
        aux[j].push_back(i);
    };
};

/**
 * @brief Block Cut Tree
 */
