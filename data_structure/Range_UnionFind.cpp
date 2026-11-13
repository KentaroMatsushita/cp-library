namespace RangeUF_impl {

struct UnionFind {
    vector<int> data;
    UnionFind(int n) { data.assign(n, -1); }

    bool unite(int x, int y) {
        x = find(x), y = find(y);
        if(x == y) return false;
        if(data[x] > data[y]) swap(x, y);
        data[x] += data[y];
        data[y] = x;
        return true;
    }

    bool same(int x, int y) { return find(x) == find(y); }

    vector<int> _v;
    int find(int x) {
        if(data[x] < 0) return x;
        _v.clear();
        while(data[x] >= 0) {
            _v.emplace_back(x);
            x = data[x];
        }
        for(auto e : _v) data[e] = x;
        return x;
    }
    int size(int x) { return -data[find(x)]; }
    const int operator[](const int k) { return find(k); }
    vector<vector<int>> belong() {
        vector<vector<int>> res(data.size());
        for(int i = 0; i < data.size(); i++) res[find(i)].emplace_back(i);
        return res;
    }
    void reset() { data.assign(si(data), -1); }
};

struct RangeUF {
    int n, lg, num;
    UnionFind uf;
    int log_floor(int n) {
        int res = 0;
        while((1 << res) * 2 <= n) res++;
        return res;
    }
    RangeUF(int n) : n(n), num(n), lg(log_floor(n)), uf((log_floor(n) + 1) * n) {}

  private:
    queue<tuple<int, int, int>> q;
    constexpr int id(int t, int x) { return t * n + x; }

    template <typename F> bool _unite(int u, int v, int l, const F &f) {
        if(l == 0) return false;
        int t = topbit(l);
        if(uf.unite(id(t, u), id(t, v))) q.emplace(u, v, t);
        if(l > (1 << t) and uf.unite(id(t, u + l - (1 << t)), id(t, v + l - (1 << t)))) q.emplace(u + l - (1 << t), v + l - (1 << t), t);
        bool flag = false;
        while(!empty(q)) {
            auto [u, v, t] = q.front();
            q.pop();
            if(t) {
                if(uf.unite(id(t - 1, u), id(t - 1, v))) q.emplace(u, v, t - 1);
                int nu = u + (1 << t - 1), nv = v + (1 << t - 1);
                if(uf.unite(id(t - 1, nu), id(t - 1, nv))) q.emplace(nu, nv, t - 1);
            }
            if(t == 0) {
                flag = true;
                num--;
                f(u, v);
            }
        }
        return flag;
    }

  public:
    // unite(u, v), unite(u + 1, v + 1), ..., unite(u + l - 1, v + l - 1)
    bool unite(int u, int v, int len) {
        return _unite(u, v, len, [](int, int) {});
    }
    template <typename F> bool unite(int u, int v, int len, const F &f) { return _unite(u, v, len, f); }

    bool same(int u, int v) { return uf.same(u, v); }
    int size(int u) { return uf.size(u); }
    vector<vector<int>> belong() {
        vector<vector<int>> res(n);
        rep(i, n) res[uf[i]].emplace_back(i);
        return res;
    }
    int component_number() { return num; }
    int par(int u) { return uf[u]; }
};

struct RangeUF2Way {
    int n, lg;
    // vector<UnionFind> uf;
    int num;
    UnionFind uf;
    int log_floor(int n) {
        int res = 0;
        while((1 << res) * 2 <= n) res++;
        return res;
    }
    RangeUF2Way(int _n) : n(_n), num(n), lg(log_floor(_n)), uf((log_floor(_n) + 1) * (_n * 2) - _n) {}

  private:
    constexpr int id(int t, int x) {
        if(t == 0) return (x < n ? x : rev(x));
        return t * (2 * n) - n + x;
    }
    queue<tuple<int, int, int>> q;
    constexpr int rev(int u) { return 2 * n - 1 - u; }

    template <typename F> bool _unite(int u, int v, int l, const F &f) {
        if(l == 0) return false;
        int t = topbit(l);
        if(uf.unite(id(t, u), id(t, v))) q.emplace(u, v, t);
        if(l > (1 << t) and uf.unite(id(t, u + l - (1 << t)), id(t, v + l - (1 << t)))) q.emplace(u + l - (1 << t), v + l - (1 << t), t);
        bool flag = false;
        while(!empty(q)) {
            auto [u, v, t] = q.front();
            q.pop();
            if(t == 0) {
                if(u >= n) u = rev(u);
                if(v >= n) v = rev(v);
            }
            if(t) {
                if(uf.unite(id(t - 1, u), id(t - 1, v))) q.emplace(u, v, t - 1);
                int nu = u + (1 << t - 1), nv = v + (1 << t - 1);
                if(uf.unite(id(t - 1, nu), id(t - 1, nv))) q.emplace(nu, nv, t - 1);
            }
            if(t == 0) {
                flag = true;
                num--;
                f(u, v);
            }
        }
        return flag;
    }

  public:
    // unite(u, v), unite(u + udir, v + vdir), ..., unite(u + (l - 1) * udir, v + (l - 1) * vdir)
    bool unite(int u, int udir, int v, int vdir, int len) {
        return _unite(udir == 1 ? u : rev(u), vdir == 1 ? v : rev(v), len, [](int, int) {});
    }
    template <typename F> bool unite(int u, int udir, int v, int vdir, int len, const F &f) {
        return _unite(udir == 1 ? u : rev(u), vdir == 1 ? v : rev(v), len, f);
    }

    bool same(int u, int v) { return uf.same(u, v); }
    int size(int u) { return uf.size(u); }
    vector<vector<int>> belong() {
        vector<vector<int>> res(n);
        for(int i = 0; i < n; i++) {
            int x = uf[i];
            if(!(x < n)) x = rev(x);
            res[x].emplace_back(i);
        }
        return res;
    }
    int component_number() { return num; }
    int par(int u) {
        int x = uf[u];
        if(x >= n) x = rev(x);
        return x;
    }
};
} // namespace RangeUF_impl
using RangeUF_impl::RangeUF;
using RangeUF_impl::RangeUF2Way;
