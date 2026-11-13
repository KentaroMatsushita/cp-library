
struct TreeReduction {
    vector<int> root, sz, d, par;
    Graph g;
    struct UnionFind {
        using T = int;
        vector<T> data;
        vector<int> sz;
        int num;
        function<int(int)> f;
        UnionFind(int n, function<int(int)> f) : num(n), sz(n, -1), data(n), f(f) {}

        bool unite(int x, int y) {
            x = find(x), y = find(y);
            if(x == y) return false;
            num--;
            if(-sz[x] < -sz[y]) swap(x, y);
            if(f(data[x]) > f(data[y])) data[x] = data[y];
            sz[x] += sz[y];
            sz[y] = x;
            return true;
        }

        bool same(int x, int y) { return find(x) == find(y); }

        int find(int x) {
            if(sz[x] < 0) return x;
            return (sz[x] = find(sz[x]));
        }
        int size(int x) { return -sz[find(x)]; }
        T &get(int x) { return data[find(x)]; }
        const int operator[](const int k) { return find(k); }
        vector<vector<int>> belong() {
            vector<vector<int>> res(data.size());
            for(int i = 0; i < data.size(); i++) res[find(i)].emplace_back(i);
            return res;
        }
    };
    UnionFind uf;

    TreeReduction(int n) : g(n), root(n), sz(n, 1), d(n), par(n, -1), uf(n, [&](int x) { return d[x]; }) {
        rep(i, n) uf.get(i) = i;
        rep(i, n) root[i] = i;
    }
    void dfs(int x, int p) {
        sz[x] = 1;
        for(auto &&e : g[x]) {
            if(e == p) continue;
            root[e] = root[x], d[e] = d[x] + 1, par[e] = x;
            dfs(e, x);
            sz[x] += sz[e];
        }
    }
    TreeReduction(Graph g) : g(g), root(g.size(), -1), sz(g.size()), d(g.size()), par(g.size(), -1), uf(g.size(), [&](int x) { return d[x]; }) {
        rep(i, g.size()) uf.get(i) = i;
        rep(i, g.size()) {
            if(par[i] == -1) {
                root[i] = i;
                dfs(i, -1);
            }
        }
    }
    void add(int x, int y) {
        if(root[x] != root[y]) {
            g[x].emplace_back(y), g[y].emplace_back(x);
            if(sz[root[x]] < sz[root[y]]) swap(x, y);
            sz[root[x]] += sz[root[y]];
            queue<int> q;
            q.emplace(y);
            d[y] = d[x] + 1, par[y] = x;
            vi w;
            while(!empty(q)) {
                auto a = q.front();
                w.eb(a);
                q.pop();
                root[a] = root[x];
                for(auto &&e : g[a]) {
                    if(e != par[a]) {
                        q.emplace(e);
                        d[e] = d[a] + 1;
                        par[e] = a;
                    }
                }
            }
            for(auto &&e : w) {
                if(d[uf.get(e)] > d[e]) uf.get(e) = e;
            }
        } else {
            x = uf.get(x), y = uf.get(y);
            while(x != y) {
                if(d[x] < d[y]) swap(x, y);
                uf.unite(x, par[x]);
                x = uf.get(x);
            }
        }
    }
};
