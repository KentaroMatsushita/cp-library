template <typename key_t, typename sum_t = key_t> struct ReRooting {
    struct Edge {
        int from, to;
        key_t cost;
        sum_t dp, ndp;
    };
    using F = function<sum_t(sum_t, sum_t)>;
    using G = function<sum_t(sum_t, Edge)>;

    vector<vector<Edge>> g;
    const F f;
    const G gg;
    const sum_t id;
    vector<sum_t> subdp, dp;

    void add_edge(int u, int v, const key_t &d) {
        g[u].emplace_back(Edge{u, v, d, id, id});
        g[v].emplace_back(Edge{v, u, d, id, id});
    }
    void add_edge(int u, int v, const key_t &d, const key_t &rev_d) {
        g[u].emplace_back(Edge{u, v, d, id, id});
        g[v].emplace_back(Edge{v, u, rev_d, id, id});
    }
    void add_directed_edge(int u, int v, const key_t &d) { g[u].emplace_back(Edge{u, v, d, id, id}); }

    ReRooting(int n, const G &g, const F &f, const sum_t &id = sum_t{}) : g(n), f(f), gg(g), id(id), subdp(n, id), dp(n, id) {}
    ReRooting(Graph &graph, const G &g, const F &f, const sum_t &id = sum_t{})
        : g(graph.size()), f(f), gg(g), id(id), subdp(graph.size(), id), dp(graph.size(), id) {
        for(int i = 0; i < graph.size(); i++)
            for(auto e : graph[i]) add_directed_edge(i, e, {});
    }
    template <typename Weight>
    ReRooting(Wgraph<Weight> &graph, const G &g, const F &f, const sum_t &id = sum_t{})
        : g(graph.size()), f(f), gg(g), id(id), subdp(graph.size(), id), dp(graph.size(), id) {
        for(int i = 0; i < graph.size(); i++)
            for(auto e : graph[i]) add_directed_edge(i, e, e.cost);
    }

    void dfs_sub(int x, int p) {
        for(auto &e : g[x]) {
            if(e.to == p) continue;
            dfs_sub(e.to, x);
            e.dp = gg(subdp[e.to], e);
            subdp[x] = f(subdp[x], e.dp);
        }
    }

    void dfs_all(int x, int p, sum_t top) {
        sum_t now = id;
        for(int i = 0; i < (int)g[x].size(); i++) {
            auto &e = g[x][i];
            e.ndp = now;
            if(e.to == p) e.dp = gg(top, e);
            now = f(now, e.dp);
        }
        dp[x] = now;
        now = id;
        for(int i = (int)g[x].size() - 1; i >= 0; i--) {
            auto &e = g[x][i];
            if(e.to != p) dfs_all(e.to, x, f(e.ndp, now));
            now = f(now, e.dp);
        }
    }

    vector<sum_t> build() {
        fill(subdp.begin(), subdp.end(), id);
        fill(dp.begin(), dp.end(), id);
        if(g.empty()) return dp;
        dfs_sub(0, -1);
        dfs_all(0, -1, id);
        return dp;
    }
};

/*
key_t : 辺に持たせるもの
sum_t : DP の結果
ff: sum_t 同士をまとめる
gg: 頂点での DP の結果を、辺を伝って上に伝える(from が自分、to が伝える先)

ex. 最遠頂点と距離

auto ff = [&](pii x, pii y) {
    return (x.fi > y.fi ? x : y)
};
auto gg = [&](pii x, auto e) {
    if(x.fi == 0) x.se = e.from;
    return pii(x.fi + e.cost, x.se);
};
ReRooting<int, pii> rr(g, gg, ff, pii());
auto res = rr.build();

*/
