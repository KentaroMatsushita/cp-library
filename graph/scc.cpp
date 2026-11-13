#pragma once
template <typename G = Graph> struct SCC {
    G g;
    Graph rg;
    vector<int> comp, ord;
    vector<bool> used;
    int num; // 連結成分の数

    SCC(G &g) : g(g), rg(g.size()), comp(g.size(), -1), used(g.size()) {
        rep(i, g.size()) for(auto &e : g[i]) rg[e].emplace_back(i);
        ord.reserve(g.size());
        build();
    };
    SCC(int n) : g(n), rg(n), comp(n, -1), used(n) { ord.reserve(n); };
    inline void add(int a, int b) {
        g[a].emplace_back(b);
        rg[b].emplace_back(a);
    }
    int operator[](int k) { return comp[k]; }
    void dfs(int x) {
        if(used[x]) return;
        used[x] = true;
        for(auto &e : g[x])
            if(!used[e]) dfs(e);
        ord.emplace_back(x);
    }
    void rdfs(int x, int cnt) {
        if(comp[x] != -1) return;
        comp[x] = cnt;
        for(int &e : rg[x])
            if(comp[e] == -1) rdfs(e, cnt);
    }
    void build() {
        fill(used.begin(), used.end(), false);
        fill(comp.begin(), comp.end(), -1);
        ord.clear();
        rep(i, g.size()) dfs(i);
        reverse(all(ord));
        num = 0;
        for(int &i : ord)
            if(comp[i] == -1) { rdfs(i, num), num++; }
    }
    // 強連結成分を潰した DAG を返す
    Graph getGraph() {
        Graph res(num);
        rep(i, g.size()) {
            for(auto &e : g[i]) {
                if(comp[e] == comp[i]) continue;
                res[comp[i]].emplace_back(comp[e]);
            }
        }
        rep(i, num) UNIQUE(res[i]);
        return res;
    }
    // 強連結成分ごとに属する頂点を返す
    vector<vector<int>> belong() {
        vector<vector<int>> res(num);
        rep(i, g.size()) res[comp[i]].emplace_back(i);
        return res;
    }
};
