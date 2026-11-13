// 重み付き木の直径とそのパスを求める
template <typename T> struct Diameter {
    Wgraph<T> &g;
    Diameter(Wgraph<T> &g) : g(g) {}
    vector<int> path;
    vector<int> to;

    pair<T, int> dfs(int x, int p) {
        pair<T, int> ret(0, x);
        for(auto &e : g[x]) {
            if(e != p) {
                auto cost = dfs(e, x);
                cost.first += e.cost;
                if(ret < cost) {
                    ret = cost;
                    to[x] = e.to;
                }
            }
        }
        return ret;
    }
    T build() {
        path.clear();
        if(g.empty()) return T{};
        to.assign(g.size(), -1);
        auto p = dfs(0, -1);
        auto q = dfs(p.second, -1);
        int now = p.second;
        while(now != q.second) {
            path.emplace_back(now);
            now = to[now];
        }
        path.emplace_back(now);
        return q.first;
    }
};
