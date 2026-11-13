// 木の直径とそのパスを求める
struct Diameter {
    Graph &g;
    Diameter(Graph &g) : g(g) {}
    vector<int> path;
    vector<int> to;

    pair<int, int> dfs(int x, int p) {
        pair<int, int> ret(0, x);
        for(auto &e : g[x]) {
            if(e != p) {
                auto cost = dfs(e, x);
                cost.first++;
                if(ret < cost) {
                    ret = cost;
                    to[x] = e;
                }
            }
        }
        return ret;
    }
    int build() {
        path.clear();
        if(g.empty()) return 0;
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
