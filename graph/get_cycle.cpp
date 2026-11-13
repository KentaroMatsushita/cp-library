// 有向グラフ
struct DirectedCycle {
    Graph &g;
    vector<int> vis;
    vector<int> res;
    DirectedCycle(Graph &g) : g(g), vis(g.size()) {}
    int start;
    bool dfs(int x) {
        vis[x] = 1;
        for(auto e : g[x]) {
            if(vis[e] == 1) {
                res.emplace_back(x);
                start = e;
                return start != x;
            } else if(vis[e] == 2)
                continue;
            if(dfs(e)) {
                res.emplace_back(x);
                if(start == x) return false;
                return true;
            } else if(!res.empty())
                return false;
        }
        vis[x] = 2;
        return false;
    }
    vector<int> build() {
        fill(vis.begin(),vis.end(),0);
        res.clear();
        start = -1;
        rep(i, g.size()) {
            if(!vis[i]) {
                dfs(i);
                if(!res.empty()) {
                    reverse(res.begin(), res.end());
                    return res;
                }
            }
        }
        return vector<int>();
    }
};

// 無向グラフ
struct UndirectedCycle {
    Graph &g;
    vector<int> vis;
    vector<int> res;
    UndirectedCycle(Graph &g) : g(g), vis(g.size()) {}
    int start;
    bool dfs(int x, int p = -1) {
        vis[x] = 1;
        bool skipped_parent = false;
        for(auto e : g[x]) {
            if(e == p && !exchange(skipped_parent,true)) continue;
            if(vis[e] == 1) {
                res.emplace_back(x);
                start = e;
                return start != x;
            } else if(vis[e] == 2)
                continue;
            if(dfs(e, x)) {
                res.emplace_back(x);
                if(start == x) return false;
                return true;
            } else if(!res.empty())
                return false;
        }
        vis[x] = 2;
        return false;
    }
    vector<int> build() {
        fill(vis.begin(),vis.end(),0);
        res.clear();
        start = -1;
        rep(i, g.size()) {
            if(!vis[i]) {
                dfs(i);
                if(!res.empty()) {
                    reverse(res.begin(), res.end());
                    return res;
                }
            }
        }
        return vector<int>();
    }
};
using GetCycle = DirectedCycle;
