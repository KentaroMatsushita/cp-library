// 連結成分ごとに (黒, 白) を分ける
template <typename G> vector<pair<vector<int>, vector<int>>> bipartite_decompose(const G &g) {
    int n = g.size();
    vector<int> vis(n), col(n);
    vector<pair<vector<int>, vector<int>>> res;
    for(int i = 0; i < n; i++) {
        if(!vis[i]) {
            vector<int> b, w;
            queue<int> q({i});
            vis[i] = true;
            while(!q.empty()) {
                int x = q.front();
                (col[x] ? w : b).emplace_back(x);
                q.pop();
                for(auto &&e : g[x]) {
                    if(!vis[e]) {
                        vis[e] = true;
                        col[e] = col[x] ^ 1;
                        q.emplace(e);
                    } else if(col[e] == col[x]) {
                        return {};
                    }
                }
            }
            res.emplace_back(b, w);
        }
    }
    return res;
}