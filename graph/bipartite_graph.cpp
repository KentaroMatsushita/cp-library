// 二部グラフでないときは空の vector を返す
template <typename G> vector<bool> color_two(G &g) {
    int n = g.size();
    vector<bool> used(n), res(n);
    queue<int> q;
    for(int i = 0; i < n; ++i) {
        if(!used[i]) used[i] = true, q.emplace(i);
        while(!q.empty()) {
            int now = q.front();
            q.pop();
            for(auto &e : g[now]) {
                if(!used[e])
                    used[e] = true, res[e] = res[now] ^ 1, q.emplace(e);
                else if(res[e] == res[now])
                    return vector<bool>();
            }
        }
    }
    return res;
}