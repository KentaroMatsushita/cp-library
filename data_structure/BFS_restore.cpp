vector<int> BFS(const Graph &g, int s, int t) {
    int n = g.size();
    vector<int> d(n, INT_MAX), pre(n);
    d[s] = 0;
    queue<int> q({s});
    while(!q.empty()) {
        int a = q.front();
        q.pop();
        for(auto e : g[a]) {
            if(d[e] > d[a] + 1) {
                d[e] = d[a] + 1;
                pre[e] = a;
                q.emplace(e);
            }
        }
    }
    if(d[t] == INT_MAX) return {};
    vector<int> res;
    for(int now=t;;now=pre[now]) {
        res.push_back(now);
        if(now==s) break;
    }
    reverse(res.begin(),res.end());
    return res;
}
