pair<int, vector<int>> BFS(const Wgraph<int> &g, int s, int t) {
    int n = g.size();
    vector<int> d(n, INT_MAX), pre(n), used(n);
    d[s] = 0;
    deque<int> q({s});
    while(!q.empty()) {
        int a = q.front();
        q.pop_front();
        if(used[a]++) continue;
        for(auto e : g[a]) {
            if(d[e] > d[a] + e.cost) {
                d[e] = d[a] + e.cost;
                pre[e] = a;
                if(e.cost)
                    q.emplace_back(e);
                else
                    q.emplace_front(e);
            }
        }
    }
    if(d[t] > n) return pair(INT_MAX, vector<int>());
    int now = t;
    vector<int> res{t};
    while(now != s) {
        now = pre[now];
        res.eb(now);
    }
    reverse(begin(res), end(res));
    return {d[t], res};
}