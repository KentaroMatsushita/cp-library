
template <typename T> pair<vector<int>, T> dijkstra(Wgraph<T> &g, int s, int t, const T INF = numeric_limits<T>::max()) {
    vector<T> d(g.size(), INF);
    vector<int> pre(g.size(), -1);
    using P = pair<T, int>;
    priority_queue<P, vector<P>, greater<P>> pq;
    d[s] = 0;
    pq.emplace(0, s);
    while(!pq.empty()) {
        auto [cost, now] = pq.top();
        pq.pop();
        if(d[now] < cost) continue;
        for(auto &e : g[now]) {
            if(d[e] <= d[now] + e.cost) continue;
            d[e] = d[now] + e.cost;
            pre[e] = now;
            pq.emplace(d[e], e.to);
        }
    }
    if(d[t] == INF) return {{}, -1};
    vector<int> res{t};
    while(res.back() != s) { res.emplace_back(pre[res.back()]); }
    reverse(all(res));
    return {res, d[t]};
}
