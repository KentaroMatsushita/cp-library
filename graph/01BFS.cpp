vector<int> BFS(const Wgraph<int> &g, int s) {
    int n = g.size();
    vector<int> d(n, INT_MAX);
    d[s] = 0;
    deque<int> q({s});
    while(!q.empty()) {
        int a = q.front();
        q.pop_front();
        for(auto e : g[a]) {
            if(d[e] > d[a] + e.cost) {
                d[e] = d[a] + e.cost;
                if(e.cost)
                    q.emplace_back(e);
                else
                    q.emplace_front(e);
            }
        }
    }
    return d;
}