vector<int> BFS(const Graph &g, int s) {
    int n = g.size();
    vector<int> d(n, INT_MAX);
    d[s] = 0;
    queue<int> q({s});
    while(!q.empty()) {
        int a = q.front();
        q.pop();
        for(auto e : g[a]) {
            if(d[e] > d[a] + 1) {
                d[e] = d[a] + 1;
                q.emplace(e);
            }
        }
    }
    return d;
}