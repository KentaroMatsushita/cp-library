struct oneDGraph {
    int n;
    vector<pair<int, int>> edges;
    vector<int> to;
    vector<int> start;
    oneDGraph(int _n) : n(_n), start(_n + 1) {}
    void reserve(int N) { edges.reserve(N); }
    void add(int u, int v) { edges.eb(u, v); }
    void build() {
        for(auto &[u, v] : edges) { start[u + 1]++; }
        for(int i = 0; i < n; i++) start[i + 1] += start[i];
        to.resize(edges.size());
        auto counter = start;
        for(auto &[u, v] : edges) to[counter[u]++] = v;
    }
};
