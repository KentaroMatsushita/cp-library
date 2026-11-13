template <class G> vector<int> road(int u, int v, const G &g) {
    int n = g.size();
    vector<int> a(n);
    queue<pii> q;
    swap(u, v);
    q.emplace(u, -1);
    while(!empty(q)) {
        auto [x, p] = q.front();
        q.pop();
        for(auto e : g[x])
            if(e != p) q.emplace(e, x), a[e] = x;
    }
    vector<int> res{v};
    while(res.back() != u) { res.emplace_back(a[res.back()]); }
    return res;
}