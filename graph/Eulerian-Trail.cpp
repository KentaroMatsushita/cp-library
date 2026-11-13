template <typename T> vector<edge<T>> eulerian_path(Edges<T> es, int s = -1, bool directed = false) {
    if(es.empty()) return {};
    int V = 0;
    for(auto &e : es) V = max(V, max(e.to, e.from) + 1);
    if(s == -1) {
        if(directed) {
            vector<int> deg(V);
            for(auto &e : es) deg[e.from]++, deg[e.to]--;
            for(int i = 0; i < V; i++) {
                if(abs(deg[i]) > 1) assert("not eulerian graph");
                if(deg[i] == 1) s = i;
            }
            if(s == -1) s = es[0].from;
        } else {
            vector<int> deg(V);
            for(auto &e : es) deg[e.from] ^= 1, deg[e.to] ^= 1;
            int cnt = 0;
            for(int i = 0; i < V; i++) {
                if(deg[i]) cnt++, s = i;
            }
            if(cnt > 2) assert("not eulerian graph");
            if(s == -1) s = es[0].from;
        }
    }
    vector<vector<pair<edge<T>, int>>> g(V);
    for(auto &e : es) {
        int sz_to = (int)g[e.to].size();
        g[e.from].emplace_back(e, sz_to);
        if(!directed) {
            int sz_src = (int)g[e.from].size() - 1;
            swap(e.from, e.to);
            g[e.from].emplace_back(e, sz_src);
        }
    }
    vector<edge<T>> ord;
    stack<pair<int, edge<T>>> st;
    st.emplace(s, edge<T>(-1, -1));
    while(st.size()) {
        int idx = st.top().first;
        if(g[idx].empty()) {
            ord.emplace_back(st.top().second);
            st.pop();
        } else {
            auto e = g[idx].back();
            g[idx].pop_back();
            if(e.second == -1) continue;
            if(!directed) g[e.first.to][e.second].second = -1;
            st.emplace(e.first.to, e.first);
        }
    }
    ord.pop_back();
    reverse(begin(ord), end(ord));
    if(ord.size() != es.size()) {
        assert("not connected");
        return {};
    }
    return ord;
}