/**
 * @brief Topological-Sort(トポロジカルソート)
 * @docs docs/topological-sort.md
 */
vector<int> topological_sort(const auto &g) {
    const int N = (int)g.size();
    vector<int> deg(N);
    for(int i = 0; i < N; i++) {
        for(auto &to : g[i]) ++deg[to];
    }
    stack<int> st;
    for(int i = 0; i < N; i++) {
        if(deg[i] == 0) st.emplace(i);
    }
    vector<int> ord;
    while(!st.empty()) {
        auto p = st.top();
        st.pop();
        ord.emplace_back(p);
        for(auto &to : g[p]) {
            if(--deg[to] == 0) st.emplace(to);
        }
    }
    return ord;
}
/**
 * @brief Offline-Dag-Reachability(DAGの到達可能性クエリ)
 * @docs docs/offline-dag-reachability.md
 */
vector<int> offline_dag_reachability(const auto &g, vector<pair<int, int>> &qs) {
    const int N = (int)g.size();
    const int Q = (int)qs.size();
    auto ord = topological_sort(g);
    vector<int> ans(Q);
    for(int l = 0; l < Q; l += 64) {
        int r = min(Q, l + 64);
        vector<int64_t> dp(N);
        for(int k = l; k < r; k++) { dp[qs[k].first] |= int64_t(1) << (k - l); }
        for(auto &idx : ord) {
            for(auto &to : g[idx]) dp[to] |= dp[idx];
        }
        for(int k = l; k < r; k++) { ans[k] = (dp[qs[k].second] >> (k - l)) & 1; }
    }
    return ans;
}