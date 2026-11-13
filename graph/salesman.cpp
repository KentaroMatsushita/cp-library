template <typename T> T TSP(const vector<vector<T>> &d, bool go_back = true, int start = -1) {
    int n = si(d);
    if(n == 0) return T{};
    const T INF = numeric_limits<T>::max() / 2;
    vector<vector<T>> dp(1 << n, vector<T>(n, INF));
    if(!go_back and start != -1)
        dp[1 << start][start] = 0;
    else if(!go_back)
        rep(i, n) dp[1 << i][i] = 0;
    else
        dp[1][0] = 0;
    rep(mask, 1 << n) rep(s, n) {
        if((mask & 1 << s) && dp[mask][s] != INF) rep(t, n) {
                if((~mask & 1 << t) && d[s][t] != INF) { chmin(dp[mask | 1 << t][t], dp[mask][s] + d[s][t]); }
            }
    }
    T res = INF;
    rep(i, n) {
        if(go_back) {
            if(dp[(1 << n) - 1][i] != INF) chmin(res, dp[(1 << n) - 1][i] + d[i][0]);
        } else
            chmin(res, dp[(1 << n) - 1][i]);
    }
    if(res == INF) res = -1;
    return res;
}
