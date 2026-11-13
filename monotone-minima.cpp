// monotone 行列の各行について、最小値を取る場所を返す
// f(i, j, k) : A[i][j] >= A[i][k] ? を返す関数 (j < k)
template <typename F> vector<int> monotone_minima(int h, int w, const F &f) {
    vector<int> dp(h, -1);
    auto rec = [&](auto &&rec, int u, int d, int l, int r) {
        if(u > d) return;
        int mid = u + d >> 1;
        auto &idx = dp[mid];
        idx = l;
        rep(i, l + 1, r + 1) if(f(mid, idx, i)) idx = i;
        rec(rec, u, mid - 1, l, idx);
        rec(rec, mid + 1, d, idx, r);
    };
    rec(rec, 0, h - 1, 0, w - 1);
    return dp;
}