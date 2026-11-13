vector<vector<ll>> calc_binomail(int n, int m) {
    vector<vector<ll>> res(n + 1, vector<ll>(m + 1));
    res[0][0] = 1;
    for(int i = 1; i <= n; i++) {
        res[i][0] = 1;
        for(int j = 1; j <= min(i, m); j++) res[i][j] = res[i - 1][j - 1] + res[i - 1][j];
    }
    return res;
}