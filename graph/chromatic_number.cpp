int chromatic_number(const vector<vector<int>> &g) {
    const int n = g.size();
    if(!n) return 0;
    vector<int> es(n);
    for(int i = 0; i < n; i++)
        for(auto e : g[i]) es[i] |= 1 << e;
    const int mod = 1e9 + 7;
    vector<int> I(1 << n), f(1 << n, 1);
    I[0] = 1;
    for(int b = 1; b < 1 << n; b++) {
        int t = __builtin_ctz(b);
        I[b] = I[b ^ 1 << t] + I[(b ^ (1 << t)) & ~es[t]];
    }

    for(int i = 1; i < n; i++) {
        ll sum = 0;
        for(int j = 0; j < 1 << n; j++) {
            f[j] = (1LL * f[j] * I[j]) % mod;
            sum += (n - popcount(j) & 1 ? mod - f[j] : f[j]);
        }
        if(sum % mod) return i;
    }
    return n;
}