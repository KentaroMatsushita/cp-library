struct BIT01 {
    int n;
    using u64 = uint64_t;
    static constexpr int B = 64;
    vector<int> a;
    vector<u64> v;
    BIT01(int sz) : a((sz + B - 1) / B + 1), v((sz + B - 1) / B + 1), n(sz) {}
    int sum(int k) {
        if(k < 0) return 0;
        if(k >= n) k = n - 1;
        int r = (k + 1) / B;
        int res = __builtin_popcountll(v[r] & ((1ULL << (k + 1 - r * B)) - 1));
        for(r; r > 0; r -= r & -r) res += a[r];
        return res;
    }
    // y = 0, 1
    void set(int x, int y) {
        int r = x / B;
        if(y ^ ((v[r] >> (x % B)) & 1)) {
            v[r] ^= 1ULL << (x % B);
            for(++r; r < a.size(); r += r & -r) a[r] += (y ? 1 : -1);
        }
    }
    int query(int l, int r) {
        if(l >= r) return 0;
        return sum(r - 1) - sum(l - 1);
    }
};