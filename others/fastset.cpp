// 0 以上 2 ^ (6D) 未満の整数を取れる
// D : 3 -> 2.6 * 10 ^ 5
// D : 4 -> 1.6 * 10 ^ 7
// D : 5 -> 1.0 * 10 ^ 9
// D : 6 -> int の範囲全部
template <const int D = 6> struct FastSet {
    int n;
    using u64 = uint64_t;
    vector<u64> a[D];
    FastSet(int n_ = (1 << std::min(30, D * 6))) : n(n_) {
        static_assert(D >= 1 && D <= 6);
        if(n_ < 0 || uint64_t(n_) > (uint64_t(1) << (6 * D))) throw std::invalid_argument("FastSet capacity exceeds depth");
        for(int i = 0; i < D; i++) {
            n_ = std::max(1, int((uint64_t(n_) + 63) >> 6));
            a[i].assign(n_, 0);
        }
    }
    bool empty() const { return !a[D - 1][0]; }
    bool contains(int x) const { return x >= 0 && x < n && ((a[0][x >> 6] >> (x & 63)) & 1); }
    void insert(int x) {
        assert(0 <= x && x < n);
        for(int i = 0; i < D; i++) {
            const int y = x & 63;
            x >>= 6;
            a[i][x] |= 1ULL << y;
        }
    }
    void erase(int x) {
        assert(0 <= x && x < n);
        for(int i = 0; i < D; i++) {
            const int y = x & 63;
            x >>= 6;
            if((a[i][x] &= ~(1ULL << y))) break;
        }
    }
    int next(int x) const {
        x = std::max(x, 0);
        if(x >= n) return n;
        for(int i = 0; i < D; i++) {
            const int k = x >> 6, y = x & 63;
            if(k >= a[i].size()) return n;
            const u64 top = a[i][k] >> y;
            if(top) {
                x += __builtin_ctzll(top);
                for(int j = i - 1; j >= 0; --j) x = x << 6 | __builtin_ctzll(a[j][x]);
                return x;
            }
            x = k + 1;
        }
        return n;
    }
    int prev(int x) const {
        x = std::min(x, n - 1);
        for(int i = 0; i < D; ++i) {
            if(x < 0) return -1;
            const int k = x >> 6, y = x & 63;
            const u64 bot = a[i][k] << (63 - y);
            if(bot) {
                x -= __builtin_clzll(bot);
                for(int j = i - 1; j >= 0; --j) x = x << 6 | (63 - __builtin_clzll(a[j][x]));
                return x;
            }
            x = k - 1;
        }
        return -1;
    }
    int max() const { return prev(n - 1); }
    int min() const { return next(0); }
};
