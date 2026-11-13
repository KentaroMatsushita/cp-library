
// 0 以上 2 ^ (6D) 未満の整数を取れる
// D : 3 -> 2.6 * 10 ^ 5
// D : 4 -> 1.6 * 10 ^ 7
// D : 5 -> 1.0 * 10 ^ 9
// D : 6 -> int の範囲全部
template <const int D = 6> struct FastSet {
    int n;
    using u64 = uint64_t;
    vector<u64> a[D];
    FastSet(int n_ = (1 << std::min(31, D * 6))) : n(n_) {
        for(int i = 0; i < D; i++) {
            n_ = (n_ + 63) >> 6;
            a[i].assign(n_, 0);
        }
    }
    bool empty() const { return !a[D - 1][0]; }
    bool contains(int x) const { return (a[0][x >> 6] >> (x & 63)) & 1; }
    void insert(int x) {
        for(int i = 0; i < D; i++) {
            const int y = x & 63;
            x >>= 6;
            a[i][x] |= 1ULL << y;
        }
    }
    void erase(int x) {
        for(int i = 0; i < D; i++) {
            const int y = x & 63;
            x >>= 6;
            if((a[i][x] &= ~(1ULL << y))) break;
        }
    }
    int next(int x) const {
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
    int max() const { return prev(n); }
    int min() const { return next(0); }
};

constexpr pair<ll, int> calc_base_info(ll b) {
    ll B = 1, lg = 0;
    while(B <= LLONG_MAX / b) B *= b, lg++;
    return {B, lg};
}

template <int b, long long B = calc_base_info(b).first, int lg = calc_base_info(b).second> struct Base_b_Number {
    int l, r, n;
    vector<long long> a;
    FastSet<4> fs;
    array<long long, lg> ps;
    // b^(_l) ~ b^(_r)
    Base_b_Number(int _l, int _r) : l(floor(_l, lg)), r(ceil(_r, lg) + 2), n(r - l), a(r - l), fs(n) {
        ps[0] = 1;
        rep(i, 1, lg) ps[i] = ps[i - 1] * b;
    }
    Base_b_Number(int _n) : Base_b_Number(0, _n) {}
    // add b^i * x
    void add(int i, i128 x) {
        x *= ps[i - floor(i, lg) * lg];
        i = floor(i, lg);
        i -= l;
        while(a[i] >= B or x) {
            int at = bool(a[i]);
            if(x) {
                a[i] += x % B;
                x /= B;
            }
            if(a[i] >= B) {
                a[i] -= B;
                x++;
            }
            if(at != bool(a[i])) {
                if(at)
                    fs.erase(i);
                else
                    fs.insert(i);
            }
            i++;
        }
    }

    // sub b^i * x
    void sub(int i, i128 x) {
        x *= ps[i - floor(i, lg) * lg];
        i = floor(i, lg);
        i -= l;
        while(a[i] <= -2 or x) {
            int at = bool(a[i]);
            if(x) {
                a[i] -= x % B;
                x /= B;
            }
            if(a[i] <= -2) {
                a[i] += B;
                x++;
            }
            if(at != bool(a[i])) {
                if(at)
                    fs.erase(i);
                else
                    fs.insert(i);
            }
            i++;
        }
    }

    int sgn() const {
        int i = fs.max();
        if(i == -1) return 0;
        return a[i] > 0 ? 1 : -1;
    }

    void set_zero() {
        for(int i = fs.min(); i < n; i = fs.min()) {
            a[i] = 0;
            fs.erase(i);
        }
    }

    // constexpr ostream &operator<<(ostream &os) const { return os << a; }
};

// 基数動的
struct Base_b_Number {
    int b, lg;
    long long B;
    int l, r, n;
    vector<long long> a;
    FastSet<4> fs;
    vector<long long> ps;
    // b^(_l) ~ b^(_r)
    Base_b_Number(int b, int _l, int _r)
        : b(b), B(calc_base_info(b).first), lg(calc_base_info(b).second), l(floor(_l, lg)), r(ceil(_r, lg) + 2), n(r - l), a(r - l), fs(n), ps(lg) {
        ps[0] = 1;
        rep(i, 1, lg) ps[i] = ps[i - 1] * b;
    }
    Base_b_Number(int b, int _n) : Base_b_Number(b, 0, _n) {}
    // add b^i * x
    void add(int i, i128 x) {
        x *= ps[i - floor(i, lg) * lg];
        i = floor(i, lg);
        i -= l;
        while(a[i] >= B or x) {
            int at = bool(a[i]);
            if(x) {
                a[i] += x % B;
                x /= B;
            }
            if(a[i] >= B) {
                a[i] -= B;
                x++;
            }
            if(at != bool(a[i])) {
                if(at)
                    fs.erase(i);
                else
                    fs.insert(i);
            }
            i++;
        }
    }

    // sub b^i * x
    void sub(int i, i128 x) {
        x *= ps[i - floor(i, lg) * lg];
        i = floor(i, lg);
        i -= l;
        while(a[i] <= -2 or x) {
            int at = bool(a[i]);
            if(x) {
                a[i] -= x % B;
                x /= B;
            }
            if(a[i] <= -2) {
                a[i] += B;
                x++;
            }
            if(at != bool(a[i])) {
                if(at)
                    fs.erase(i);
                else
                    fs.insert(i);
            }
            i++;
        }
    }

    int sgn() const {
        int i = fs.max();
        if(i == -1) return 0;
        return a[i] > 0 ? 1 : -1;
    }

    void set_zero() {
        for(int i = fs.min(); i < n; i = fs.min()) {
            a[i] = 0;
            fs.erase(i);
        }
    }

    // constexpr ostream &operator<<(ostream &os) const { return os << a; }
};