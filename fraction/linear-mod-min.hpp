namespace fast_gcd {
using u64 = uint64_t;
using u32 = uint32_t;

__attribute__((target("bmi"))) constexpr u64 binary_gcd(u64 a, u64 b) {
    if(a == 0 || b == 0) return a + b;
    int n = __builtin_ctzll(a);
    int m = __builtin_ctzll(b);
    a >>= n;
    b >>= m;
    while(a != b) {
        int m = __builtin_ctzll(a - b);
        bool f = a > b;
        u64 c = f ? a : b;
        b = f ? b : a;
        a = (c - b) >> m;
    }
    return a << min(n, m);
}
} // namespace fast_gcd

using fast_gcd::binary_gcd;
// q / p
struct frac {
    template <typename T> static constexpr inline T gcd(const T a, const T b) { return (b == 0) ? a : gcd(b, a % b); }
    ll q, p;
    constexpr void simplify() {
        if(p < 0) {
            p *= -1;
            q *= -1;
        }
        ll g = binary_gcd(q < 0 ? -q : q, p);
        if(g) {
            p /= g;
            q /= g;
            if(p == 0 and q < 0) q = 1;
        }
    }
    constexpr frac(ll q = 0, ll p = 1) noexcept : q(q), p(p) { simplify(); }
    static frac raw(ll _q, ll _p) {
        frac f;
        f.q = _q, f.p = _p;
        return f;
    }

    constexpr bool operator<(const frac &r) const { return q * r.p < p * r.q; }
    constexpr bool operator>(const frac &r) const { return q * r.p > p * r.q; }
    constexpr bool operator<=(const frac &r) const { return q * r.p <= p * r.q; }
    constexpr bool operator>=(const frac &r) const { return q * r.p >= p * r.q; }
    template <class T> constexpr bool operator<(const T &r) const { return *this < frac(r); }
    template <class T> constexpr bool operator>(const T &r) const { return *this > frac(r); }
    constexpr bool operator==(const frac &r) const { return q * r.p == p * r.q; }
    constexpr bool operator!=(const frac &r) const { return !((*this) == r); }
    constexpr frac operator+() const noexcept { return *this; }
    constexpr frac operator-() const noexcept { return frac(-q, p); }
    constexpr frac operator+(const frac r) const noexcept { return frac(*this) += r; }
    constexpr frac operator-(const frac r) const noexcept { return frac(*this) -= r; }
    constexpr frac operator*(const frac r) const noexcept { return frac(*this) *= r; }
    constexpr frac operator/(const frac r) const noexcept { return frac(*this) /= r; }
    constexpr frac &operator+=(const frac &r) noexcept {
        ll g = binary_gcd(p, r.p);
        q = r.p / g * q + p / g * r.q;
        p = p / g * r.p;
        (*this).simplify();
        return *this;
    }
    constexpr frac &operator-=(const frac &r) noexcept {
        ll g = binary_gcd(p, r.p);
        q = r.p / g * q - p / g * r.q;
        p = p / g * r.p;
        (*this).simplify();
        return *this;
    }
    constexpr frac &operator*=(const frac &r) noexcept {
        q *= r.q;
        p *= r.p;
        (*this).simplify();
        return *this;
    }
    constexpr frac &operator/=(const frac &r) noexcept {
        q *= r.p;
        p *= r.q;
        (*this).simplify();
        return *this;
    }

    void print() {
        long double tmp = (long double)q / (long double)p;
        cout << tmp;
    }
};
istream &operator>>(istream &is, frac &p) {
    ll a;
    is >> a;
    p = frac(a, 1);
    return is;
}
ostream &operator<<(ostream &os, const frac &r) { return os << 1.0 * r.q / r.p; }

struct SBTnode;
struct raw_sbtnode {
    long long dep;
    frac l, r;
    raw_sbtnode() {
        dep = 0;
        l = frac::raw(0, 1), r = frac::raw(1, 0);
    }
    raw_sbtnode(const SBTnode &node);
    void ladd(long long y) {
        l.q += r.q * y;
        l.p += r.p * y;
    }
    void radd(long long y) {
        r.q += l.q * y;
        r.p += l.p * y;
    }
    void move(bool dir, long long d) {
        if(dir == false)
            radd(d);
        else
            ladd(d);
    }
    frac val() { return frac::raw(l.q + r.q, l.p + r.p); }
    operator frac() { return val(); }
};
struct SBTnode {
    static constexpr bool L = false, R = true;
    using Path = vector<pair<bool, long long>>;
    Path path;
    long long dep = 0;
    frac lval = frac::raw(0, 1), rval = frac::raw(1, 0);

    void ladd(long long y) {
        lval.q += rval.q * y;
        lval.p += rval.p * y;
    }
    void radd(long long y) {
        rval.q += lval.q * y;
        rval.p += lval.p * y;
    }

    frac val() { return frac::raw(lval.q + rval.q, lval.p + rval.p); }
    operator frac() { return val(); }

    // a / b
    SBTnode(long long a = 1, long long b = 1) {
        bool dir = a < b ? L : R;
        if(dir == L) swap(a, b);
        for(; b; dir = !dir) {
            long long q = a / b, r = a % b;
            move_down(dir, q - (!r));
            a = b, b = r;
        }
    }
    SBTnode(const frac &f) : SBTnode(f.q, f.p) {}
    SBTnode(const Path &p) : SBTnode() {
        for(const auto &[dir, d] : p) move_down(dir, d);
    }

    raw_sbtnode move_virtual(bool dir, long long d) {
        raw_sbtnode res(*this);
        res.move(dir, d);
        return res;
    }

    void move_down(const pair<bool, long long> &p) { return move_down(p.first, p.second); }
    void move_down(bool dir, long long d) {
        if(!d) return;
        if(empty(path) or path.back().first != dir) path.emplace_back(dir, 0);
        path.back().second += d, dep += d;
        if(dir == L)
            radd(d);
        else
            ladd(d);
    }

    raw_sbtnode move_up_virtual(long long d) {
        assert(d <= dep);
        raw_sbtnode res(*this);
        if(d == 0) return res;
        for(int i = (int)path.size() - 1; i >= 0; i--) {
            auto &[dir, t] = path[i];
            long long m = min(t, d);
            res.move(dir, -m);
            if(!(d -= m)) return res;
        }
        assert(0);
    }
    void move_up(long long d) {
        assert(d <= dep);
        while(d) {
            auto [dir, t] = path.back();
            path.pop_back();
            ll m = min(d, t);
            dep -= m, d -= m;
            if(dir == L)
                radd(-m);
            else
                ladd(-m);
            if(t > m) path.emplace_back(dir, t - m);
        }
    }

    static SBTnode lca(const SBTnode &x, const SBTnode &y) {
        SBTnode res;
        const auto &px = x.path;
        const auto &py = y.path;
        for(int i = 0; i < px.size(); i++) {
            if(i == py.size()) break;
            if(px[i] == py[i])
                res.move_down(px[i]);
            else {
                if(px[i].first == py[i].first) res.move_down(px[i].first, min(px[i].second, py[i].second));
                break;
            }
        }
        return res;
    }

    template <typename F> void bin_search(const F &f, long long max_p, long long max_q = -1) {
        auto check = [&](bool dir, long long d) {
            auto nxt = move_virtual(dir, d);
            auto nl = (dir == L ? nxt.r : nxt.l);
            return !(nl.p > max_p or (max_q != -1 and nl.q > max_q));
        };

        while(true) {
            bool flag = f(*this);
            bool dir = flag ? R : L;
            if(!check(dir, 1)) break;
            for(long long b = 1;; b <<= 1) {
                if(!check(dir, b) or f(move_virtual(dir, b)) != flag) {
                    b /= 2;
                    long long c = 0;
                    while(b) {
                        if(check(dir, b + c) and f(move_virtual(dir, b + c)) == flag) c += b;
                        b /= 2;
                    }
                    move_down(dir, c + (check(dir, c + 1)));
                    break;
                }
            }
        }
    }
};
raw_sbtnode::raw_sbtnode(const SBTnode &node) : dep(node.dep), l(node.lval), r(node.rval) {}

// f(0 / 1) : true のとき、分母が max_p, 分子が max_q 以下(設定しなかったら制約は分母のみ) で (true, false) の組になるような分数 が lval, rval
template <typename F> SBTnode sbt_bin_search(const F &f, long long max_p, long long max_q = -1) {
    SBTnode node;
    node.bin_search(f, max_p, max_q);
    return node;
}

// ax + b mod m の prefix minimum な x の情報
// x は log(m) 個の等差数列をなすので、変化点の座標と、その間の交差の vector を返す
// https://maspypy.com/library-checker-min-of-mod-of-linear
pair<vector<long long>, vector<long long>> mod_linear_prefix_minimum(long long a, long long b, long long m) {
    if(a >= m) a %= m;
    if(b >= m) b %= m;
    if((a = m - a) == m) a = 0;
    vector<long long> x, dx;
    x.emplace_back(0);
    long long X = 0, Y = b;
    long long g = gcd(a, m);
    SBTnode node = sbt_bin_search([&](frac f) { return f <= frac::raw(a, m); }, m / g);
    auto path = node.path;

    SBTnode tmp;
    if(long long c = a * tmp.lval.p - m * tmp.lval.q; c <= Y) {
        if(c == 0) return {x, dx};
        long long k = Y / c;
        x.emplace_back(X += tmp.lval.p * k);
        dx.emplace_back(tmp.lval.p);
        Y -= c * k;
    }
    for(int i = 0; i < path.size(); i++) {
        auto &[dir, d] = path[i];
        if(dir == SBTnode::L)
            tmp.move_down(dir, d);
        else {
            if(frac f = tmp.move_virtual(SBTnode::R, d).l; a * f.p - m * f.q <= Y) {
                long long t = ceil(a * tmp.lval.p - m * tmp.lval.q - Y, m * tmp.rval.q - a * tmp.rval.p);
                if(t < 0) assert(0);
                tmp.move_down(dir, t);

                long long c = a * tmp.lval.p - m * tmp.lval.q;
                if(c == 0) break;
                if(c < 0) assert(0);
                long long k = Y / c;
                x.emplace_back(X += tmp.lval.p * k);
                dx.emplace_back(tmp.lval.p);
                Y -= c * k;
                if(t < d) {
                    d -= t;
                    i--;
                }
            } else {
                tmp.move_down(dir, d);
            }
        }
    }
    return {x, dx};
}

// min ax + b mod m  (0 <= x < n)
long long min_mod_linear(long long a, long long b, long long m, long long n) {
    auto [x, dx] = mod_linear_prefix_minimum(a, b, m);
    ll tx = 0;
    if(x.back() < n)
        tx = x.back();
    else {
        int i = lower_bound(begin(x), end(x), n) - begin(x) - 1;
        tx = x[i] + (n - 1 - x[i]) / dx[i] * dx[i];
    }

    return (a * tx + b) % m;
}

pair<vector<long long>, vector<long long>> mod_linear_prefix_maximum(long long a, long long b, long long m) {
    long long ma = (m - a % m) % m;
    long long mb = (m - 1 - b % m);
    return mod_linear_prefix_minimum(ma, mb, m);
}
