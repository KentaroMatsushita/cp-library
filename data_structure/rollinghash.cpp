namespace RollingHashes {
constexpr ll MOD = (1LL << 61) - 1;

struct modint {

    ll a;

    constexpr modint(const ll x = 0) noexcept : a((x % MOD + MOD) % MOD) {}
    constexpr ll &value() noexcept { return a; }
    constexpr const ll &value() const noexcept { return a; }
    constexpr modint operator-() const noexcept { return modint() - *this; }
    constexpr modint operator+() const noexcept { return *this; }
    constexpr modint &operator++() noexcept {
        if(++a == MOD) a = 0;
        return *this;
    }
    constexpr modint &operator--() noexcept {
        if(!a) a = MOD;
        a--;
        return *this;
    }
    constexpr modint operator++(int) {
        modint res = *this;
        ++*this;
        return res;
    }
    constexpr modint operator--(int) {
        modint res = *this;
        --*this;
        return res;
    }
    constexpr modint &operator+=(const modint rhs) noexcept {
        a += rhs.a;
        if(a >= MOD) { a -= MOD; }
        return *this;
    }
    constexpr modint &operator-=(const modint rhs) noexcept {
        if(a < rhs.a) { a += MOD; }
        a -= rhs.a;
        return *this;
    }
    constexpr modint &operator*=(const modint rhs) noexcept {
        i128 t = (i128)(a) * (i128)(rhs.a);
        t = (t >> 61) + (t & MOD);
        if(t >= MOD) t -= MOD;
        a = t;
        return *this;
    }
    constexpr modint pow(long long n) const noexcept {
        if(n < 0) {
            n %= MOD - 1;
            n = (MOD - 1) + n;
        }
        modint x = *this, r = 1;
        while(n) {
            if(n & 1) r *= x;
            x *= x;
            n >>= 1;
        }
        return r;
    }
    constexpr modint &operator/=(const modint rhs) noexcept { return (*this *= rhs.pow(MOD - 2)); }
    constexpr modint inv() const noexcept { return pow(MOD - 2); }
    constexpr friend modint operator+(const modint &lhs, const modint &rhs) { return modint(lhs) += modint(rhs); }
    constexpr friend modint operator-(const modint &lhs, const modint &rhs) { return modint(lhs) -= modint(rhs); }
    constexpr friend modint operator*(const modint &lhs, const modint &rhs) { return modint(lhs) *= modint(rhs); }
    constexpr friend modint operator/(const modint &lhs, const modint &rhs) { return modint(lhs) /= modint(rhs); }
    constexpr friend bool operator==(const modint &lhs, const modint &rhs) { return lhs.a == rhs.a; }
    constexpr friend bool operator!=(const modint &lhs, const modint &rhs) { return lhs.a != rhs.a; }
    // constexpr friend modint operator^=(const modint &lhs, const modint &rhs) { return modint(lhs) ^= modint(rhs); }
};
using mint = modint;
using Hash = mint;
static mint base = 5454547;
vector<mint> p{1};
void calc_pow(int N) {
    while(si(p) < N) p.eb(p.back() * base);
}
struct RollingHash {
    vector<mint> h;
    template <typename T> RollingHash(const vector<T> &s) {
        int sz = (int)s.size();
        h.assign(sz + 1, 0);
        p[0] = 1;
        if(sz + 1 > (int)p.size()) {
            int psz = p.size();
            p.resize(sz + 1);
            for(int i = psz; i < sz + 1; i++) p[i] = p[i - 1] * base;
        }
        for(int i = 0; i < sz; i++) { h[i + 1] = h[i] * base + s[i]; }
    }
    RollingHash(const string &s) : RollingHash(vector<char>(begin(s), end(s))) {}

    ll get(int l, int r) const {
        mint res = h[r] + MOD - h[l] * p[r - l];
        return res.a;
    }

    ll connect(ull h1, ull h2, int h2len) const {
        mint res = h1 * p[h2len] + h2;
        return res.a;
    }

    int LCP(const RollingHash &b, int l1, int r1, int l2, int r2) {
        int len = min(r1 - l1, r2 - l2);
        int low = -1, high = len + 1;
        while(high - low > 1) {
            int mid = (low + high) / 2;
            if(get(l1, l1 + mid) == b.get(l2, l2 + mid))
                low = mid;
            else
                high = mid;
        }
        return (low);
    }
    RollingHash() = default;
};
template <class T> ll get_hash(const vector<T> &v) {
    mint res;
    rep(i, si(v)) res += p[i] * v[i];
    return res.a;
}
} // namespace RollingHashes
using namespace RollingHashes;
