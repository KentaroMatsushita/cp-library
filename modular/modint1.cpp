namespace modular {
constexpr int MOD = 1000000007;
const int MAXN = 11000000;
template <int Modulus> class modint;
#define mint modint<MOD>
#define vmint vector<mint>
vector<mint> Inv;
mint inv(int x);
template <int Modulus> class modint {

  public:
    static constexpr int mod() { return Modulus; }
    int a;

    constexpr modint(const ll x = 0) noexcept : a(((x % Modulus) + Modulus) % Modulus) {}
    constexpr int &val() noexcept { return a; }
    constexpr const int &val() const noexcept { return a; }
    constexpr modint operator-() const noexcept { return modint() - *this; }
    constexpr modint operator+() const noexcept { return *this; }
    constexpr modint &operator++() noexcept {
        if(++a == Modulus) a = 0;
        return *this;
    }
    constexpr modint &operator--() noexcept {
        if(!a) a = Modulus;
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
        if(a >= Modulus) { a -= Modulus; }
        return *this;
    }
    constexpr modint &operator-=(const modint rhs) noexcept {
        if(a < rhs.a) { a += Modulus; }
        a -= rhs.a;
        return *this;
    }
    constexpr modint &operator*=(const modint rhs) noexcept {
        a = (long long)a * rhs.a % Modulus;
        return *this;
    }
    constexpr modint &operator/=(const modint rhs) noexcept {
        if constexpr(Modulus == MOD) a = (long long)a * modular::inv(rhs.a).a % Modulus;
        else a = (long long)a * rhs.inv().a % Modulus;
        return *this;
    }
    constexpr modint pow(long long n) const noexcept {
        if(n < 0) {
            n %= Modulus - 1;
            n = (Modulus - 1) + n;
        }
        modint x = *this, r = 1;
        while(n) {
            if(n & 1) r *= x;
            x *= x;
            n >>= 1;
        }
        return r;
    }
    constexpr modint inv() const noexcept { return pow(Modulus - 2); }
    constexpr friend modint operator+(const modint &lhs, const modint &rhs) { return modint(lhs) += modint(rhs); }
    constexpr friend modint operator-(const modint &lhs, const modint &rhs) { return modint(lhs) -= modint(rhs); }
    constexpr friend modint operator*(const modint &lhs, const modint &rhs) { return modint(lhs) *= modint(rhs); }
    constexpr friend modint operator/(const modint &lhs, const modint &rhs) { return modint(lhs) /= modint(rhs); }
    constexpr friend bool operator==(const modint &lhs, const modint &rhs) { return lhs.a == rhs.a; }
    constexpr friend bool operator!=(const modint &lhs, const modint &rhs) { return lhs.a != rhs.a; }
    // constexpr friend modint operator^=(const modint &lhs, const modint &rhs) { return modint(lhs) ^= modint(rhs); }
};
vmint Fact{1, 1}, Ifact{1, 1};
mint inv(int n) {
    if(n > MAXN) return (mint(n)).pow(MOD - 2);
    if(Inv.empty()) Inv.emplace_back(0), Inv.emplace_back(1);
    if(Inv.size() > n)
        return Inv[n];
    else {
        for(int i = Inv.size(); i <= n; ++i) {
            auto [y, x] = div(int(MOD), i);
            Inv.emplace_back(Inv[x] * (-y));
        }
        return Inv[n];
    }
}
mint fact(int n) {
    if(Fact.size() > n)
        return Fact[n];
    else
        for(int i = Fact.size(); i <= n; ++i) Fact.emplace_back(Fact[i - 1] * i);
    return Fact[n];
}
mint ifact(int n) {
    if(Ifact.size() > n)
        return Ifact[n];
    else
        for(int i = Ifact.size(); i <= n; ++i) Ifact.emplace_back(Ifact[i - 1] * inv(i));
    return Ifact[n];
}
mint modpow(ll a, ll n) { return mint(a).pow(n); }
mint inv(mint a) { return inv(a.a); }
mint ifact(mint a) { return ifact(a.a); }
mint fact(mint a) { return fact(a.a); }
mint modpow(mint a, ll n) { return modpow(a.a, n); }
mint C(int a, int b) {
    if(a < 0 || b < 0) return 0;
    if(a < b) return 0;
    if(a > MAXN) {
        mint res = 1;
        rep(i, b) res *= a - i, res /= i + 1;
        return res;
    }
    return fact(a) * ifact(b) * ifact(a - b);
}
mint P(int a, int b) {
    if(a < 0 || b < 0) return 0;
    if(a < b) return 0;
    if(a > MAXN) {
        mint res = 1;
        rep(i, b) res *= a - i;
        return res;
    }
    return fact(a) * ifact(a - b);
}
ostream &operator<<(ostream &os, mint a) {
    os << a.a;
    return os;
}
istream &operator>>(istream &is, mint &a) {
    ll x;
    is >> x;
    a = x;
    return is;
}
using Poly = vmint;
template <int mod> void FMT(vector<modint<mod>> &f, bool inv = false) {
    using V = vector<modint<mod>>;
    static V g(30), ig(30);
    if(g.front().a == 0) {
        modint<mod> root = 2;
        while((root.pow((mod - 1) / 2)).a == 1) root += 1;
        rep(i, 30) g[i] = -(root.pow((mod - 1) >> (i + 2))), ig[i] = g[i].inv();
    }
    int n = size(f);
    if(!inv) {
        for(int m = n; m >>= 1;) {
            modint<mod> w = 1;
            for(int s = 0, k = 0; s < n; s += 2 * m) {
                for(int i = s, j = s + m; i < s + m; ++i, ++j) {
                    auto x = f[i], y = f[j] * w;
                    // if(x.a >= mod) x.a -= mod;
                    f[i] = x + y, f[j] = x - y;
                    // f[i].a = x.a + y.a, f[j].a = x.a + (mod - y.a);
                    // if(f[i].a >= mod) f[i].a -= mod;
                    // if(f[j].a >= mod) f[j].a -= mod;
                }
                w *= g[__builtin_ctz(++k)];
            }
        }
    } else {
        for(int m = 1; m < n; m *= 2) {
            modint<mod> w = 1;
            for(int s = 0, k = 0; s < n; s += 2 * m) {
                for(int i = s, j = s + m; i < s + m; ++i, ++j) {
                    auto x = f[i], y = f[j];
                    f[i] = x + y, f[j] = x - y, f[j] *= w;
                    // f[i] = x + y, f[j].a = x.a + (mod - y.a), f[j] *= w;
                    // if(f[j].a >= mod) f[j].a -= mod;
                }
                w *= ig[__builtin_ctz(++k)];
            }
        }
    }
    modint<mod> c;
    if(inv)
        c = modint<mod>(n).inv();
    else
        c = 1;
    for(auto &&e : f) e *= c;
}
Poly operator-(Poly f) {
    for(auto &&e : f) e = -e;
    return f;
}
Poly &operator+=(Poly &l, const Poly &r) {
    l.resize(max(l.size(), r.size()));
    rep(i, r.size()) l[i] += r[i];
    return l;
}
Poly operator+(Poly l, const Poly &r) { return l += r; }
Poly &operator-=(Poly &l, const Poly &r) {
    l.resize(max(l.size(), r.size()));
    rep(i, r.size()) l[i] -= r[i];
    return l;
}
Poly operator-(Poly l, const Poly &r) { return l -= r; }
Poly &operator<<=(Poly &f, size_t n) { return f.insert(f.begin(), n, 0), f; }
Poly operator<<(Poly f, size_t n) { return f <<= n; }
Poly &operator>>=(Poly &f, size_t n) { return f.erase(f.begin(), f.begin() + min(f.size(), n)), f; }
Poly operator>>(Poly f, size_t n) { return f >>= n; }

namespace NTT {
constexpr int32_t mod0 = 167772161;
constexpr int32_t mod1 = 469762049;
constexpr int32_t mod2 = 754974721;
using M0 = modint<mod0>;
using M1 = modint<mod1>;
using M2 = modint<mod2>;

template <int mod> void mul(vector<modint<mod>> &l, vector<modint<mod>> &r) {
    int n = size(l), m = size(r), sz = 1 << __lg(2 * (n + m - 1) - 1);
    l.resize(sz), FMT<mod>(l);
    r.resize(sz), FMT<mod>(r);
    rep(i, sz) l[i] *= r[i];
    FMT<mod>(l, true);
    l.resize(n + m - 1);
}
Poly mulraw(const Poly &l, const Poly &r) {
    if(l.empty() or r.empty()) return Poly();
    if(min(si(l), si(r)) < 100) {
        Poly res(si(l) + si(r) - 1);
        rep(i, si(l)) rep(j, si(r)) res[i + j] += l[i] * r[j];
        return res;
    }
    int n = size(l), m = size(r);
    vector<M0> l0(n), r0(m);
    vector<M1> l1(n), r1(m);
    vector<M2> l2(n), r2(m);
    rep(i, n) l0[i] = l[i].a, l1[i] = l[i].a, l2[i] = l[i].a;
    rep(i, m) r0[i] = r[i].a, r1[i] = r[i].a, r2[i] = r[i].a;
    mul<mod0>(l0, r0), mul<mod1>(l1, r1), mul<mod2>(l2, r2);
    Poly res(n + m - 1);
    // garner
    static constexpr int r01 = M1(mod0).inv().a;
    static constexpr int r02 = M2(mod0).inv().a;
    static constexpr int r12 = M2(mod1).inv().a;
    static constexpr int r02r12 = (ll)(r02)*r12 % mod2;
    static constexpr ll w1 = mod0;
    static constexpr ll w2 = (ll)(mod0)*mod1;

    using i64 = int64_t;
    rep(i, n + m - 1) {
        int n1 = l1[i].a, n2 = l2[i].a, a = l0[i].a;
        int b = i64(n1 + mod1 - a) * r01 % mod1;
        int c = (i64(n2 + mod2 - a) * r02r12 + i64(mod2 - b) * r12) % mod2;
        res[i] = (__int128(a) + __int128(b) * w1 + __int128(c) * w2) % MOD;
    }
    return res;
}
} // namespace NTT
Poly operator*(const Poly &l, const Poly &r) { return NTT::mulraw(l, r); }
Poly &operator*=(Poly &l, const Poly &r) { return l = l * r; }
Poly integ(const Poly &f) {
    Poly res(f.size() + 1);
    for(int i = 1; i < (int)res.size(); ++i) res[i] = f[i - 1] / i;
    return res;
}
struct Prd {
    deque<Poly> deq;
    Prd() = default;
    void emplace(const Poly &f) { deq.emplace_back(f); }
    Poly calc() {
        if(deq.empty()) return {1};
        sort(all(deq), [&](const Poly &f, const Poly &g) { return si(f) < si(g); });
        while(deq.size() > 1) {
            deq.emplace_back(deq[0] * deq[1]);
            for(int i = 0; i < 2; ++i) deq.pop_front();
        }
        return deq.front();
    }
};
Poly prd(vector<Poly> &v) {
    Prd p;
    for(auto &e : v) p.emplace(e);
    return p.calc();
}
// Poly deriv(const Poly &f) {
//     if(f.size() == 0) return Poly();
//     Poly res(f.size() - 1);
//     rep(i, res.size()) res[i] = f[i + 1] * (i + 1);
//     return res;
// }
ostream &operator<<(ostream &os, const Poly &a) {
    for(auto e : a) cout << e.a << " ";
    return os;
}

vmint power_table(mint x, int len) {
    vmint res(len + 1);
    res[0] = 1;
    rep(i, len) res[i + 1] = res[i] * x;
    return res;
}

// calc f(x + a)
Poly TaylorShift(Poly f, mint a) {
    int n = f.size();
    rep(i, n) f[i] *= fact(i);
    reverse(all(f));
    Poly g(n, 1);
    rep(i, 1, n) g[i] = g[i - 1] * a * inv(i);
    f = (f * g);
    f.resize(n);
    reverse(begin(f), end(f));

    rep(i, n) f[i] *= ifact(i);
    return f;
}

// ボールの数、一個以上必要な数、入っていなくてもいい数(区別あり)
mint choose(int num, int a, int b = 0) {
    if(num == 0) return !a;
    return C(num + b - 1, a + b - 1);
}

// +1 n 個 -1 m 個で累積和 >= 0
mint Catalan(int n, int m) { return C(n + m, m) - C(n + m, m - 1); }

// +1 n 個 -1 m 個で累積和 > -k
mint Catalan2(int n, int m, int k) {
    if(m < k) return C(n + m, m);
    if(m < n + k) return C(n + m, m) - C(n + m, m - k);
    return 0;
}

// +1 n 個 -1 m 個で累積和 < +k
mint Catalan3(int n, int m, int k) { return Catalan2(m, n, k); }

string to_fraction(mint x) {
    static const int M = sqrtl(MOD);
    rep(i, 1, M + 1) {
        if((x * i).a < M) return (i > 1 ? to_string((x * i).a) + " / " + to_string(i) : to_string((x * i).a));
        if(MOD - (x * i).a < M) return "-" + (i > 1 ? to_string(MOD - (x * i).a) + " / " + to_string(i) : to_string(MOD - (x * i).a));
    }
    return "?";
}

string to_fraction(vmint x) {
    string res;
    bool st = true;
    rep(i, si(x)) {
        if(x[i].a == 0) continue;
        auto c = to_fraction(x[i]);
        if(c[0] == '-') {
            res += " " + c + (i ? " x^" + to_string(i) : "");
            st = false;
        } else {
            if(!st)
                res += " +";
            else
                st = false;
            res += c + (i ? " x^" + to_string(i) : "");
        }
    }
    return res;
}

mint eval(const vmint &f, mint x) {
    mint res, t = 1;
    fore(e, f) res += t * e, t *= x;
    return res;
}
} // namespace modular
using namespace modular;
