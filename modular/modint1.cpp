#pragma once
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <deque>
#include <istream>
#include <iterator>
#include <ostream>
#include <string>
#include <vector>

// Prime-modulus arithmetic and polynomial convolution. Self-contained C++17.
namespace modular {
constexpr int MOD = 1000000007;
const int MAXN = 11000000;
template <int Modulus> class modint;
using mint = modint<MOD>;
using vmint = std::vector<mint>;
inline mint inv(int x);
template <int Modulus> class modint {

  public:
    static constexpr int mod() { return Modulus; }
    int a;

    constexpr modint(const long long x = 0) noexcept : a(((x % Modulus) + Modulus) % Modulus) {}
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
inline vmint Inv, Fact{1, 1}, Ifact{1, 1};
inline mint inv(int n) {
    if(n > MAXN) return (mint(n)).pow(MOD - 2);
    if(Inv.empty()) Inv.emplace_back(0), Inv.emplace_back(1);
    if(Inv.size() > n)
        return Inv[n];
    else {
        for(int i = Inv.size(); i <= n; ++i) {
            auto [y, x] = std::div(int(MOD), i);
            Inv.emplace_back(Inv[x] * (-y));
        }
        return Inv[n];
    }
}
inline mint fact(int n) {
    if(Fact.size() > n)
        return Fact[n];
    else
        for(int i = Fact.size(); i <= n; ++i) Fact.emplace_back(Fact[i - 1] * i);
    return Fact[n];
}
inline mint ifact(int n) {
    if(Ifact.size() > n)
        return Ifact[n];
    else
        for(int i = Ifact.size(); i <= n; ++i) Ifact.emplace_back(Ifact[i - 1] * inv(i));
    return Ifact[n];
}
inline mint modpow(long long a, long long n) { return mint(a).pow(n); }
inline mint inv(mint a) { return inv(a.a); }
inline mint ifact(mint a) { return ifact(a.a); }
inline mint fact(mint a) { return fact(a.a); }
inline mint modpow(mint a, long long n) { return modpow(a.a, n); }
inline mint C(int a, int b) {
    if(a < 0 || b < 0) return 0;
    if(a < b) return 0;
    if(a > MAXN) {
        mint res = 1;
        for(int i = 0; i < b; ++i) res *= a - i, res /= i + 1;
        return res;
    }
    return fact(a) * ifact(b) * ifact(a - b);
}
inline mint P(int a, int b) {
    if(a < 0 || b < 0) return 0;
    if(a < b) return 0;
    if(a > MAXN) {
        mint res = 1;
        for(int i = 0; i < b; ++i) res *= a - i;
        return res;
    }
    return fact(a) * ifact(a - b);
}
inline std::ostream &operator<<(std::ostream &os, mint a) {
    os << a.a;
    return os;
}
inline std::istream &operator>>(std::istream &is, mint &a) {
    long long x;
    is >> x;
    a = x;
    return is;
}
using Poly = vmint;
template <int mod> void FMT(std::vector<modint<mod>> &f, bool inv = false) {
    using V = std::vector<modint<mod>>;
    static V g(30), ig(30);
    if(g.front().a == 0) {
        modint<mod> root = 2;
        while((root.pow((mod - 1) / 2)).a == 1) root += 1;
        for(int i = 0; i < 30; ++i) g[i] = -(root.pow((mod - 1) >> (i + 2))), ig[i] = g[i].inv();
    }
    int n = static_cast<int>(f.size());
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
inline Poly operator-(Poly f) {
    for(auto &&e : f) e = -e;
    return f;
}
inline Poly &operator+=(Poly &l, const Poly &r) {
    l.resize(std::max(l.size(), r.size()));
    for(std::size_t i = 0; i < r.size(); ++i) l[i] += r[i];
    return l;
}
inline Poly operator+(Poly l, const Poly &r) { return l += r; }
inline Poly &operator-=(Poly &l, const Poly &r) {
    l.resize(std::max(l.size(), r.size()));
    for(std::size_t i = 0; i < r.size(); ++i) l[i] -= r[i];
    return l;
}
inline Poly operator-(Poly l, const Poly &r) { return l -= r; }
inline Poly &operator<<=(Poly &f, std::size_t n) { return f.insert(f.begin(), n, 0), f; }
inline Poly operator<<(Poly f, std::size_t n) { return f <<= n; }
inline Poly &operator>>=(Poly &f, std::size_t n) { return f.erase(f.begin(), f.begin() + std::min(f.size(), n)), f; }
inline Poly operator>>(Poly f, std::size_t n) { return f >>= n; }

namespace NTT {
constexpr std::int32_t mod0 = 167772161;
constexpr std::int32_t mod1 = 469762049;
constexpr std::int32_t mod2 = 754974721;
using M0 = modint<mod0>;
using M1 = modint<mod1>;
using M2 = modint<mod2>;

template <int mod> void mul(std::vector<modint<mod>> &l, std::vector<modint<mod>> &r) {
    int n = static_cast<int>(l.size()), m = static_cast<int>(r.size()), sz = 1;
    if(n == 0 || m == 0) { l.clear(); return; }
    while(sz < n + m - 1) sz *= 2;
    l.resize(sz), FMT<mod>(l);
    r.resize(sz), FMT<mod>(r);
    for(int i = 0; i < sz; ++i) l[i] *= r[i];
    FMT<mod>(l, true);
    l.resize(n + m - 1);
}
inline Poly mulraw(const Poly &l, const Poly &r) {
    if(l.empty() or r.empty()) return Poly();
    if(std::min(l.size(), r.size()) < 100) {
        Poly res(l.size() + r.size() - 1);
        for(std::size_t i = 0; i < l.size(); ++i)
            for(std::size_t j = 0; j < r.size(); ++j) res[i + j] += l[i] * r[j];
        return res;
    }
    int n = static_cast<int>(l.size()), m = static_cast<int>(r.size());
    std::vector<M0> l0(n), r0(m);
    std::vector<M1> l1(n), r1(m);
    std::vector<M2> l2(n), r2(m);
    for(int i = 0; i < n; ++i) l0[i] = l[i].a, l1[i] = l[i].a, l2[i] = l[i].a;
    for(int i = 0; i < m; ++i) r0[i] = r[i].a, r1[i] = r[i].a, r2[i] = r[i].a;
    mul<mod0>(l0, r0), mul<mod1>(l1, r1), mul<mod2>(l2, r2);
    Poly res(n + m - 1);
    // garner
    static constexpr int r01 = M1(mod0).inv().a;
    static constexpr int r02 = M2(mod0).inv().a;
    static constexpr int r12 = M2(mod1).inv().a;
    static constexpr int r02r12 = (long long)(r02)*r12 % mod2;
    static constexpr long long w1 = mod0;
    static constexpr long long w2 = (long long)(mod0)*mod1;

    using i64 = std::int64_t;
    for(int i = 0; i < n + m - 1; ++i) {
        int n1 = l1[i].a, n2 = l2[i].a, a = l0[i].a;
        int b = i64(n1 + mod1 - a) * r01 % mod1;
        int c = (i64(n2 + mod2 - a) * r02r12 + i64(mod2 - b) * r12) % mod2;
        res[i] = (__int128(a) + __int128(b) * w1 + __int128(c) * w2) % MOD;
    }
    return res;
}
} // namespace NTT
inline Poly operator*(const Poly &l, const Poly &r) { return NTT::mulraw(l, r); }
inline Poly &operator*=(Poly &l, const Poly &r) { return l = l * r; }
inline Poly integ(const Poly &f) {
    Poly res(f.size() + 1);
    for(int i = 1; i < (int)res.size(); ++i) res[i] = f[i - 1] / i;
    return res;
}
struct Prd {
    std::deque<Poly> deq;
    Prd() = default;
    void emplace(const Poly &f) { deq.emplace_back(f); }
    Poly calc() {
        if(deq.empty()) return {1};
        std::sort(deq.begin(), deq.end(), [&](const Poly &f, const Poly &g) { return f.size() < g.size(); });
        while(deq.size() > 1) {
            deq.emplace_back(deq[0] * deq[1]);
            for(int i = 0; i < 2; ++i) deq.pop_front();
        }
        return deq.front();
    }
};
inline Poly prd(std::vector<Poly> &v) {
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
inline std::ostream &operator<<(std::ostream &os, const Poly &a) {
    for(auto e : a) os << e.a << " ";
    return os;
}

inline vmint power_table(mint x, int len) {
    vmint res(len + 1);
    res[0] = 1;
    for(int i = 0; i < len; ++i) res[i + 1] = res[i] * x;
    return res;
}

// calc f(x + a)
inline Poly TaylorShift(Poly f, mint a) {
    int n = f.size();
    for(int i = 0; i < n; ++i) f[i] *= fact(i);
    std::reverse(f.begin(), f.end());
    Poly g(n, 1);
    for(int i = 1; i < n; ++i) g[i] = g[i - 1] * a * inv(i);
    f = (f * g);
    f.resize(n);
    std::reverse(f.begin(), f.end());

    for(int i = 0; i < n; ++i) f[i] *= ifact(i);
    return f;
}

// ボールの数、一個以上必要な数、入っていなくてもいい数(区別あり)
inline mint choose(int num, int a, int b = 0) {
    if(num == 0) return !a;
    return C(num + b - 1, a + b - 1);
}

// +1 n 個 -1 m 個で累積和 >= 0
inline mint Catalan(int n, int m) { return C(n + m, m) - C(n + m, m - 1); }

// +1 n 個 -1 m 個で累積和 > -k
inline mint Catalan2(int n, int m, int k) {
    if(m < k) return C(n + m, m);
    if(m < n + k) return C(n + m, m) - C(n + m, m - k);
    return 0;
}

// +1 n 個 -1 m 個で累積和 < +k
inline mint Catalan3(int n, int m, int k) { return Catalan2(m, n, k); }

inline std::string to_fraction(mint x) {
    static const int M = std::sqrt(static_cast<long double>(MOD));
    for(int i = 1; i <= M; ++i) {
        if((x * i).a < M) return (i > 1 ? std::to_string((x * i).a) + " / " + std::to_string(i) : std::to_string((x * i).a));
        if(MOD - (x * i).a < M) return "-" + (i > 1 ? std::to_string(MOD - (x * i).a) + " / " + std::to_string(i) : std::to_string(MOD - (x * i).a));
    }
    return "?";
}

inline std::string to_fraction(vmint x) {
    std::string res;
    bool st = true;
    for(std::size_t i = 0; i < x.size(); ++i) {
        if(x[i].a == 0) continue;
        auto c = to_fraction(x[i]);
        if(c[0] == '-') {
            res += " " + c + (i ? " x^" + std::to_string(i) : "");
            st = false;
        } else {
            if(!st)
                res += " +";
            else
                st = false;
            res += c + (i ? " x^" + std::to_string(i) : "");
        }
    }
    return res;
}

inline mint eval(const vmint &f, mint x) {
    mint res, t = 1;
    for(const auto &e : f) res += t * e, t *= x;
    return res;
}
} // namespace modular
using namespace modular;
