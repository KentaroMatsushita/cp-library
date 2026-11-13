namespace Modular998 {
template <typename mint> vector<mint> BerlekampMassey(const vector<mint> &s) {
    const int N = (int)s.size();
    vector<mint> b, c;
    b.reserve(N + 1);
    c.reserve(N + 1);
    b.push_back(mint(1));
    c.push_back(mint(1));
    mint y = mint(1);
    for(int ed = 1; ed <= N; ed++) {
        int l = int(c.size()), m = int(b.size());
        mint x = 0;
        for(int i = 0; i < l; i++) x += c[i] * s[ed - l + i];
        b.emplace_back(mint(0));
        m++;
        if(x == mint(0)) continue;
        mint freq = x / y;
        if(l < m) {
            auto tmp = c;
            c.insert(begin(c), m - l, mint(0));
            for(int i = 0; i < m; i++) c[m - 1 - i] -= freq * b[m - 1 - i];
            b = tmp;
            y = x;
        } else {
            for(int i = 0; i < m; i++) c[l - 1 - i] -= freq * b[m - 1 - i];
        }
    }
    reverse(begin(c), end(c));
    return c;
}

template <typename T> struct Binomial {
    vector<T> f, g, h;
    Binomial(int MAX = 0) : f(1, T(1)), g(1, T(1)), h(1, T(1)) {
        while(MAX >= (int)f.size()) extend();
    }

    void extend() {
        int n = f.size();
        int m = n * 2;
        f.resize(m);
        g.resize(m);
        h.resize(m);
        for(int i = n; i < m; i++) f[i] = f[i - 1] * T(i);
        g[m - 1] = f[m - 1].inverse();
        h[m - 1] = g[m - 1] * f[m - 2];
        for(int i = m - 2; i >= n; i--) {
            g[i] = g[i + 1] * T(i + 1);
            h[i] = g[i] * f[i - 1];
        }
    }

    T fac(int i) {
        if(i < 0) return T(0);
        while(i >= (int)f.size()) extend();
        return f[i];
    }

    T finv(int i) {
        if(i < 0) return T(0);
        while(i >= (int)g.size()) extend();
        return g[i];
    }

    T inv(int i) {
        if(i < 0) return -inv(-i);
        while(i >= (int)h.size()) extend();
        return h[i];
    }

    T C(int n, int r) {
        if(n < 0 || n < r || r < 0) return T(0);
        return fac(n) * finv(n - r) * finv(r);
    }

    inline T operator()(int n, int r) { return C(n, r); }

    template <typename I> T multinomial(const vector<I> &r) {
        static_assert(is_integral<I>::value == true);
        int n = 0;
        for(auto &x : r) {
            if(x < 0) return T(0);
            n += x;
        }
        T res = fac(n);
        for(auto &x : r) res *= finv(x);
        return res;
    }

    template <typename I> T operator()(const vector<I> &r) { return multinomial(r); }

    T C_naive(int n, int r) {
        if(n < 0 || n < r || r < 0) return T(0);
        T ret = T(1);
        r = min(r, n - r);
        for(int i = 1; i <= r; ++i) ret *= inv(i) * (n--);
        return ret;
    }

    T P(int n, int r) {
        if(n < 0 || n < r || r < 0) return T(0);
        return fac(n) * finv(n - r);
    }

    T H(int n, int r) {
        if(n < 0 || r < 0) return T(0);
        return r == 0 ? 1 : C(n + r - 1, r);
    }
};

template <uint32_t mod> struct LazyMontgomeryModInt {
    using mint = LazyMontgomeryModInt;
    using i32 = int32_t;
    using u32 = uint32_t;
    using u64 = uint64_t;

    static constexpr u32 get_r() {
        u32 ret = mod;
        for(i32 i = 0; i < 4; ++i) ret *= 2 - mod * ret;
        return ret;
    }

    static constexpr u32 r = get_r();
    static constexpr u32 n2 = -u64(mod) % mod;
    static_assert(r * mod == 1, "invalid, r * mod != 1");
    static_assert(mod < (1 << 30), "invalid, mod >= 2 ^ 30");
    static_assert((mod & 1) == 1, "invalid, mod % 2 == 0");

    u32 a;

    constexpr LazyMontgomeryModInt() : a(0) {}
    constexpr LazyMontgomeryModInt(const int64_t &b) : a(reduce(u64(b % mod + mod) * n2)) {};

    static constexpr u32 reduce(const u64 &b) { return (b + u64(u32(b) * u32(-r)) * mod) >> 32; }

    constexpr mint &operator+=(const mint &b) {
        if(i32(a += b.a - 2 * mod) < 0) a += 2 * mod;
        return *this;
    }

    constexpr mint &operator-=(const mint &b) {
        if(i32(a -= b.a) < 0) a += 2 * mod;
        return *this;
    }

    constexpr mint &operator*=(const mint &b) {
        a = reduce(u64(a) * b.a);
        return *this;
    }

    constexpr mint &operator/=(const mint &b) {
        *this *= b.inverse();
        return *this;
    }

    constexpr mint operator+(const mint &b) const { return mint(*this) += b; }
    constexpr mint operator-(const mint &b) const { return mint(*this) -= b; }
    constexpr mint operator*(const mint &b) const { return mint(*this) *= b; }
    constexpr mint operator/(const mint &b) const { return mint(*this) /= b; }
    constexpr bool operator==(const mint &b) const { return (a >= mod ? a - mod : a) == (b.a >= mod ? b.a - mod : b.a); }
    constexpr bool operator!=(const mint &b) const { return (a >= mod ? a - mod : a) != (b.a >= mod ? b.a - mod : b.a); }
    constexpr mint operator-() const { return mint() - mint(*this); }

    constexpr mint pow(u64 n) const {
        mint ret(1), mul(*this);
        while(n > 0) {
            if(n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }

    constexpr mint inverse() const { return pow(mod - 2); }

    friend ostream &operator<<(ostream &os, const mint &b) { return os << b.get(); }

    friend istream &operator>>(istream &is, mint &b) {
        int64_t t;
        is >> t;
        b = LazyMontgomeryModInt<mod>(t);
        return (is);
    }

    constexpr u32 get() const {
        u32 ret = reduce(a);
        return ret >= mod ? ret - mod : ret;
    }

    static constexpr u32 get_mod() { return mod; }
};
#include <immintrin.h>

__attribute__((target("sse4.2"))) inline __m128i my128_mullo_epu32(const __m128i &a, const __m128i &b) { return _mm_mullo_epi32(a, b); }

__attribute__((target("sse4.2"))) inline __m128i my128_mulhi_epu32(const __m128i &a, const __m128i &b) {
    __m128i a13 = _mm_shuffle_epi32(a, 0xF5);
    __m128i b13 = _mm_shuffle_epi32(b, 0xF5);
    __m128i prod02 = _mm_mul_epu32(a, b);
    __m128i prod13 = _mm_mul_epu32(a13, b13);
    __m128i prod = _mm_unpackhi_epi64(_mm_unpacklo_epi32(prod02, prod13), _mm_unpackhi_epi32(prod02, prod13));
    return prod;
}

__attribute__((target("sse4.2"))) inline __m128i montgomery_mul_128(const __m128i &a, const __m128i &b, const __m128i &r, const __m128i &m1) {
    return _mm_sub_epi32(_mm_add_epi32(my128_mulhi_epu32(a, b), m1), my128_mulhi_epu32(my128_mullo_epu32(my128_mullo_epu32(a, b), r), m1));
}

__attribute__((target("sse4.2"))) inline __m128i montgomery_add_128(const __m128i &a, const __m128i &b, const __m128i &m2, const __m128i &m0) {
    __m128i ret = _mm_sub_epi32(_mm_add_epi32(a, b), m2);
    return _mm_add_epi32(_mm_and_si128(_mm_cmpgt_epi32(m0, ret), m2), ret);
}

__attribute__((target("sse4.2"))) inline __m128i montgomery_sub_128(const __m128i &a, const __m128i &b, const __m128i &m2, const __m128i &m0) {
    __m128i ret = _mm_sub_epi32(a, b);
    return _mm_add_epi32(_mm_and_si128(_mm_cmpgt_epi32(m0, ret), m2), ret);
}

__attribute__((target("avx2"))) inline __m256i my256_mullo_epu32(const __m256i &a, const __m256i &b) { return _mm256_mullo_epi32(a, b); }

__attribute__((target("avx2"))) inline __m256i my256_mulhi_epu32(const __m256i &a, const __m256i &b) {
    __m256i a13 = _mm256_shuffle_epi32(a, 0xF5);
    __m256i b13 = _mm256_shuffle_epi32(b, 0xF5);
    __m256i prod02 = _mm256_mul_epu32(a, b);
    __m256i prod13 = _mm256_mul_epu32(a13, b13);
    __m256i prod = _mm256_unpackhi_epi64(_mm256_unpacklo_epi32(prod02, prod13), _mm256_unpackhi_epi32(prod02, prod13));
    return prod;
}

__attribute__((target("avx2"))) inline __m256i montgomery_mul_256(const __m256i &a, const __m256i &b, const __m256i &r, const __m256i &m1) {
    return _mm256_sub_epi32(_mm256_add_epi32(my256_mulhi_epu32(a, b), m1), my256_mulhi_epu32(my256_mullo_epu32(my256_mullo_epu32(a, b), r), m1));
}

__attribute__((target("avx2"))) inline __m256i montgomery_add_256(const __m256i &a, const __m256i &b, const __m256i &m2, const __m256i &m0) {
    __m256i ret = _mm256_sub_epi32(_mm256_add_epi32(a, b), m2);
    return _mm256_add_epi32(_mm256_and_si256(_mm256_cmpgt_epi32(m0, ret), m2), ret);
}

__attribute__((target("avx2"))) inline __m256i montgomery_sub_256(const __m256i &a, const __m256i &b, const __m256i &m2, const __m256i &m0) {
    __m256i ret = _mm256_sub_epi32(a, b);
    return _mm256_add_epi32(_mm256_and_si256(_mm256_cmpgt_epi32(m0, ret), m2), ret);
}
namespace ntt_inner {
using u64 = uint64_t;
constexpr uint32_t get_pr(uint32_t mod) {
    if(mod == 2) return 1;
    u64 ds[32] = {};
    int idx = 0;
    u64 m = mod - 1;
    for(u64 i = 2; i * i <= m; ++i) {
        if(m % i == 0) {
            ds[idx++] = i;
            while(m % i == 0) m /= i;
        }
    }
    if(m != 1) ds[idx++] = m;

    uint32_t pr = 2;
    while(1) {
        int flg = 1;
        for(int i = 0; i < idx; ++i) {
            u64 a = pr, b = (mod - 1) / ds[i], r = 1;
            while(b) {
                if(b & 1) r = r * a % mod;
                a = a * a % mod;
                b >>= 1;
            }
            if(r == 1) {
                flg = 0;
                break;
            }
        }
        if(flg == 1) break;
        ++pr;
    }
    return pr;
}

constexpr int SZ_FFT_BUF = 1 << 23;
uint32_t _buf1[SZ_FFT_BUF] __attribute__((aligned(64)));
uint32_t _buf2[SZ_FFT_BUF] __attribute__((aligned(64)));
} // namespace ntt_inner

template <typename mint> struct NTT {
    static constexpr uint32_t mod = mint::get_mod();
    static constexpr uint32_t pr = ntt_inner::get_pr(mint::get_mod());
    static constexpr int level = __builtin_ctzll(mod - 1);
    mint dw[level], dy[level];
    mint *buf1, *buf2;

    constexpr NTT() {
        setwy(level);
        union raw_cast {
            mint dat;
            uint32_t _;
        };
        buf1 = &(((raw_cast *)(ntt_inner::_buf1))->dat);
        buf2 = &(((raw_cast *)(ntt_inner::_buf2))->dat);
    }

    constexpr void setwy(int k) {
        mint w[level], y[level];
        w[k - 1] = mint(pr).pow((mod - 1) / (1 << k));
        y[k - 1] = w[k - 1].inverse();
        for(int i = k - 2; i > 0; --i) w[i] = w[i + 1] * w[i + 1], y[i] = y[i + 1] * y[i + 1];
        dw[0] = dy[0] = w[1] * w[1];
        dw[1] = w[1], dy[1] = y[1], dw[2] = w[2], dy[2] = y[2];
        for(int i = 3; i < k; ++i) {
            dw[i] = dw[i - 1] * y[i - 2] * w[i];
            dy[i] = dy[i - 1] * w[i - 2] * y[i];
        }
    }

    __attribute__((target("avx2"))) void ntt(mint *a, int n) {
        int k = n ? __builtin_ctz(n) : 0;
        if(k == 0) return;
        if(k == 1) {
            mint a1 = a[1];
            a[1] = a[0] - a[1];
            a[0] = a[0] + a1;
            return;
        }
        if(k & 1) {
            int v = 1 << (k - 1);
            if(v < 8) {
                for(int j = 0; j < v; ++j) {
                    mint ajv = a[j + v];
                    a[j + v] = a[j] - ajv;
                    a[j] += ajv;
                }
            } else {
                const __m256i m0 = _mm256_set1_epi32(0);
                const __m256i m2 = _mm256_set1_epi32(mod + mod);
                int j0 = 0;
                int j1 = v;
                for(; j0 < v; j0 += 8, j1 += 8) {
                    __m256i T0 = _mm256_loadu_si256((__m256i *)(a + j0));
                    __m256i T1 = _mm256_loadu_si256((__m256i *)(a + j1));
                    __m256i naj = montgomery_add_256(T0, T1, m2, m0);
                    __m256i najv = montgomery_sub_256(T0, T1, m2, m0);
                    _mm256_storeu_si256((__m256i *)(a + j0), naj);
                    _mm256_storeu_si256((__m256i *)(a + j1), najv);
                }
            }
        }
        int u = 1 << (2 + (k & 1));
        int v = 1 << (k - 2 - (k & 1));
        mint one = mint(1);
        mint imag = dw[1];
        while(v) {
            if(v == 1) {
                mint ww = one, xx = one, wx = one;
                for(int jh = 0; jh < u;) {
                    ww = xx * xx, wx = ww * xx;
                    mint t0 = a[jh + 0], t1 = a[jh + 1] * xx;
                    mint t2 = a[jh + 2] * ww, t3 = a[jh + 3] * wx;
                    mint t0p2 = t0 + t2, t1p3 = t1 + t3;
                    mint t0m2 = t0 - t2, t1m3 = (t1 - t3) * imag;
                    a[jh + 0] = t0p2 + t1p3, a[jh + 1] = t0p2 - t1p3;
                    a[jh + 2] = t0m2 + t1m3, a[jh + 3] = t0m2 - t1m3;
                    xx *= dw[__builtin_ctz((jh += 4))];
                }
            } else if(v == 4) {
                const __m128i m0 = _mm_set1_epi32(0);
                const __m128i m1 = _mm_set1_epi32(mod);
                const __m128i m2 = _mm_set1_epi32(mod + mod);
                const __m128i r = _mm_set1_epi32(mint::r);
                const __m128i Imag = _mm_set1_epi32(imag.a);
                mint ww = one, xx = one, wx = one;
                for(int jh = 0; jh < u;) {
                    if(jh == 0) {
                        int j0 = 0;
                        int j1 = v;
                        int j2 = j1 + v;
                        int j3 = j2 + v;
                        int je = v;
                        for(; j0 < je; j0 += 4, j1 += 4, j2 += 4, j3 += 4) {
                            const __m128i T0 = _mm_loadu_si128((__m128i *)(a + j0));
                            const __m128i T1 = _mm_loadu_si128((__m128i *)(a + j1));
                            const __m128i T2 = _mm_loadu_si128((__m128i *)(a + j2));
                            const __m128i T3 = _mm_loadu_si128((__m128i *)(a + j3));
                            const __m128i T0P2 = montgomery_add_128(T0, T2, m2, m0);
                            const __m128i T1P3 = montgomery_add_128(T1, T3, m2, m0);
                            const __m128i T0M2 = montgomery_sub_128(T0, T2, m2, m0);
                            const __m128i T1M3 = montgomery_mul_128(montgomery_sub_128(T1, T3, m2, m0), Imag, r, m1);
                            _mm_storeu_si128((__m128i *)(a + j0), montgomery_add_128(T0P2, T1P3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j1), montgomery_sub_128(T0P2, T1P3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j2), montgomery_add_128(T0M2, T1M3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j3), montgomery_sub_128(T0M2, T1M3, m2, m0));
                        }
                    } else {
                        ww = xx * xx, wx = ww * xx;
                        const __m128i WW = _mm_set1_epi32(ww.a);
                        const __m128i WX = _mm_set1_epi32(wx.a);
                        const __m128i XX = _mm_set1_epi32(xx.a);
                        int j0 = jh * v;
                        int j1 = j0 + v;
                        int j2 = j1 + v;
                        int j3 = j2 + v;
                        int je = j1;
                        for(; j0 < je; j0 += 4, j1 += 4, j2 += 4, j3 += 4) {
                            const __m128i T0 = _mm_loadu_si128((__m128i *)(a + j0));
                            const __m128i T1 = _mm_loadu_si128((__m128i *)(a + j1));
                            const __m128i T2 = _mm_loadu_si128((__m128i *)(a + j2));
                            const __m128i T3 = _mm_loadu_si128((__m128i *)(a + j3));
                            const __m128i MT1 = montgomery_mul_128(T1, XX, r, m1);
                            const __m128i MT2 = montgomery_mul_128(T2, WW, r, m1);
                            const __m128i MT3 = montgomery_mul_128(T3, WX, r, m1);
                            const __m128i T0P2 = montgomery_add_128(T0, MT2, m2, m0);
                            const __m128i T1P3 = montgomery_add_128(MT1, MT3, m2, m0);
                            const __m128i T0M2 = montgomery_sub_128(T0, MT2, m2, m0);
                            const __m128i T1M3 = montgomery_mul_128(montgomery_sub_128(MT1, MT3, m2, m0), Imag, r, m1);
                            _mm_storeu_si128((__m128i *)(a + j0), montgomery_add_128(T0P2, T1P3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j1), montgomery_sub_128(T0P2, T1P3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j2), montgomery_add_128(T0M2, T1M3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j3), montgomery_sub_128(T0M2, T1M3, m2, m0));
                        }
                    }
                    xx *= dw[__builtin_ctz((jh += 4))];
                }
            } else {
                const __m256i m0 = _mm256_set1_epi32(0);
                const __m256i m1 = _mm256_set1_epi32(mod);
                const __m256i m2 = _mm256_set1_epi32(mod + mod);
                const __m256i r = _mm256_set1_epi32(mint::r);
                const __m256i Imag = _mm256_set1_epi32(imag.a);
                mint ww = one, xx = one, wx = one;
                for(int jh = 0; jh < u;) {
                    if(jh == 0) {
                        int j0 = 0;
                        int j1 = v;
                        int j2 = j1 + v;
                        int j3 = j2 + v;
                        int je = v;
                        for(; j0 < je; j0 += 8, j1 += 8, j2 += 8, j3 += 8) {
                            const __m256i T0 = _mm256_loadu_si256((__m256i *)(a + j0));
                            const __m256i T1 = _mm256_loadu_si256((__m256i *)(a + j1));
                            const __m256i T2 = _mm256_loadu_si256((__m256i *)(a + j2));
                            const __m256i T3 = _mm256_loadu_si256((__m256i *)(a + j3));
                            const __m256i T0P2 = montgomery_add_256(T0, T2, m2, m0);
                            const __m256i T1P3 = montgomery_add_256(T1, T3, m2, m0);
                            const __m256i T0M2 = montgomery_sub_256(T0, T2, m2, m0);
                            const __m256i T1M3 = montgomery_mul_256(montgomery_sub_256(T1, T3, m2, m0), Imag, r, m1);
                            _mm256_storeu_si256((__m256i *)(a + j0), montgomery_add_256(T0P2, T1P3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j1), montgomery_sub_256(T0P2, T1P3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j2), montgomery_add_256(T0M2, T1M3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j3), montgomery_sub_256(T0M2, T1M3, m2, m0));
                        }
                    } else {
                        ww = xx * xx, wx = ww * xx;
                        const __m256i WW = _mm256_set1_epi32(ww.a);
                        const __m256i WX = _mm256_set1_epi32(wx.a);
                        const __m256i XX = _mm256_set1_epi32(xx.a);
                        int j0 = jh * v;
                        int j1 = j0 + v;
                        int j2 = j1 + v;
                        int j3 = j2 + v;
                        int je = j1;
                        for(; j0 < je; j0 += 8, j1 += 8, j2 += 8, j3 += 8) {
                            const __m256i T0 = _mm256_loadu_si256((__m256i *)(a + j0));
                            const __m256i T1 = _mm256_loadu_si256((__m256i *)(a + j1));
                            const __m256i T2 = _mm256_loadu_si256((__m256i *)(a + j2));
                            const __m256i T3 = _mm256_loadu_si256((__m256i *)(a + j3));
                            const __m256i MT1 = montgomery_mul_256(T1, XX, r, m1);
                            const __m256i MT2 = montgomery_mul_256(T2, WW, r, m1);
                            const __m256i MT3 = montgomery_mul_256(T3, WX, r, m1);
                            const __m256i T0P2 = montgomery_add_256(T0, MT2, m2, m0);
                            const __m256i T1P3 = montgomery_add_256(MT1, MT3, m2, m0);
                            const __m256i T0M2 = montgomery_sub_256(T0, MT2, m2, m0);
                            const __m256i T1M3 = montgomery_mul_256(montgomery_sub_256(MT1, MT3, m2, m0), Imag, r, m1);
                            _mm256_storeu_si256((__m256i *)(a + j0), montgomery_add_256(T0P2, T1P3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j1), montgomery_sub_256(T0P2, T1P3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j2), montgomery_add_256(T0M2, T1M3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j3), montgomery_sub_256(T0M2, T1M3, m2, m0));
                        }
                    }
                    xx *= dw[__builtin_ctz((jh += 4))];
                }
            }
            u <<= 2;
            v >>= 2;
        }
    }

    __attribute__((target("avx2"))) void intt(mint *a, int n, int normalize = true) {
        int k = n ? __builtin_ctz(n) : 0;
        if(k == 0) return;
        if(k == 1) {
            mint a1 = a[1];
            a[1] = a[0] - a[1];
            a[0] = a[0] + a1;
            if(normalize) {
                a[0] *= mint(2).inverse();
                a[1] *= mint(2).inverse();
            }
            return;
        }
        int u = 1 << (k - 2);
        int v = 1;
        mint one = mint(1);
        mint imag = dy[1];
        while(u) {
            if(v == 1) {
                mint ww = one, xx = one, yy = one;
                u <<= 2;
                for(int jh = 0; jh < u;) {
                    ww = xx * xx, yy = xx * imag;
                    mint t0 = a[jh + 0], t1 = a[jh + 1];
                    mint t2 = a[jh + 2], t3 = a[jh + 3];
                    mint t0p1 = t0 + t1, t2p3 = t2 + t3;
                    mint t0m1 = (t0 - t1) * xx, t2m3 = (t2 - t3) * yy;
                    a[jh + 0] = t0p1 + t2p3, a[jh + 2] = (t0p1 - t2p3) * ww;
                    a[jh + 1] = t0m1 + t2m3, a[jh + 3] = (t0m1 - t2m3) * ww;
                    xx *= dy[__builtin_ctz(jh += 4)];
                }
            } else if(v == 4) {
                const __m128i m0 = _mm_set1_epi32(0);
                const __m128i m1 = _mm_set1_epi32(mod);
                const __m128i m2 = _mm_set1_epi32(mod + mod);
                const __m128i r = _mm_set1_epi32(mint::r);
                const __m128i Imag = _mm_set1_epi32(imag.a);
                mint ww = one, xx = one, yy = one;
                u <<= 2;
                for(int jh = 0; jh < u;) {
                    if(jh == 0) {
                        int j0 = 0;
                        int j1 = v;
                        int j2 = v + v;
                        int j3 = j2 + v;
                        for(; j0 < v; j0 += 4, j1 += 4, j2 += 4, j3 += 4) {
                            const __m128i T0 = _mm_loadu_si128((__m128i *)(a + j0));
                            const __m128i T1 = _mm_loadu_si128((__m128i *)(a + j1));
                            const __m128i T2 = _mm_loadu_si128((__m128i *)(a + j2));
                            const __m128i T3 = _mm_loadu_si128((__m128i *)(a + j3));
                            const __m128i T0P1 = montgomery_add_128(T0, T1, m2, m0);
                            const __m128i T2P3 = montgomery_add_128(T2, T3, m2, m0);
                            const __m128i T0M1 = montgomery_sub_128(T0, T1, m2, m0);
                            const __m128i T2M3 = montgomery_mul_128(montgomery_sub_128(T2, T3, m2, m0), Imag, r, m1);
                            _mm_storeu_si128((__m128i *)(a + j0), montgomery_add_128(T0P1, T2P3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j2), montgomery_sub_128(T0P1, T2P3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j1), montgomery_add_128(T0M1, T2M3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j3), montgomery_sub_128(T0M1, T2M3, m2, m0));
                        }
                    } else {
                        ww = xx * xx, yy = xx * imag;
                        const __m128i WW = _mm_set1_epi32(ww.a);
                        const __m128i XX = _mm_set1_epi32(xx.a);
                        const __m128i YY = _mm_set1_epi32(yy.a);
                        int j0 = jh * v;
                        int j1 = j0 + v;
                        int j2 = j1 + v;
                        int j3 = j2 + v;
                        int je = j1;
                        for(; j0 < je; j0 += 4, j1 += 4, j2 += 4, j3 += 4) {
                            const __m128i T0 = _mm_loadu_si128((__m128i *)(a + j0));
                            const __m128i T1 = _mm_loadu_si128((__m128i *)(a + j1));
                            const __m128i T2 = _mm_loadu_si128((__m128i *)(a + j2));
                            const __m128i T3 = _mm_loadu_si128((__m128i *)(a + j3));
                            const __m128i T0P1 = montgomery_add_128(T0, T1, m2, m0);
                            const __m128i T2P3 = montgomery_add_128(T2, T3, m2, m0);
                            const __m128i T0M1 = montgomery_mul_128(montgomery_sub_128(T0, T1, m2, m0), XX, r, m1);
                            __m128i T2M3 = montgomery_mul_128(montgomery_sub_128(T2, T3, m2, m0), YY, r, m1);
                            _mm_storeu_si128((__m128i *)(a + j0), montgomery_add_128(T0P1, T2P3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j2), montgomery_mul_128(montgomery_sub_128(T0P1, T2P3, m2, m0), WW, r, m1));
                            _mm_storeu_si128((__m128i *)(a + j1), montgomery_add_128(T0M1, T2M3, m2, m0));
                            _mm_storeu_si128((__m128i *)(a + j3), montgomery_mul_128(montgomery_sub_128(T0M1, T2M3, m2, m0), WW, r, m1));
                        }
                    }
                    xx *= dy[__builtin_ctz(jh += 4)];
                }
            } else {
                const __m256i m0 = _mm256_set1_epi32(0);
                const __m256i m1 = _mm256_set1_epi32(mod);
                const __m256i m2 = _mm256_set1_epi32(mod + mod);
                const __m256i r = _mm256_set1_epi32(mint::r);
                const __m256i Imag = _mm256_set1_epi32(imag.a);
                mint ww = one, xx = one, yy = one;
                u <<= 2;
                for(int jh = 0; jh < u;) {
                    if(jh == 0) {
                        int j0 = 0;
                        int j1 = v;
                        int j2 = v + v;
                        int j3 = j2 + v;
                        for(; j0 < v; j0 += 8, j1 += 8, j2 += 8, j3 += 8) {
                            const __m256i T0 = _mm256_loadu_si256((__m256i *)(a + j0));
                            const __m256i T1 = _mm256_loadu_si256((__m256i *)(a + j1));
                            const __m256i T2 = _mm256_loadu_si256((__m256i *)(a + j2));
                            const __m256i T3 = _mm256_loadu_si256((__m256i *)(a + j3));
                            const __m256i T0P1 = montgomery_add_256(T0, T1, m2, m0);
                            const __m256i T2P3 = montgomery_add_256(T2, T3, m2, m0);
                            const __m256i T0M1 = montgomery_sub_256(T0, T1, m2, m0);
                            const __m256i T2M3 = montgomery_mul_256(montgomery_sub_256(T2, T3, m2, m0), Imag, r, m1);
                            _mm256_storeu_si256((__m256i *)(a + j0), montgomery_add_256(T0P1, T2P3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j2), montgomery_sub_256(T0P1, T2P3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j1), montgomery_add_256(T0M1, T2M3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j3), montgomery_sub_256(T0M1, T2M3, m2, m0));
                        }
                    } else {
                        ww = xx * xx, yy = xx * imag;
                        const __m256i WW = _mm256_set1_epi32(ww.a);
                        const __m256i XX = _mm256_set1_epi32(xx.a);
                        const __m256i YY = _mm256_set1_epi32(yy.a);
                        int j0 = jh * v;
                        int j1 = j0 + v;
                        int j2 = j1 + v;
                        int j3 = j2 + v;
                        int je = j1;
                        for(; j0 < je; j0 += 8, j1 += 8, j2 += 8, j3 += 8) {
                            const __m256i T0 = _mm256_loadu_si256((__m256i *)(a + j0));
                            const __m256i T1 = _mm256_loadu_si256((__m256i *)(a + j1));
                            const __m256i T2 = _mm256_loadu_si256((__m256i *)(a + j2));
                            const __m256i T3 = _mm256_loadu_si256((__m256i *)(a + j3));
                            const __m256i T0P1 = montgomery_add_256(T0, T1, m2, m0);
                            const __m256i T2P3 = montgomery_add_256(T2, T3, m2, m0);
                            const __m256i T0M1 = montgomery_mul_256(montgomery_sub_256(T0, T1, m2, m0), XX, r, m1);
                            const __m256i T2M3 = montgomery_mul_256(montgomery_sub_256(T2, T3, m2, m0), YY, r, m1);
                            _mm256_storeu_si256((__m256i *)(a + j0), montgomery_add_256(T0P1, T2P3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j2), montgomery_mul_256(montgomery_sub_256(T0P1, T2P3, m2, m0), WW, r, m1));
                            _mm256_storeu_si256((__m256i *)(a + j1), montgomery_add_256(T0M1, T2M3, m2, m0));
                            _mm256_storeu_si256((__m256i *)(a + j3), montgomery_mul_256(montgomery_sub_256(T0M1, T2M3, m2, m0), WW, r, m1));
                        }
                    }
                    xx *= dy[__builtin_ctz(jh += 4)];
                }
            }
            u >>= 4;
            v <<= 2;
        }
        if(k & 1) {
            v = 1 << (k - 1);
            if(v < 8) {
                for(int j = 0; j < v; ++j) {
                    mint ajv = a[j] - a[j + v];
                    a[j] += a[j + v];
                    a[j + v] = ajv;
                }
            } else {
                const __m256i m0 = _mm256_set1_epi32(0);
                const __m256i m2 = _mm256_set1_epi32(mod + mod);
                int j0 = 0;
                int j1 = v;
                for(; j0 < v; j0 += 8, j1 += 8) {
                    const __m256i T0 = _mm256_loadu_si256((__m256i *)(a + j0));
                    const __m256i T1 = _mm256_loadu_si256((__m256i *)(a + j1));
                    __m256i naj = montgomery_add_256(T0, T1, m2, m0);
                    __m256i najv = montgomery_sub_256(T0, T1, m2, m0);
                    _mm256_storeu_si256((__m256i *)(a + j0), naj);
                    _mm256_storeu_si256((__m256i *)(a + j1), najv);
                }
            }
        }
        if(normalize) {
            mint invn = mint(n).inverse();
            for(int i = 0; i < n; i++) a[i] *= invn;
        }
    }

    __attribute__((target("avx2"))) void inplace_multiply(int l1, int l2, int zero_padding = true) {
        int l = l1 + l2 - 1;
        int M = 4;
        while(M < l) M <<= 1;
        if(zero_padding) {
            for(int i = l1; i < M; i++) ntt_inner::_buf1[i] = 0;
            for(int i = l2; i < M; i++) ntt_inner::_buf2[i] = 0;
        }
        const __m256i m0 = _mm256_set1_epi32(0);
        const __m256i m1 = _mm256_set1_epi32(mod);
        const __m256i r = _mm256_set1_epi32(mint::r);
        const __m256i N2 = _mm256_set1_epi32(mint::n2);
        for(int i = 0; i < l1; i += 8) {
            __m256i a = _mm256_loadu_si256((__m256i *)(ntt_inner::_buf1 + i));
            __m256i b = montgomery_mul_256(a, N2, r, m1);
            _mm256_storeu_si256((__m256i *)(ntt_inner::_buf1 + i), b);
        }
        for(int i = 0; i < l2; i += 8) {
            __m256i a = _mm256_loadu_si256((__m256i *)(ntt_inner::_buf2 + i));
            __m256i b = montgomery_mul_256(a, N2, r, m1);
            _mm256_storeu_si256((__m256i *)(ntt_inner::_buf2 + i), b);
        }
        ntt(buf1, M);
        ntt(buf2, M);
        for(int i = 0; i < M; i += 8) {
            __m256i a = _mm256_loadu_si256((__m256i *)(ntt_inner::_buf1 + i));
            __m256i b = _mm256_loadu_si256((__m256i *)(ntt_inner::_buf2 + i));
            __m256i c = montgomery_mul_256(a, b, r, m1);
            _mm256_storeu_si256((__m256i *)(ntt_inner::_buf1 + i), c);
        }
        intt(buf1, M, false);
        const __m256i INVM = _mm256_set1_epi32((mint(M).inverse()).a);
        for(int i = 0; i < l; i += 8) {
            __m256i a = _mm256_loadu_si256((__m256i *)(ntt_inner::_buf1 + i));
            __m256i b = montgomery_mul_256(a, INVM, r, m1);
            __m256i c = my256_mulhi_epu32(my256_mullo_epu32(b, r), m1);
            __m256i d = _mm256_and_si256(_mm256_cmpgt_epi32(c, m0), m1);
            __m256i e = _mm256_sub_epi32(d, c);
            _mm256_storeu_si256((__m256i *)(ntt_inner::_buf1 + i), e);
        }
    }

    void ntt(vector<mint> &a) {
        int M = (int)a.size();
        for(int i = 0; i < M; i++) buf1[i].a = a[i].a;
        ntt(buf1, M);
        for(int i = 0; i < M; i++) a[i].a = buf1[i].a;
    }

    void intt(vector<mint> &a) {
        int M = (int)a.size();
        for(int i = 0; i < M; i++) buf1[i].a = a[i].a;
        intt(buf1, M, true);
        for(int i = 0; i < M; i++) a[i].a = buf1[i].a;
    }

    vector<mint> multiply(const vector<mint> &a, const vector<mint> &b) {
        if(a.size() == 0 || b.size() == 0) return vector<mint>{};
        int l = a.size() + b.size() - 1;
        if(min<int>(a.size(), b.size()) <= 40) {
            vector<mint> s(l);
            for(int i = 0; i < (int)a.size(); ++i)
                for(int j = 0; j < (int)b.size(); ++j) s[i + j] += a[i] * b[j];
            return s;
        }
        assert(l <= ntt_inner::SZ_FFT_BUF);
        int M = 4;
        while(M < l) M <<= 1;
        for(int i = 0; i < (int)a.size(); ++i) buf1[i].a = a[i].a;
        for(int i = (int)a.size(); i < M; ++i) buf1[i].a = 0;
        for(int i = 0; i < (int)b.size(); ++i) buf2[i].a = b[i].a;
        for(int i = (int)b.size(); i < M; ++i) buf2[i].a = 0;
        ntt(buf1, M);
        ntt(buf2, M);
        for(int i = 0; i < M; ++i) buf1[i].a = mint::reduce(uint64_t(buf1[i].a) * buf2[i].a);
        intt(buf1, M, false);
        vector<mint> s(l);
        mint invm = mint(M).inverse();
        for(int i = 0; i < l; ++i) s[i] = buf1[i] * invm;
        return s;
    }

    vector<mint> pow2(const vector<mint> &a) {
        if(a.size() == 0) return vector<mint>{};
        int l = a.size() + a.size() - 1;
        if(min<int>(a.size(), a.size()) <= 40) {
            vector<mint> s(l);
            for(int i = 0; i < (int)a.size(); ++i)
                for(int j = 0; j < (int)a.size(); ++j) s[i + j] += a[i] * a[j];
            return s;
        }
        assert(l <= ntt_inner::SZ_FFT_BUF);
        int M = 4;
        while(M < l) M <<= 1;
        for(int i = 0; i < (int)a.size(); ++i) buf1[i].a = a[i].a;
        for(int i = (int)a.size(); i < M; ++i) buf1[i].a = 0;
        // for(int i = 0; i < (int)b.size(); ++i) buf2[i].a = b[i].a;
        // for(int i = (int)b.size(); i < M; ++i) buf2[i].a = 0;
        ntt(buf1, M);
        // ntt(buf2, M);
        for(int i = 0; i < M; ++i) buf1[i].a = mint::reduce(uint64_t(buf1[i].a) * buf1[i].a);
        intt(buf1, M, false);
        vector<mint> s(l);
        mint invm = mint(M).inverse();
        for(int i = 0; i < l; ++i) s[i] = buf1[i] * invm;
        return s;
    }

    void ntt_doubling(vector<mint> &a) {
        int M = (int)a.size();
        for(int i = 0; i < M; i++) buf1[i].a = a[i].a;
        intt(buf1, M);
        mint r = 1, zeta = mint(pr).pow((mint::get_mod() - 1) / (M << 1));
        for(int i = 0; i < M; i++) buf1[i] *= r, r *= zeta;
        ntt(buf1, M);
        a.resize(2 * M);
        for(int i = 0; i < M; i++) a[M + i].a = buf1[i].a;
    }
};
template <typename mint> struct FormalPowerSeries : vector<mint> {
    using vector<mint>::vector;
    using FPS = FormalPowerSeries;

    FPS &operator+=(const FPS &r) {
        if(r.size() > this->size()) this->resize(r.size());
        for(int i = 0; i < (int)r.size(); i++) (*this)[i] += r[i];
        return *this;
    }

    FPS &operator+=(const mint &r) {
        if(this->empty()) this->resize(1);
        (*this)[0] += r;
        return *this;
    }

    FPS &operator-=(const FPS &r) {
        if(r.size() > this->size()) this->resize(r.size());
        for(int i = 0; i < (int)r.size(); i++) (*this)[i] -= r[i];
        return *this;
    }

    FPS &operator-=(const mint &r) {
        if(this->empty()) this->resize(1);
        (*this)[0] -= r;
        return *this;
    }

    FPS &operator*=(const mint &v) {
        for(int k = 0; k < (int)this->size(); k++) (*this)[k] *= v;
        return *this;
    }

    FPS &operator/=(const FPS &r) {
        if(this->size() < r.size()) {
            this->clear();
            return *this;
        }
        int n = this->size() - r.size() + 1;
        if((int)r.size() <= 64) {
            FPS f(*this), g(r);
            g.shrink();
            mint coeff = g.back().inverse();
            for(auto &x : g) x *= coeff;
            int deg = (int)f.size() - (int)g.size() + 1;
            int gs = g.size();
            FPS quo(deg);
            for(int i = deg - 1; i >= 0; i--) {
                quo[i] = f[i + gs - 1];
                for(int j = 0; j < gs; j++) f[i + j] -= quo[i] * g[j];
            }
            *this = quo * coeff;
            this->resize(n, mint(0));
            return *this;
        }
        return *this = ((*this).rev().pre(n) * r.rev().inv(n)).pre(n).rev();
    }

    FPS &operator%=(const FPS &r) {
        *this -= *this / r * r;
        shrink();
        return *this;
    }

    FPS operator+(const FPS &r) const { return FPS(*this) += r; }
    FPS operator+(const mint &v) const { return FPS(*this) += v; }
    FPS operator-(const FPS &r) const { return FPS(*this) -= r; }
    FPS operator-(const mint &v) const { return FPS(*this) -= v; }
    FPS operator*(const FPS &r) const { return FPS(*this) *= r; }
    FPS operator*(const mint &v) const { return FPS(*this) *= v; }
    FPS operator/(const FPS &r) const { return FPS(*this) /= r; }
    FPS operator%(const FPS &r) const { return FPS(*this) %= r; }
    FPS operator-() const {
        FPS ret(this->size());
        for(int i = 0; i < (int)this->size(); i++) ret[i] = -(*this)[i];
        return ret;
    }

    void shrink() {
        while(this->size() && this->back() == mint(0)) this->pop_back();
    }

    FPS rev() const {
        FPS ret(*this);
        reverse(begin(ret), end(ret));
        return ret;
    }

    FPS dot(FPS r) const {
        FPS ret(min(this->size(), r.size()));
        for(int i = 0; i < (int)ret.size(); i++) ret[i] = (*this)[i] * r[i];
        return ret;
    }

    FPS pre(int sz) const { return FPS(begin(*this), begin(*this) + min((int)this->size(), sz)); }

    FPS operator>>(int sz) const {
        if((int)this->size() <= sz) return {};
        FPS ret(*this);
        ret.erase(ret.begin(), ret.begin() + sz);
        return ret;
    }

    FPS operator<<(int sz) const {
        FPS ret(*this);
        ret.insert(ret.begin(), sz, mint(0));
        return ret;
    }

    FPS diff() const {
        const int n = (int)this->size();
        FPS ret(max(0, n - 1));
        mint one(1), coeff(1);
        for(int i = 1; i < n; i++) {
            ret[i - 1] = (*this)[i] * coeff;
            coeff += one;
        }
        return ret;
    }

    FPS integral() const {
        const int n = (int)this->size();
        FPS ret(n + 1);
        ret[0] = mint(0);
        if(n > 0) ret[1] = mint(1);
        auto mod = mint::get_mod();
        for(int i = 2; i <= n; i++) ret[i] = (-ret[mod % i]) * (mod / i);
        for(int i = 0; i < n; i++) ret[i + 1] *= (*this)[i];
        return ret;
    }

    mint eval(mint x) const {
        mint r = 0, w = 1;
        for(auto &v : *this) r += w * v, w *= x;
        return r;
    }

    FPS log(int deg = -1) const {
        assert((*this)[0] == mint(1));
        if(deg == -1) deg = (int)this->size();
        return (this->diff() * this->inv(deg)).pre(deg - 1).integral();
    }

    FPS pow(int64_t k, int deg = -1) const {
        const int n = (int)this->size();
        if(deg == -1) deg = n;
        if(k == 0) { FPS one(deg, mint(0)); if(deg) one[0] = 1; return one; }
        for(int i = 0; i < n; i++) {
            if((*this)[i] != mint(0)) {
                if(__int128(i) * k >= deg) return FPS(deg, mint(0));
                mint rev = mint(1) / (*this)[i];
                FPS ret = (((*this * rev) >> i).log(deg) * k).exp(deg) * ((*this)[i].pow(k));
                ret = (ret << (i * k)).pre(deg);
                if((int)ret.size() < deg) ret.resize(deg, mint(0));
                return ret;
            }
        }
        return FPS(deg, mint(0));
    }

    static void *ntt_ptr;
    static void set_fft();
    FPS &operator*=(const FPS &r);
    FPS pow2();
    void ntt();
    void intt();
    void ntt_doubling();
    static int ntt_pr();
    FPS inv(int deg = -1) const;
    FPS exp(int deg = -1) const;
};
template <typename mint> void *FormalPowerSeries<mint>::ntt_ptr = nullptr;

template <typename mint> void FormalPowerSeries<mint>::set_fft() {
    if(!ntt_ptr) ntt_ptr = new NTT<mint>;
}

template <typename mint> FormalPowerSeries<mint> &FormalPowerSeries<mint>::operator*=(const FormalPowerSeries<mint> &r) {
    if(this->empty() || r.empty()) {
        this->clear();
        return *this;
    }
    set_fft();
    auto ret = static_cast<NTT<mint> *>(ntt_ptr)->multiply(*this, r);
    return *this = FormalPowerSeries<mint>(ret.begin(), ret.end());
}

template <typename mint> FormalPowerSeries<mint> FormalPowerSeries<mint>::pow2() {
    if(this->empty()) {
        this->clear();
        return *this;
    }
    set_fft();
    auto ret = static_cast<NTT<mint> *>(ntt_ptr)->pow2(*this);
    return FormalPowerSeries<mint>(ret.begin(), ret.end());
}

template <typename mint> void FormalPowerSeries<mint>::ntt() {
    set_fft();
    static_cast<NTT<mint> *>(ntt_ptr)->ntt(*this);
}

template <typename mint> void FormalPowerSeries<mint>::intt() {
    set_fft();
    static_cast<NTT<mint> *>(ntt_ptr)->intt(*this);
}

template <typename mint> void FormalPowerSeries<mint>::ntt_doubling() {
    set_fft();
    static_cast<NTT<mint> *>(ntt_ptr)->ntt_doubling(*this);
}

template <typename mint> int FormalPowerSeries<mint>::ntt_pr() {
    set_fft();
    return static_cast<NTT<mint> *>(ntt_ptr)->pr;
}

template <typename mint> FormalPowerSeries<mint> FormalPowerSeries<mint>::inv(int deg) const {
    assert((*this)[0] != mint(0));
    if(deg == -1) deg = (int)this->size();
    FormalPowerSeries<mint> res(deg);
    res[0] = {mint(1) / (*this)[0]};
    for(int d = 1; d < deg; d <<= 1) {
        FormalPowerSeries<mint> f(2 * d), g(2 * d);
        for(int j = 0; j < min((int)this->size(), 2 * d); j++) f[j] = (*this)[j];
        for(int j = 0; j < d; j++) g[j] = res[j];
        f.ntt();
        g.ntt();
        for(int j = 0; j < 2 * d; j++) f[j] *= g[j];
        f.intt();
        for(int j = 0; j < d; j++) f[j] = 0;
        f.ntt();
        for(int j = 0; j < 2 * d; j++) f[j] *= g[j];
        f.intt();
        for(int j = d; j < min(2 * d, deg); j++) res[j] = -f[j];
    }
    return res.pre(deg);
}

template <typename mint> FormalPowerSeries<mint> FormalPowerSeries<mint>::exp(int deg) const {
    using fps = FormalPowerSeries<mint>;
    assert((*this).size() == 0 || (*this)[0] == mint(0));
    if(deg == -1) deg = this->size();

    fps inv;
    inv.reserve(deg + 1);
    inv.push_back(mint(0));
    inv.push_back(mint(1));

    auto inplace_integral = [&](fps &F) -> void {
        const int n = (int)F.size();
        auto mod = mint::get_mod();
        while((int)inv.size() <= n) {
            int i = inv.size();
            inv.push_back((-inv[mod % i]) * (mod / i));
        }
        F.insert(begin(F), mint(0));
        for(int i = 1; i <= n; i++) F[i] *= inv[i];
    };

    auto inplace_diff = [](fps &F) -> void {
        if(F.empty()) return;
        F.erase(begin(F));
        mint coeff = 1, one = 1;
        for(int i = 0; i < (int)F.size(); i++) {
            F[i] *= coeff;
            coeff += one;
        }
    };

    fps b{1, 1 < (int)this->size() ? (*this)[1] : 0}, c{1}, z1, z2{1, 1};
    for(int m = 2; m < deg; m *= 2) {
        auto y = b;
        y.resize(2 * m);
        y.ntt();
        z1 = z2;
        fps z(m);
        for(int i = 0; i < m; ++i) z[i] = y[i] * z1[i];
        z.intt();
        fill(begin(z), begin(z) + m / 2, mint(0));
        z.ntt();
        for(int i = 0; i < m; ++i) z[i] *= -z1[i];
        z.intt();
        c.insert(end(c), begin(z) + m / 2, end(z));
        z2 = c;
        z2.resize(2 * m);
        z2.ntt();
        fps x(begin(*this), begin(*this) + min<int>(this->size(), m));
        x.resize(m);
        inplace_diff(x);
        x.push_back(mint(0));
        x.ntt();
        for(int i = 0; i < m; ++i) x[i] *= y[i];
        x.intt();
        x -= b.diff();
        x.resize(2 * m);
        for(int i = 0; i < m - 1; ++i) x[m + i] = x[i], x[i] = mint(0);
        x.ntt();
        for(int i = 0; i < 2 * m; ++i) x[i] *= z2[i];
        x.intt();
        x.pop_back();
        inplace_integral(x);
        for(int i = m; i < min<int>(this->size(), 2 * m); ++i) x[i] += (*this)[i];
        fill(begin(x), begin(x) + m, mint(0));
        x.ntt();
        for(int i = 0; i < 2 * m; ++i) x[i] *= y[i];
        x.intt();
        b.insert(end(b), begin(x) + m, end(x));
    }
    return fps{begin(b), begin(b) + deg};
}

template <typename mint> vector<mint> FastMultiEval(const FormalPowerSeries<mint> &f, const vector<mint> &xs) {
    using fps = FormalPowerSeries<mint>;
    int s = xs.size();
    int N = 1 << (32 - __builtin_clz((int)xs.size() - 1));
    if(f.empty() || xs.empty()) return vector<mint>(s, mint(0));
    vector<FormalPowerSeries<mint>> buf(2 * N);
    for(int i = 0; i < N; i++) {
        mint n = mint{i < s ? -xs[i] : mint(0)};
        buf[i + N] = fps{n + 1, n - 1};
    }
    for(int i = N - 1; i > 0; i--) {
        fps &g(buf[(i << 1) | 0]), &h(buf[(i << 1) | 1]);
        int n = g.size();
        int m = n << 1;
        buf[i].reserve(m);
        buf[i].resize(n);
        for(int j = 0; j < n; j++) buf[i][j] = g[j] * h[j] - mint(1);
        if(i != 1) {
            buf[i].ntt_doubling();
            for(int j = 0; j < m; j++) buf[i][j] += j < n ? mint(1) : -mint(1);
        }
    }

    int fs = f.size();
    fps root = buf[1];
    root.intt();
    root.push_back(1);
    reverse(begin(root), end(root));
    root = root.inv(fs).rev() * f;
    root.erase(begin(root), begin(root) + fs - 1);
    root.resize(N, mint(0));

    vector<mint> ans(s);

    auto calc = [&](auto rec, int i, int l, int r, fps g) -> void {
        if(i >= N) {
            ans[i - N] = g[0];
            return;
        }
        int len = g.size(), m = (l + r) >> 1;
        g.ntt();
        fps tmp = buf[i * 2 + 1];
        for(int j = 0; j < len; j++) tmp[j] *= g[j];
        tmp.intt();
        rec(rec, i * 2 + 0, l, m, fps{begin(tmp) + (len >> 1), end(tmp)});
        if(m >= s) return;
        tmp = buf[i * 2 + 0];
        for(int j = 0; j < len; j++) tmp[j] *= g[j];
        tmp.intt();
        rec(rec, i * 2 + 1, m, r, fps{begin(tmp) + (len >> 1), end(tmp)});
    };
    calc(calc, 1, 0, N, root);
    return ans;
}

struct ArbitraryLazyMontgomeryModInt {
    using mint = ArbitraryLazyMontgomeryModInt;
    using i32 = int32_t;
    using u32 = uint32_t;
    using u64 = uint64_t;

    static u32 mod;
    static u32 r;
    static u32 n2;

    static u32 get_r() {
        u32 ret = mod;
        for(i32 i = 0; i < 4; ++i) ret *= 2 - mod * ret;
        return ret;
    }

    static void set_mod(u32 m) {
        assert(m < (1 << 30));
        assert((m & 1) == 1);
        mod = m;
        n2 = -u64(m) % m;
        r = get_r();
        assert(r * mod == 1);
    }

    u32 a;

    ArbitraryLazyMontgomeryModInt() : a(0) {}
    ArbitraryLazyMontgomeryModInt(const int64_t &b) : a(reduce(u64(b % mod + mod) * n2)) {};

    static u32 reduce(const u64 &b) { return (b + u64(u32(b) * u32(-r)) * mod) >> 32; }

    mint &operator+=(const mint &b) {
        if(i32(a += b.a - 2 * mod) < 0) a += 2 * mod;
        return *this;
    }

    mint &operator-=(const mint &b) {
        if(i32(a -= b.a) < 0) a += 2 * mod;
        return *this;
    }

    mint &operator*=(const mint &b) {
        a = reduce(u64(a) * b.a);
        return *this;
    }

    mint &operator/=(const mint &b) {
        *this *= b.inverse();
        return *this;
    }

    mint operator+(const mint &b) const { return mint(*this) += b; }
    mint operator-(const mint &b) const { return mint(*this) -= b; }
    mint operator*(const mint &b) const { return mint(*this) *= b; }
    mint operator/(const mint &b) const { return mint(*this) /= b; }
    bool operator==(const mint &b) const { return (a >= mod ? a - mod : a) == (b.a >= mod ? b.a - mod : b.a); }
    bool operator!=(const mint &b) const { return (a >= mod ? a - mod : a) != (b.a >= mod ? b.a - mod : b.a); }
    mint operator-() const { return mint() - mint(*this); }

    mint pow(u64 n) const {
        mint ret(1), mul(*this);
        while(n > 0) {
            if(n & 1) ret *= mul;
            mul *= mul;
            n >>= 1;
        }
        return ret;
    }

    friend ostream &operator<<(ostream &os, const mint &b) { return os << b.get(); }

    friend istream &operator>>(istream &is, mint &b) {
        int64_t t;
        is >> t;
        b = ArbitraryLazyMontgomeryModInt(t);
        return (is);
    }

    mint inverse() const { return pow(mod - 2); }

    u32 get() const {
        u32 ret = reduce(a);
        return ret >= mod ? ret - mod : ret;
    }

    static u32 get_mod() { return mod; }
};
typename ArbitraryLazyMontgomeryModInt::u32 ArbitraryLazyMontgomeryModInt::mod;
typename ArbitraryLazyMontgomeryModInt::u32 ArbitraryLazyMontgomeryModInt::r;
typename ArbitraryLazyMontgomeryModInt::u32 ArbitraryLazyMontgomeryModInt::n2;

int64_t mod_sqrt(const int64_t &a, const int64_t &p) {
    assert(0 <= a && a < p);
    if(a < 2) return a;
    using Mint = ArbitraryLazyMontgomeryModInt;
    Mint::set_mod(p);
    if(Mint(a).pow((p - 1) >> 1) != 1) return -1;
    Mint b = 1, one = 1;
    while(b.pow((p - 1) >> 1) == 1) b += one;
    int64_t m = p - 1, e = 0;
    while(m % 2 == 0) m >>= 1, e += 1;
    Mint x = Mint(a).pow((m - 1) >> 1);
    Mint y = Mint(a) * x * x;
    x *= a;
    Mint z = Mint(b).pow(m);
    while(y != 1) {
        int64_t j = 0;
        Mint t = y;
        while(t != one) {
            j += 1;
            t *= t;
        }
        z = z.pow(int64_t(1) << (e - j - 1));
        x *= z;
        z *= z;
        y *= z;
        e = j;
    }
    return x.get();
}

template <typename mint> FormalPowerSeries<mint> sqrt(const FormalPowerSeries<mint> &f, int deg = -1) {
    if(deg == -1) deg = (int)f.size();
    if((int)f.size() == 0) return FormalPowerSeries<mint>(deg, 0);
    if(f[0] == mint(0)) {
        for(int i = 1; i < (int)f.size(); i++) {
            if(f[i] != mint(0)) {
                if(i & 1) return {};
                if(deg - i / 2 <= 0) break;
                auto ret = sqrt(f >> i, deg - i / 2);
                if(ret.empty()) return {};
                ret = ret << (i / 2);
                if((int)ret.size() < deg) ret.resize(deg, mint(0));
                return ret;
            }
        }
        return FormalPowerSeries<mint>(deg, 0);
    }

    int64_t sqr = mod_sqrt(f[0].get(), mint::get_mod());
    if(sqr == -1) return {};
    assert(sqr * sqr % mint::get_mod() == f[0].get());
    FormalPowerSeries<mint> ret = {mint(sqr)};
    mint inv2 = mint(2).inverse();
    for(int i = 1; i < deg; i <<= 1) { ret = (ret + f.pre(i << 1) * ret.inv(i << 1)) * inv2; }
    return ret.pre(deg);
}

template <typename mint> mint LinearRecurrence(long long k, FormalPowerSeries<mint> Q, FormalPowerSeries<mint> P) {
    Q.shrink();
    mint ret = 0;
    if(P.size() >= Q.size()) {
        auto R = P / Q;
        P -= R * Q;
        P.shrink();
        if(k < (int)R.size()) ret += R[k];
    }
    if((int)P.size() == 0) return ret;

    FormalPowerSeries<mint>::set_fft();
    if(FormalPowerSeries<mint>::ntt_ptr == nullptr) {
        P.resize((int)Q.size() - 1);
        while(k) {
            auto Q2 = Q;
            for(int i = 1; i < (int)Q2.size(); i += 2) Q2[i] = -Q2[i];
            auto S = P * Q2;
            auto T = Q * Q2;
            if(k & 1) {
                for(int i = 1; i < (int)S.size(); i += 2) P[i >> 1] = S[i];
                for(int i = 0; i < (int)T.size(); i += 2) Q[i >> 1] = T[i];
            } else {
                for(int i = 0; i < (int)S.size(); i += 2) P[i >> 1] = S[i];
                for(int i = 0; i < (int)T.size(); i += 2) Q[i >> 1] = T[i];
            }
            k >>= 1;
        }
        return ret + P[0];
    } else {
        int N = 1;
        while(N < (int)Q.size()) N <<= 1;

        P.resize(2 * N);
        Q.resize(2 * N);
        P.ntt();
        Q.ntt();
        vector<mint> S(2 * N), T(2 * N);

        vector<int> btr(N);
        for(int i = 0, logn = __builtin_ctz(N); i < (1 << logn); i++) { btr[i] = (btr[i >> 1] >> 1) + ((i & 1) << (logn - 1)); }
        mint dw = mint(FormalPowerSeries<mint>::ntt_pr()).inverse().pow((mint::get_mod() - 1) / (2 * N));

        while(k) {
            mint inv2 = mint(2).inverse();

            // even degree of Q(x)Q(-x)
            T.resize(N);
            for(int i = 0; i < N; i++) T[i] = Q[(i << 1) | 0] * Q[(i << 1) | 1];

            S.resize(N);
            if(k & 1) {
                // odd degree of P(x)Q(-x)
                for(auto &i : btr) {
                    S[i] = (P[(i << 1) | 0] * Q[(i << 1) | 1] - P[(i << 1) | 1] * Q[(i << 1) | 0]) * inv2;
                    inv2 *= dw;
                }
            } else {
                // even degree of P(x)Q(-x)
                for(int i = 0; i < N; i++) { S[i] = (P[(i << 1) | 0] * Q[(i << 1) | 1] + P[(i << 1) | 1] * Q[(i << 1) | 0]) * inv2; }
            }
            swap(P, S);
            swap(Q, T);
            k >>= 1;
            if(k < N) break;
            P.ntt_doubling();
            Q.ntt_doubling();
        }
        P.intt();
        Q.intt();
        return ret + (P * (Q.inv()))[k];
    }
}

template <typename mint> mint kitamasa(long long N, FormalPowerSeries<mint> Q, FormalPowerSeries<mint> a) {
    assert(!Q.empty() && Q[0] != 0);
    if(N < (int)a.size()) return a[N];
    assert((int)a.size() >= int(Q.size()) - 1);
    auto P = a.pre((int)Q.size() - 1) * Q;
    P.resize(Q.size() - 1);
    dump(Q, P);
    return LinearRecurrence<mint>(N, Q, P);
}

// given  : y(x=0) , y(x=1) , ... , y(k)
// return : y(x)
template <typename mint> mint lagrange_interpolation(const vector<mint> &y, long long x, Binomial<mint> &C) {
    int N = (int)y.size() - 1;
    if(x <= N) return y[x];
    mint ret = 0;
    vector<mint> dp(N + 1, 1), pd(N + 1, 1);
    mint a = x, one = 1;
    for(int i = 0; i < N; i++) dp[i + 1] = dp[i] * a, a -= one;
    for(int i = N; i > 0; i--) pd[i - 1] = pd[i] * a, a += one;
    for(int i = 0; i <= N; i++) {
        mint tmp = y[i] * dp[i] * pd[i] * C.finv(i) * C.finv(N - i);
        ret += ((N - i) & 1) ? -tmp : tmp;
    }
    return ret;
}

template <typename mint> mint nth_term(long long n, const vector<mint> &s) {
    using fps = FormalPowerSeries<mint>;
    auto bm = BerlekampMassey<mint>(s);
    return kitamasa(n, fps{begin(bm), end(bm)}, fps{begin(s), end(s)});
}

using mint = LazyMontgomeryModInt<998244353>;
using fps = FormalPowerSeries<mint>;
using vmint = vector<mint>;
Binomial<mint> binomial;
mint inv(int i) { return binomial.inv(i); }
mint C(int r, int c) { return binomial.C(r, c); }
mint P(int r, int c) { return binomial.P(r, c); }
mint fact(int r) { return binomial.fac(r); }
mint ifact(int r) { return binomial.finv(r); }

struct Prd {
    deque<fps> deq;
    Prd() = default;
    void emplace(const fps &f) { deq.emplace_back(f); }
    fps calc() {
        if(deq.empty()) return {1};
        sort(all(deq), [&](const fps &f, const fps &g) { return si(f) < si(g); });
        while(deq.size() > 1) {
            deq.emplace_back(deq[0] * deq[1]);
            for(int i = 0; i < 2; ++i) deq.pop_front();
        }
        return deq.front();
    }
};
fps prd(vector<fps> &v) {
    Prd p;
    for(auto &e : v) p.emplace(e);
    return p.calc();
}

vmint power_table(mint x, int len) {
    vmint res(len + 1);
    res[0] = 1;
    rep(i, len) res[i + 1] = res[i] * x;
    return res;
}

// calc f(x + a)
fps TaylorShift(fps f, mint a) {
    int n = f.size();
    rep(i, n) f[i] *= fact(i);
    reverse(all(f));
    fps g(n, 1);
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
    static constexpr int MOD = 998244353;
    static const int M = sqrtl(MOD);
    rep(i, 1, M + 1) {
        if((x * i).a < M) return (i > 1 ? to_string((x * i).a) + " / " + to_string(i) : to_string((x * i).a));
        if(MOD - (x * i).a < M) return (i > 1 ? to_string(MOD - (x * i).a) + " / " + to_string(i) : to_string(MOD - (x * i).a));
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
            res += " " + c + " x^" + to_string(i);
            st = false;
        } else {
            if(!st)
                res += " +";
            else
                st = false;
            res += c + " x^" + to_string(i);
        }
    }
    return res;
}

} // namespace Modular998

using namespace Modular998;