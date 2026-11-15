#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

struct WaveletMatrixBitVector {
    using u32 = uint32_t;
    using i64 = int64_t;
    using u64 = uint64_t;

    static constexpr u32 w = 64;
    std::vector<u64> block;
    std::vector<u32> count;
    u32 n = 0, zeros = 0;

    inline u32 get(u32 i) const { return u32(block[i / w] >> (i % w)) & 1u; }
    inline void set(u32 i) { block[i / w] |= 1ULL << (i % w); }

    WaveletMatrixBitVector() : block(1), count(1) {}
    WaveletMatrixBitVector(int _n) { init(_n); }
    void init(int _n) {
        if(_n < 0) throw std::invalid_argument("negative Wavelet Matrix size");
        n = zeros = _n;
        block.assign(n / w + 1, 0);
        count.assign(block.size(), 0);
    }

    void build() {
        for(u32 i = 1; i < block.size(); ++i) count[i] = count[i - 1] + __builtin_popcountll(block[i - 1]);
        zeros = rank0(n);
    }

    inline u32 rank0(u32 i) const { return i - rank1(i); }
    inline u32 rank1(u32 i) const {
        const u32 offset = i % w;
        const u64 mask = offset == 0 ? 0 : (u64(1) << offset) - 1;
        return count[i / w] + __builtin_popcountll(block[i / w] & mask);
    }
};

template <typename T, const int lg = 30> struct WaveletMatrix {
    static_assert(std::is_integral<T>::value && lg > 0 && lg <= std::numeric_limits<T>::digits, "invalid Wavelet Matrix value type or bit width");
    using u32 = uint32_t;
    using i64 = int64_t;
    using u64 = uint64_t;

    int n = 0;
    std::vector<T> a;
    std::array<WaveletMatrixBitVector, lg> bv;

    WaveletMatrix() = default;
    WaveletMatrix(u32 _n) : n(_n), a(n) {}
    WaveletMatrix(const std::vector<T> &_a) : n(_a.size()), a(_a) { build(); }

    void build() {
        validate_values();
        // lg = __lg(max<T>(*max_element(std::begin(a), std::end(a)), 1)) + 1;
        for(int i = 0; i < lg; i++) bv[i] = WaveletMatrixBitVector(n);
        std::vector<T> cur = a, nxt(n);
        for(int h = lg - 1; h >= 0; --h) {
            for(int i = 0; i < n; ++i)
                if((cur[i] >> h) & 1) bv[h].set(i);
            bv[h].build();
            std::array<decltype(std::begin(nxt)), 2> it{std::begin(nxt), std::begin(nxt) + bv[h].zeros};
            for(int i = 0; i < n; ++i) *it[bv[h].get(i)]++ = cur[i];
            std::swap(cur, nxt);
        }
        return;
    }

    void set(u32 i, const T &x) {
        validate_value(x);
        a[i] = x;
    }

    static void validate_value(T x) {
        if constexpr(std::is_signed<T>::value) {
            if(x < 0) throw std::invalid_argument("Wavelet Matrix values must be nonnegative");
        }
        if constexpr(lg < std::numeric_limits<T>::digits) {
            if(x >= (T(1) << lg)) throw std::out_of_range("Wavelet Matrix value exceeds bit width");
        }
    }
    void validate_values() {
        n = int(a.size());
        for(const T &x : a) validate_value(x);
    }

    inline std::pair<u32, u32> succ0(int l, int r, int h) const { return std::make_pair(bv[h].rank0(l), bv[h].rank0(r)); }

    inline std::pair<u32, u32> succ1(int l, int r, int h) const {
        u32 l0 = bv[h].rank0(l);
        u32 r0 = bv[h].rank0(r);
        u32 zeros = bv[h].zeros;
        return std::make_pair(l + zeros - l0, r + zeros - r0);
    }

    // return a[k]
    T access(u32 k) const {
        T ret = 0;
        for(int h = lg - 1; h >= 0; --h) {
            u32 f = bv[h].get(k);
            ret |= f ? T(1) << h : 0;
            k = f ? bv[h].rank1(k) + bv[h].zeros : bv[h].rank0(k);
        }
        return ret;
    }

    // k-th (0-indexed) smallest number in a[l, r)
    T kth_smallest(u32 l, u32 r, u32 k) const {
        T res = 0;
        for(int h = lg - 1; h >= 0; --h) {
            u32 l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if(k < r0 - l0)
                l = l0, r = r0;
            else {
                k -= r0 - l0;
                res |= (T)1 << h;
                l += bv[h].zeros - l0;
                r += bv[h].zeros - r0;
            }
        }
        return res;
    }

    // k-th (0-indexed) largest number in a[l, r)
    T kth_largest(int l, int r, int k) { return kth_smallest(l, r, r - l - k - 1); }

    // count i s.t. (l <= i < r) && (v[i] < upper)
    int range_freq(int l, int r, T upper) {
        if(upper <= 0) return 0;
        if constexpr(lg < std::numeric_limits<T>::digits) {
            if(upper >= (T(1) << lg)) return r - l;
        }
        int ret = 0;
        for(int h = lg - 1; h >= 0; --h) {
            bool f = (upper >> h) & 1;
            u32 l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if(f) {
                ret += r0 - l0;
                l += bv[h].zeros - l0;
                r += bv[h].zeros - r0;
            } else {
                l = l0;
                r = r0;
            }
        }
        return ret;
    }

    int range_freq(int l, int r, T lower, T upper) { return range_freq(l, r, upper) - range_freq(l, r, lower); }

    // max v[i] s.t. (l <= i < r) && (v[i] < upper)
    T prev_value(int l, int r, T upper) {
        int cnt = range_freq(l, r, upper);
        return cnt == 0 ? T(-1) : kth_smallest(l, r, cnt - 1);
    }

    // min v[i] s.t. (l <= i < r) && (lower <= v[i])
    T next_value(int l, int r, T lower) {
        int cnt = range_freq(l, r, lower);
        return cnt == r - l ? T(-1) : kth_smallest(l, r, cnt);
    }
};

/*
 * @brief Wavelet Matrix
 * @docs docs/data-structure-2d/wavelet-matrix.md
 */

template <typename T, typename sum_t = T, const int lg = 31> struct WaveletMatrix_sum : WaveletMatrix<T, lg> {
    using u32 = int32_t;
    std::array<std::vector<sum_t>, lg> sums;
    std::vector<sum_t> acc;
    using WM = WaveletMatrix<T, lg>;
    using WM::a;
    using WM::bv;
    using WM::n;
    WaveletMatrix_sum() : acc(1) {}
    explicit WaveletMatrix_sum(u32 size) : WM(size), acc(size + 1) {}
    WaveletMatrix_sum(const std::vector<T> &_a) {
        a = _a;
        build2();
    }
    void build() { build2(); }
    void build2() {
        this->validate_values();
        acc.assign(n + 1, sum_t());
        for(int i = 0; i < n; i++) acc[i + 1] = acc[i] + a[i];
        for(int i = 0; i < lg; i++) bv[i] = WaveletMatrixBitVector(n), sums[i].assign(n + 1, sum_t());
        std::vector<T> cur = a, nxt(n);
        for(int h = lg - 1; h >= 0; --h) {
            for(int i = 0; i < n; ++i)
                if((cur[i] >> h) & 1) bv[h].set(i);
            bv[h].build();
            std::array<decltype(std::begin(nxt)), 2> it{std::begin(nxt), std::begin(nxt) + bv[h].zeros};
            for(int i = 0; i < n; ++i) *it[bv[h].get(i)]++ = cur[i];
            std::swap(cur, nxt);
            for(int i = 0; i < n; i++) sums[h][i + 1] = sums[h][i] + cur[i];
        }
    }

    sum_t bottom_k_sum(int l, int r, int k) {
        sum_t res = 0;
        for(int h = lg - 1; h >= 0; --h) {
            u32 l0 = bv[h].rank0(l), r0 = bv[h].rank0(r);
            if(k < r0 - l0) {
                l = l0, r = r0;
            } else {
                res += sums[h][r0] - sums[h][l0];
                k -= r0 - l0;
                l += bv[h].zeros - l0;
                r += bv[h].zeros - r0;
            }
        }
        res += sums[0][l + k] - sums[0][l];
        return res;
    }

    sum_t top_k_sum(int l, int r, int k) { return acc[r] - acc[l] - bottom_k_sum(l, r, r - l - k); }
};
