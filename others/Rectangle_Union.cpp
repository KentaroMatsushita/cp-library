// from ACL
namespace lazysegtree_impl {

namespace internal {

// @param n `0 <= n`
// @return minimum non-negative `x` s.t. `n <= 2**x`
int ceil_pow2(int n) {
    int x = 0;
    while((1U << x) < (unsigned int)(n)) x++;
    return x;
}

// @param n `1 <= n`
// @return minimum non-negative `x` s.t. `(n & (1 << x)) != 0`
constexpr int bsf_constexpr(unsigned int n) {
    int x = 0;
    while(!(n & (1 << x))) x++;
    return x;
}

// @param n `1 <= n`
// @return minimum non-negative `x` s.t. `(n & (1 << x)) != 0`
int bsf(unsigned int n) {
#ifdef _MSC_VER
    unsigned long index;
    _BitScanForward(&index, n);
    return index;
#else
    return __builtin_ctz(n);
#endif
}

} // namespace internal

} // namespace lazysegtree_impl

namespace lazysegtree_impl {
template <class S, S (*op)(S, S), S (*e)(), class F, S (*mapping)(F, S), F (*composition)(F, F), F (*id)()> struct lazy_segtree {
  public:
    lazy_segtree() : lazy_segtree(0) {}
    explicit lazy_segtree(int n) : lazy_segtree(std::vector<S>(n, e())) {}
    explicit lazy_segtree(const std::vector<S> &v) : _n(int(v.size())) {
        log = internal::ceil_pow2(_n);
        size = 1 << log;
        d = std::vector<S>(2 * size, e());
        lz = std::vector<F>(size, id());
        for(int i = 0; i < _n; i++) d[size + i] = v[i];
        for(int i = size - 1; i >= 1; i--) { update(i); }
    }

    void set(int p, S x) {
        assert(0 <= p && p < _n);
        p += size;
        for(int i = log; i >= 1; i--) push(p >> i);
        d[p] = x;
        for(int i = 1; i <= log; i++) update(p >> i);
    }

    S get(int p) {
        assert(0 <= p && p < _n);
        p += size;
        for(int i = log; i >= 1; i--) push(p >> i);
        return d[p];
    }

    S prod(int l, int r) {
        assert(0 <= l && l <= r && r <= _n);
        if(l == r) return e();

        l += size;
        r += size;

        for(int i = log; i >= 1; i--) {
            if(((l >> i) << i) != l) push(l >> i);
            if(((r >> i) << i) != r) push((r - 1) >> i);
        }

        S sml = e(), smr = e();
        while(l < r) {
            if(l & 1) sml = op(sml, d[l++]);
            if(r & 1) smr = op(d[--r], smr);
            l >>= 1;
            r >>= 1;
        }

        return op(sml, smr);
    }

    S all_prod() { return d[1]; }

    void apply(int p, F f) {
        assert(0 <= p && p < _n);
        p += size;
        for(int i = log; i >= 1; i--) push(p >> i);
        d[p] = mapping(f, d[p]);
        for(int i = 1; i <= log; i++) update(p >> i);
    }
    void apply(int l, int r, F f) {
        assert(0 <= l && l <= r && r <= _n);
        if(l == r) return;

        l += size;
        r += size;

        for(int i = log; i >= 1; i--) {
            if(((l >> i) << i) != l) push(l >> i);
            if(((r >> i) << i) != r) push((r - 1) >> i);
        }

        {
            int l2 = l, r2 = r;
            while(l < r) {
                if(l & 1) all_apply(l++, f);
                if(r & 1) all_apply(--r, f);
                l >>= 1;
                r >>= 1;
            }
            l = l2;
            r = r2;
        }

        for(int i = 1; i <= log; i++) {
            if(((l >> i) << i) != l) update(l >> i);
            if(((r >> i) << i) != r) update((r - 1) >> i);
        }
    }

    template <bool (*g)(S)> int max_right(int l) {
        return max_right(l, [](S x) { return g(x); });
    }
    template <class G> int max_right(int l, G g) {
        assert(0 <= l && l <= _n);
        assert(g(e()));
        if(l == _n) return _n;
        l += size;
        for(int i = log; i >= 1; i--) push(l >> i);
        S sm = e();
        do {
            while(l % 2 == 0) l >>= 1;
            if(!g(op(sm, d[l]))) {
                while(l < size) {
                    push(l);
                    l = (2 * l);
                    if(g(op(sm, d[l]))) {
                        sm = op(sm, d[l]);
                        l++;
                    }
                }
                return l - size;
            }
            sm = op(sm, d[l]);
            l++;
        } while((l & -l) != l);
        return _n;
    }

    template <bool (*g)(S)> int min_left(int r) {
        return min_left(r, [](S x) { return g(x); });
    }
    template <class G> int min_left(int r, G g) {
        assert(0 <= r && r <= _n);
        assert(g(e()));
        if(r == 0) return 0;
        r += size;
        for(int i = log; i >= 1; i--) push((r - 1) >> i);
        S sm = e();
        do {
            r--;
            while(r > 1 && (r % 2)) r >>= 1;
            if(!g(op(d[r], sm))) {
                while(r < size) {
                    push(r);
                    r = (2 * r + 1);
                    if(g(op(d[r], sm))) {
                        sm = op(d[r], sm);
                        r--;
                    }
                }
                return r + 1 - size;
            }
            sm = op(d[r], sm);
        } while((r & -r) != r);
        return 0;
    }

    S operator[](int k) { return get(k); }

    friend std::ostream &operator<<(std::ostream &os, lazy_segtree &ls) {
        os << "{";
        for(int i = 0; i < ls._n; i++) os << ls.get(i) << (i == ls._n - 1 ? "" : ", ");
        return os << "}";
    }

  private:
    int _n, size, log;
    std::vector<S> d;
    std::vector<F> lz;

    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
    void all_apply(int k, F f) {
        d[k] = mapping(f, d[k]);
        if(k < size) lz[k] = composition(f, lz[k]);
    }
    void push(int k) {
        all_apply(2 * k, lz[k]);
        all_apply(2 * k + 1, lz[k]);
        lz[k] = id();
    }
};

template <typename T> static constexpr T inf = numeric_limits<T>::max() / 2;
template <typename T> T min_func(T x, T y) { return min(x, y); }
template <typename T> T max_func(T x, T y) { return max(x, y); }
template <typename T> T plus_func(T x, T y) { return x + y; }
template <typename T, int ID> T set_func(T x, T y) { return (y == ID ? x : y); }
template <typename T> T min_e() { return inf<T>; }
template <typename T> T max_e() { return -inf<T>; }
template <typename T> T plus_e() { return 0; }

// range add get range min
template <typename T> using RARm = lazy_segtree<T, min_func<T>, min_e<T>, T, plus_func<T>, plus_func<T>, plus_e<T>>;
// range add get range max
template <typename T> using RARM = lazy_segtree<T, max_func<T>, max_e<T>, T, plus_func<T>, plus_func<T>, plus_e<T>>;
// range add get range sum
// first が size, second が value
template <typename T> pair<int, T> RARS_op(pair<int, T> x, pair<int, T> y) { return T{x.first + y.first, x.second + y.second}; }
template <typename T> pair<int, T> RARS_e() { return {0, 0}; }
template <typename T> pair<int, T> RARS_map(T f, pair<int, T> x) { return T{x.second + f * x.first, x.first}; }
template <typename T> using RARS = lazy_segtree<pair<int, T>, RARS_op<T>, RARS_e<T>, T, RARS_map<T>, plus_func<T>, plus_e<T>>;

// 区間の挿入、削除、union のサイズ取得
pair<int, int> segunion_merge(pair<int, int> x, pair<int, int> y) { return (x.se == y.se ? pair<int, int>{x.fi + y.fi, x.se} : x.se < y.se ? x : y); }
pair<int, int> segunion_e() { return {}; }
pair<int, int> segunion_mapp(int x, pair<int, int> p) { return {p.fi, p.se + x}; }
using SegmentUnion = lazy_segtree<pair<int, int>, segunion_merge, segunion_e, int, segunion_mapp, plus_func<int>, plus_e<int>>;
} // namespace lazysegtree_impl

using lazysegtree_impl::lazy_segtree;
using lazysegtree_impl::RARM;
using lazysegtree_impl::RARm;
using lazysegtree_impl::RARS;
using lazysegtree_impl::SegmentUnion;

namespace Rectangle_Union_impl {
template <typename T> struct S { T sz, mi, mi_cnt; };
template <typename T> S<T> f(S<T> x, S<T> y) { return S<T>{x.sz + y.sz, min(x.mi, y.mi), (x.mi <= y.mi ? x.mi_cnt : 0) + (y.mi <= x.mi ? y.mi_cnt : 0)}; }
template <typename T> S<T> e() { return S<T>{0, inf<T>, 0}; }
template <typename T> S<T> g(int x, S<T> s) {
    s.mi += x;
    return s;
}
int h(int x, int y) { return x + y; }
int id() { return 0; }
template <typename T> struct Rectangle_Union {
    vector<T> ys;
    vector<tuple<T, T, T, char>> rect;
    Rectangle_Union(int n = 0) { ys.reserve(n), rect.reserve(n); }
    void add(T lx, T ly, T rx, T ry) {
        rect.emplace_back(lx, ly, ry, 1);
        rect.emplace_back(rx, ly, ry, -1);
        ys.emplace_back(ly), ys.emplace_back(ry);
    }
    long long build() {
        if(empty(rect) or !ys.size()) return 0;
        sort(begin(ys), end(ys));
        ys.erase(unique(begin(ys), end(ys)), end(ys));
        sort(begin(rect), end(rect));
        T pre = get<0>(rect[0]);
        vector<S<T>> setter(ys.size() - 1);
        for(int i = 0; i < ys.size() - 1; i++) setter[i] = S<T>{ys[i + 1] - ys[i], 0, ys[i + 1] - ys[i]};
        lazy_segtree<S<T>, f<T>, e<T>, int, g<T>, h, id> seg(setter);
        long long res = 0;
        T len = ys.back() - ys.front();
        for(auto &&[x, l, r, v] : rect) {
            if(x != pre) {
                res += (long long)(x - pre) * (seg.all_prod().mi ? len : len - seg.all_prod().mi_cnt);
                pre = x;
            }
            seg.apply(lower_bound(begin(ys), end(ys), l) - begin(ys), lower_bound(begin(ys), end(ys), r) - begin(ys), v);
        }
        return res;
    }
};
} // namespace Rectangle_Union_impl
using Rectangle_Union_impl::Rectangle_Union;
