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
template <class S, S (*op)(S, S), S (*e)(), class F, S (*mapping)(S, F), F (*composition)(F, F), F (*id)()> struct lazy_segtree {
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
#ifdef ONLINE_JUDGE
        if(l < 0) l = 0;
        if(r > _n) r = _n;
        if(l > r) return e();
#endif
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
        d[p] = mapping(d[p], f);
        for(int i = 1; i <= log; i++) update(p >> i);
    }
    void apply(int l, int r, F f) {
#ifdef ONLINE_JUDGE
        if(l < 0) l = 0;
        if(r > _n) r = _n;
        if(l > r) return;
#endif
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
#ifdef ONLINE_JUDGE
        if(l < 0) l = 0;
        if(l > _n) l = _n;
#endif
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
#ifdef ONLINE_JUDGE
        if(r < 0) r = 0;
        if(r > _n) r = _n;
#endif
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

    int len() { return _n; }

  private:
    int _n, size, log;
    std::vector<S> d;
    std::vector<F> lz;

    void update(int k) { d[k] = op(d[2 * k], d[2 * k + 1]); }
    void all_apply(int k, F f) {
        d[k] = mapping(d[k], f);
        if(k < size) lz[k] = composition(lz[k], f);
    }
    void push(int k) {
        all_apply(2 * k, lz[k]);
        all_apply(2 * k + 1, lz[k]);
        lz[k] = id();
    }
};

template <class T> constexpr T inf = std::numeric_limits<T>::max() / 2;

template <class T, class U> constexpr std::pair<T, U> inf<std::pair<T, U>> = {inf<T>, inf<U>};

template <typename T> T min_func(T x, T y) { return min(x, y); }
template <typename T> T max_func(T x, T y) { return max(x, y); }
template <typename T> T plus_func(T x, T y) { return x + y; }
template <typename T, T ID> T set_func(T x, T y) { return (y == ID ? x : y); }
template <typename T, T ID> T set_id() { return ID; }
template <typename T> T min_e() { return inf<T>; }
template <typename T> T max_e() { return -inf<T>; }
template <typename T> T plus_e() { return T(); }

// range add get range min
template <typename T> using RARm = lazy_segtree<T, min_func<T>, min_e<T>, T, plus_func<T>, plus_func<T>, plus_e<T>>;
// range add get range max
template <typename T> using RARM = lazy_segtree<T, max_func<T>, max_e<T>, T, plus_func<T>, plus_func<T>, plus_e<T>>;

// range set get range max

template <typename T, T ID> using RSRM = lazy_segtree<T, max_func<T>, max_e<T>, T, set_func<T, ID>, set_func<T, ID>, set_id<T, ID>>;

template <typename T, T ID> using RSRm = lazy_segtree<T, min_func<T>, min_e<T>, T, set_func<T, ID>, set_func<T, ID>, set_id<T, ID>>;
// range add get range sum
// first が size, second が value
template <typename T> pair<int, T> RARS_op(pair<int, T> x, pair<int, T> y) { return {x.first + y.first, x.second + y.second}; }
template <typename T> pair<int, T> RARS_e() { return {0, 0}; }
template <typename T> pair<int, T> RARS_map(pair<int, T> x, T f) { return {x.first, x.second + f * x.first}; }

template <typename T> using _RARS = lazy_segtree<pair<int, T>, RARS_op<T>, RARS_e<T>, T, RARS_map<T>, plus_func<T>, plus_e<T>>;
template <typename T> struct RARS : _RARS<T> {
  private:
    vector<pair<int, T>> to_p(const vector<T> &v) {
        vector<pair<int, T>> res(v.size());
        for(int i = 0; i < v.size(); i++) res[i] = pair(1, v[i]);
        return res;
    }

  public:
    RARS(int n) : _RARS<T>(vector(n, pair(1, T()))) {}
    RARS(const vector<T> &v) : _RARS<T>(to_p(v)) {}
    T get(int k) { return _RARS<T>::get(k).se; }
    T prod(int l, int r) { return _RARS<T>::prod(l, r).se; }
    T operator[](int k) { return get(k); }
    friend std::ostream &operator<<(std::ostream &os, RARS &ls) {
        os << "{";
        for(int i = 0; i < ls.len(); i++) os << ls.get(i) << (i == ls.len() ? "" : ", ");
        return os << "}";
    }
};

// 区間の挿入、削除、union のサイズ取得
pair<int, int> segunion_merge(pair<int, int> x, pair<int, int> y) { return (x.se == y.se ? pair<int, int>{x.fi + y.fi, x.se} : x.se < y.se ? x : y); }
pair<int, int> segunion_e() { return {inf<int>, inf<int>}; }
pair<int, int> segunion_mapp(pair<int, int> p, int x) { return {p.fi, p.se + x}; }
using SegmentUnion = lazy_segtree<pair<int, int>, segunion_merge, segunion_e, int, segunion_mapp, plus_func<int>, plus_e<int>>;
} // namespace lazysegtree_impl

using lazysegtree_impl::lazy_segtree;
using lazysegtree_impl::RARM;
using lazysegtree_impl::RARm;
using lazysegtree_impl::RARS;
using lazysegtree_impl::RSRM;
using lazysegtree_impl::RSRm;
using lazysegtree_impl::SegmentUnion;