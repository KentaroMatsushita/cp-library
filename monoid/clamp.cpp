
namespace clamp_monoid_impl {
template <typename T> using Clamp = pair<T, T>;
template <typename T> Clamp<T> e() { return Clamp<T>(-inf<T>, inf<T>); }
template <typename T> Clamp<T> op(Clamp<T> x, Clamp<T> y) {
    Clamp<T> res;
    res.fi = (y.fi >= x.fi ? y.fi : x.fi > y.se ? y.se : x.fi);
    res.se = (y.fi >= x.se ? y.fi : x.se > y.se ? y.se : x.se);
    return res;
}
template <typename T> using segtree_clamp = segtree<Clamp<T>, op<T>, e<T>>;
} // namespace clamp_monoid_impl
using clamp_monoid_impl::Clamp;
using clamp_monoid_impl::segtree_clamp;

template <typename T> Clamp<T> e() { return Clamp<T>(-1, -1); }
template <typename T> Clamp<T> op(Clamp<T> x, Clamp<T> y) {
    if(x == e<T>()) return y;
    if(y == e<T>()) return x;
    return Clamp<T>(min(x.fi, y.fi), max(x.se, y.se));
}
template <typename T> Clamp<T> mapping(Clamp<T> f, Clamp<T> x) {
    if(x == e<T>()) return x;
    return clamp_monoid_impl::op(x, f);
}

template <typename T>
using RangeClamp_RangeOR_seg = lazy_segtree<Clamp<T>, op<T>, e<T>, Clamp<T>, mapping<T>, clamp_monoid_impl::op<T>, clamp_monoid_impl::e<T>>;
