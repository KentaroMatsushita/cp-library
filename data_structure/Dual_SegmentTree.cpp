template <typename T> struct DualSegmentTree {
    using H = function<T(T, T)>;
    int sz, height;
    vector<T> data;
    const H h;
    const T e;
    DualSegmentTree(int n, const H h, T e) : h(h), e(e) {
        sz = 1, height = 0;
        while(sz < n) sz <<= 1, height++;
        data.assign(2 * sz, e);
    }
    inline void build(vector<T> &v) {
        for(int i = 0; i < v.size(); ++i) data[sz + i] = v[i];
    }
    inline void set(int k, T x) { data[sz + k] = x; }
    DualSegmentTree(vector<T> &v, const H h, T e) : h(h), e(e) {
        sz = 1, height = 0;
        while(sz < v.size()) sz <<= 1, height++;
        data.assign(2 * sz, e);
        build(v);
    }
    inline void propagate(int k) {
        if(data[k] != e) {
            data[2 * k + 0] = h(data[2 * k + 0], data[k]);
            data[2 * k + 1] = h(data[2 * k + 1], data[k]);
            data[k] = e;
        }
    }
    inline void thrust(int k) {
        for(int i = height; i > 0; i--) propagate(k >> i);
    }
    void update(int a, int b, const T &x) {
        if(a >= b) return;
        thrust(a += sz), thrust(b += sz - 1);
        for(int l = a, r = b + 1; l < r; l >>= 1, r >>= 1) {
            if(l & 1) data[l] = h(data[l], x), ++l;
            if(r & 1) --r, data[r] = h(data[r], x);
        }
    }
    T operator[](int k) {
        thrust(k += sz);
        return data[k];
    }
    friend ostream &operator<<(ostream &os, DualSegmentTree seg) {
        for(int i = 0; i < seg.sz; ++i) { os << seg[i] << " "; }
        return os;
    }
};

// Range Add Tree
template <class T> struct RAT : DualSegmentTree<T> {
    using Seg = DualSegmentTree<T>;
    RAT(size_t n) : Seg(n, plus<T>(), T()) {}
    RAT(vector<T> &a) : Seg(a, plus<T>(), T()) {}
};

// Range Set Tree
// inf を逆元として利用
template <class T> struct RST : DualSegmentTree<T> {
    using Seg = DualSegmentTree<T>;
    RST(size_t n)
        : Seg(
              n, [](T a, T b) { return (b == numeric_limits<T>::max() ? a : b); }, numeric_limits<T>::max()) {}
    RST(vector<T> &a)
        : Seg(
              a, [](T a, T b) { return (b == numeric_limits<T>::max() ? a : b); }, numeric_limits<T>::max()) {}
};
