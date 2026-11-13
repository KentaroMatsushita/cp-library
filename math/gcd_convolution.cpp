namespace gcd_convolution {
template <typename T> void f(vector<T> &a) {
    int n = a.size();
    for(int k = 1; k < n; ++k) {
        for(int i = k * 2; i < n; i += k) a[k] += a[i];
    }
}
template <typename T> void invf(vector<T> &a) {
    int n = a.size();
    for(int k = n - 1; k >= 1; --k) {
        for(int i = k * 2; i < n; i += k) a[k] -= a[i];
    }
}
template <typename T> vector<T> mul(vector<T> a, vector<T> b) {
    f(a), f(b);
    rep(i, si(b)) a[i] *= b[i];
    invf(a);
    return a;
}
} // namespace gcd_convolution