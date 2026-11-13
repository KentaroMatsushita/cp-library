template <class T = int> struct Imos {
    int n;
    vector<T> a;
    Imos(int _n) : n(_n), a(_n + 1) {}
    void add(int l, int r, T val = 1) {
        if(l >= r) return;
        l = clamp(l, 0, n);
        r = clamp(r, 0, n + 1);
        a[l] += val;
        if(r <= n) a[r] -= val;
    }
    void build() {
        for(int i = 0; i < n; i++) a[i + 1] += a[i];
    }
    const T &operator[](int k) { return a[k]; }
};