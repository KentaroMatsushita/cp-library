template <typename T> void fwt(vector<T> &f) {
    int n = f.size();
    for(int i = 1; i < n; i <<= 1) {
        for(int j = 0; j < n; j++) {
            if((j & i) == 0) { f[j | i] += f[j]; }
        }
    }
}
template <typename T> void ifwt(vector<T> &f) {
    int n = f.size();
    for(int i = 1; i < n; i <<= 1) {
        for(int j = 0; j < n; j++) {
            if((j & i) == 0) { f[j | i] -= f[j]; }
        }
    }
}

template <typename T> vector<T> or_conv(vector<T> f, vector<T> g) {
    int n = 1;
    while(n < max(f.size(), g.size())) n <<= 1;
    f.resize(n), g.resize(n);
    fwt(f), fwt(g);
    for(int i = 0; i < n; ++i) f[i] *= g[i];
    ifwt(f);
    return f;
}