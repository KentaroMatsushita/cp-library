
// f(i) = 1^i + ... + (k-1)^i を n 乗まで
vector<mint> kpowsum(int n, ll k) {
    vmint f(n + 1), g(n + 1);
    auto ps = power_table(k, n + 2);
    rep(i, n + 1) f[i] = ps[i + 1] * ifact(i + 1);
    rep(i, n + 1) g[i] = ifact(i + 1);
    auto res = f * inv(g);
    res.resize(n + 1);
    rep(i, 1, n + 1) res[i] *= fact(i);
    res[0] = k - 1;
    return res;
}

// \sum_{i = 0} ^ {n - 1} f(i)
mint sum_eval(vmint f, int n) {
    if(n == 0 || f.empty()) return mint(0);
    mint res;
    auto kps = kpowsum(f.size(), n);
    rep(i, f.size()) { res += f[i] * kps[i]; }
    res += f[0]; // kpowsum excludes the evaluation at zero.
    return res;
}

// \sum_{i = L} ^ {R - 1} f(ix)
vmint kx_sum(vmint f, ll L, ll R) {
    vmint res(f.size());
    auto kps = kpowsum(f.size(), R) - kpowsum(f.size(), L);
    rep(i, f.size()) res[i] += f[i] * kps[i];
    return res;
}

// \prod_{i = L} ^ {R - 1} f(ix) を n 次まで
vmint kx_prod(vmint f, ll L, ll R, int n) {
    if(L == R) { vmint one(n + 1); one[0] = 1; return one; }
    mint c;
    int d = -1;
    rep(i, f.size()) {
        if(f[i].a) {
            d = i;
            c = f[i].inv();
            break;
        }
    }
    if(d == -1 or (R - L) * d > n) return {};
    mint t = f[d].pow(R - L);
    if(d) { rep(i, L, R) t *= mint(i).pow(d); }
    vmint g = log((f >> d) * c, n + 1);
    g = exp(kx_sum(g, L, R));
    g <<= (R - L) * d;
    g.resize(n + 1);
    return g * t;
}
