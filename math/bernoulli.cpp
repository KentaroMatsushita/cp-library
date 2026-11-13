Poly Bernoulli(int n) {
    Poly f(n + 1);
    rep(i, n + 1) f[i] = ifact(i + 1);
    f = inv(f);
    rep(i, n + 1) f[i] *= fact(i);
    return f;
}