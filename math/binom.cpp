long long C(int n, int m) {
    if(n < 0 || m < 0 || m > n) return 0;
    m = min(m, n - m);
    __int128_t res = 1;
    for(int i = 0; i < m; i++) {
        res *= (n - i);
        res /= (i + 1);
        if(res > std::numeric_limits<long long>::max()) throw std::overflow_error("binomial coefficient exceeds long long");
    }
    return res;
}
