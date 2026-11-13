// from https://judge.yosupo.jp/submission/5147
vector<int> prime_sieve(const int N, const int Q = 17, const int L = 1 << 15) {
    using u8 = unsigned char;
    static const int rs[] = {1, 7, 11, 13, 17, 19, 23, 29};
    struct P {
        P(int p) : p(p) {}
        int p;
        int pos[8];
    };
    auto approx_prime_count = [](const int N) -> int { return N > 60184 ? N / (log(N) - 1.1) : max(1., N / (log(N) - 1.11)) + 1; };

    const int v = sqrt(N), vv = sqrt(v);
    vector<bool> isp(v + 1, true);
    for(int i = 2; i <= vv; ++i)
        if(isp[i]) {
            for(int j = i * i; j <= v; j += i) isp[j] = false;
        }

    const int rsize = approx_prime_count(N + 30);
    vector<int> primes = {2, 3, 5};
    primes.reserve(rsize);

    vector<P> sprimes;
    size_t pbeg = 0;
    int prod = 1;
    for(int p = 7; p <= v; ++p) {
        if(!isp[p]) continue;
        if(p <= Q) prod *= p, ++pbeg, primes.push_back(p);
        auto pp = P(p);
        for(int t = 0; t < 8; ++t) {
            int j = (p <= Q) ? p : p * p;
            while(j % 30 != rs[t]) j += p << 1;
            pp.pos[t] = j / 30;
        }
        sprimes.push_back(pp);
    }

    vector<u8> pre(prod, 0xFF);
    for(size_t pi = 0; pi < pbeg; ++pi) {
        auto pp = sprimes[pi];
        const int p = pp.p;
        for(int t = 0; t < 8; ++t) {
            const u8 m = ~(1 << t);
            for(int i = pp.pos[t]; i < prod; i += p) pre[i] &= m;
        }
    }

    const int block_size = (L + prod - 1) / prod * prod;
    vector<u8> block(block_size);
    const int M = (N + 29) / 30;

    for(int beg = 0; beg < M; beg += block_size) {
        int end = min(M, beg + block_size);
        for(int i = beg; i < end; i += prod) { copy(pre.begin(), pre.end(), block.data() + (i-beg)); }
        if(beg == 0) block[0] &= 0xFE;
        for(size_t pi = pbeg; pi < sprimes.size(); ++pi) {
            auto &pp = sprimes[pi];
            const int p = pp.p;
            for(int t = 0; t < 8; ++t) {
                int i = pp.pos[t];
                const u8 m = ~(1 << t);
                for(; i < end; i += p) block[i-beg] &= m;
                pp.pos[t] = i;
            }
        }
        for(int i = beg; i < end; ++i) {
            for(int m = block[i-beg]; m > 0; m &= m - 1) { primes.push_back(i * 30 + rs[__builtin_ctz(m)]); }
        }
    }
    while(!primes.empty() && primes.back() > N) primes.pop_back();
    return primes;
}

vector<int> prime_table(int N) {
    auto r = prime_sieve(N);
    vector<int> res(N + 1);
    for(auto e : r) res[e] = 1;
    return res;
}
