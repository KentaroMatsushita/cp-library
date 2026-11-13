ll primitive_root() {
    if(MOD == 2) return 1;
    vector<int> primes;
    ll m = MOD - 1;
    for(int i = 2; i * i <= MOD; ++i) {
        if(m % i == 0) {
            while(m % i == 0) m /= i;
            primes.emplace_back(i);
        }
    }
    if(m > 1) primes.emplace_back(m);
    m = MOD - 1;
    mint g = 2;
    while(1) {
        bool is_pr = true;
        for(auto p : primes) {
            if((g ^ (m / p)) == mint(1)) is_pr = false;
        }
        if(is_pr) return g.a;
        g += 1;
    }
}
