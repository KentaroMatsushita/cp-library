uint64_t floor_root(uint64_t a, uint64_t k) {
    if(k == 0) throw std::invalid_argument("zeroth root is undefined");
    if(a <= 1 || k == 1) return a;
    if(k >= 64) return 1;
    auto fits = [&](uint64_t x) {
        uint64_t product = 1;
        for(uint64_t i=0;i<k;++i) {
            if(product > a / x) return false;
            product *= x;
        }
        return true;
    };
    uint64_t l=1,r=uint64_t(1)<<((64+k-1)/k);
    while(r-l>1) {uint64_t m=l+(r-l)/2; if(fits(m)) l=m; else r=m;}
    return l;
}
