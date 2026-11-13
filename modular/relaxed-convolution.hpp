template <typename mint> struct RelaxedConvolution {
    mint push(const mint &x, const mint &y) {
        ++n;
        f.emplace_back(x), g.emplace_back(y);
        same &= (x == y);
        for(int p = 1;; p <<= 1) {
            int l1 = n - p, r1 = n, l2 = p - 1, r2 = l2 + p;
            if(same) {
                add(l1, r1, l2, r2, (l1 == l2 ? 1 : 2));
                if(l1 == l2) break;
            } else {
                add(l1, r1, l2, r2);
                if(l1 == l2) break;
                add(l2, r2, l1, r1);
            }
            if(!(n & p)) break;
        }
        return h[n - 1];
    }
    const mint &operator[](int i) const { return h[i]; }
    vector<mint> res() const { return h; }

  private:
    bool same = true;
    int n = 0;
    vector<mint> f, g, h;
    void add(int l1, int r1, int l2, int r2, int coef = 1) {
        auto k = vector(f.begin() + l1, f.begin() + r1) * vector(g.begin() + l2, g.begin() + r2);
        int st = l1 + l2;
        if((st + k.size()) > h.size()) h.resize(r1 + r2 - 1);
        for(int i = 0; i < k.size(); i++) { h[st + i] += (coef == 1 ? k[i] : k[i] * coef); }
    }
};

// ABC 315-ex