struct Random {
  private:
    uint32_t a, b, c, d;

  public:
    Random(uint32_t seed = 517) { set_seed(seed); }
    void set_seed(uint32_t seed) {
        a = seed = 1812433253u * (seed ^ (seed >> 30));
        b = seed = 1812433253u * (seed ^ (seed >> 30)) + 1;
        c = seed = 1812433253u * (seed ^ (seed >> 30)) + 2;
        d = seed = 1812433253u * (seed ^ (seed >> 30)) + 3;
    }
    uint32_t operator()() {
        uint32_t t = (a ^ (a << 11));
        a = b;
        b = c;
        c = d;
        return d = (d ^ (d >> 19)) ^ (t ^ (t >> 8));
    }

    // [0,MAX)
    int operator()(int MAX) { return (uint64_t) operator()() * MAX >> 32; }

    // [L,R)
    uint32_t operator()(int L, int R) { return L + operator()(R - L); }

    // [0,1) の実数
    double drnd() { return (double)operator()() / (1LL << 32); }

    // 確率 p で true を返す
    bool prob(const double p) { return p >= 1.0 || (p > 0.0 && drnd() < p); }

    bool coin() { return operator()() < (1u << 31); }

    template <typename T> void shuffle(vector<T> &v) {
        for(int i = int(v.size()) - 1; i > 0; i--) { swap(v[i], v[operator()(i + 1)]); }
    }
} rnd;

// long long に対する murmur hash
struct murmur_hash_int64 {
    unsigned long long operator()(long long p) const {
        const unsigned long long m = 0xc6a4a7935bd1e995;
        auto x = static_cast<unsigned long long>(p) * m;
        unsigned long long h = (x ^ (x >> 47)) * m;
        return h = (h ^ (h >> 47)) * m, (h ^ (h >> 47));
    }
} hash64;
