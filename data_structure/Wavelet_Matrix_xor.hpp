struct bit_vector {
    using u32 = uint32_t;
    using u64 = uint64_t;

    static constexpr u32 w = 64;
    vector<u64> block;
    vector<u32> count;
    u32 n = 0, zeros = 0;

    bit_vector() = default;
    explicit bit_vector(u32 _n) { init(_n); }

    inline u32 get(u32 i) const { return (u32)((block[i / w] >> (i % w)) & 1ULL); }
    inline void set(u32 i) { block[i / w] |= (1ULL << (i % w)); }

    __attribute__((optimize("O3", "unroll-loops"))) void init(u32 _n) {
        n = zeros = _n;
        block.assign(n / w + 1, 0);
        count.assign(block.size(), 0);
    }

    __attribute__((target("popcnt"))) void build() {
        for(u32 i = 1; i < (u32)block.size(); ++i) count[i] = count[i - 1] + (u32)_mm_popcnt_u64(block[i - 1]);
        zeros = rank0(n);
    }

    inline u32 rank0(u32 i) const { return i - rank1(i); }

    __attribute__((target("bmi2,popcnt"))) inline u32 rank1(u32 i) const { return count[i / w] + (u32)_mm_popcnt_u64(_bzhi_u64(block[i / w], i % w)); }
};

// -------------------- XOR-view Wavelet Matrix (queries) --------------------
template <class T, int lg = 30> struct WaveletMatrix {
    static_assert(lg > 0 && lg < 64);
    using u32 = uint32_t;
    using u64 = uint64_t;

    int n = 0;
    vector<T> a;
    array<bit_vector, lg> bv;

    // XOR view mask (dataは変更しない)
    u64 xor_mask = 0;
    inline void set_xor(u64 x) {
        if(x >> lg) throw std::out_of_range("XOR mask exceeds Wavelet Matrix bit width");
        xor_mask = x;
    }
    inline void xor_all(u64 x) { set_xor(xor_mask ^ x); }

    WaveletMatrix() = default;
    explicit WaveletMatrix(u32 _n) : n((int)max<u32>(_n, 1)), a(n) {}
    explicit WaveletMatrix(const vector<T> &_a) : n((int)_a.size()), a(_a) { build(); }

    __attribute__((optimize("O3"))) void build() {
        for(int h = 0; h < lg; ++h) bv[h] = bit_vector((u32)n);

        vector<T> cur = a, nxt(n);
        for(int h = lg - 1; h >= 0; --h) {
            for(int i = 0; i < n; ++i)
                if(((u64)cur[i] >> h) & 1ULL) bv[h].set((u32)i);
            bv[h].build();

            auto it0 = begin(nxt);
            auto it1 = begin(nxt) + bv[h].zeros;
            for(int i = 0; i < n; ++i) {
                if(bv[h].get((u32)i))
                    *it1++ = cur[i];
                else
                    *it0++ = cur[i];
            }
            swap(cur, nxt);
        }
    }

    // (a[k] ^ xor_mask)
    T access(u32 k) const {
        u64 ret = 0;
        for(int h = lg - 1; h >= 0; --h) {
            u32 f = bv[h].get(k);
            if(f) ret |= (1ULL << h);
            k = f ? (bv[h].rank1(k) + bv[h].zeros) : bv[h].rank0(k);
        }
        return (T)(ret ^ xor_mask);
    }

    // k-th smallest (0-indexed) in (a[i]^xor_mask), i in [l,r)
    T kth_smallest(u32 l, u32 r, u32 k) const {
        u64 res = 0;
        const u64 x = xor_mask;

        for(int h = lg - 1; h >= 0; --h) {
            u32 l0 = bv[h].rank0(l);
            u32 r0 = bv[h].rank0(r);
            u32 zeros = r0 - l0;
            u32 len = r - l;

            bool xbit = (x >> h) & 1ULL;
            u32 cnt0 = xbit ? (len - zeros) : zeros; // XOR後bit=0の個数

            bool take1 = (k >= cnt0); // XOR後bit=1側へ？
            if(take1) {
                k -= cnt0;
                res |= (1ULL << h);
            }

            bool origbit = xbit ^ take1;
            if(!origbit) {
                l = l0;
                r = r0;
            } else {
                u32 Z = bv[h].zeros;
                l = l + Z - l0;
                r = r + Z - r0;
            }
        }
        return (T)res; // XOR後の値
    }

    T kth_largest(u32 l, u32 r, u32 k) const { return kth_smallest(l, r, (r - l) - k - 1); }

    // count i in [l,r) s.t. (a[i]^xor_mask) < upper
    int range_freq(int l, int r, T upper) const {
        if(l >= r) return 0;
        if(upper <= 0) return 0;
        u64 U = (u64)upper;
        if(U >= (1ULL << lg)) return r - l;

        const u64 x = xor_mask;
        int ret = 0;
        u32 L = (u32)l, R = (u32)r;

        for(int h = lg - 1; h >= 0; --h) {
            u32 l0 = bv[h].rank0(L);
            u32 r0 = bv[h].rank0(R);
            u32 zeros = r0 - l0;
            u32 len = R - L;

            bool xbit = (x >> h) & 1ULL;
            bool ubit = (U >> h) & 1ULL;

            u32 cnt0 = xbit ? (len - zeros) : zeros; // XOR後bit=0の個数

            if(ubit) {
                ret += (int)cnt0;
                bool origbit = xbit ^ 1ULL; // XOR後bit=1へ
                if(!origbit) {
                    L = l0;
                    R = r0;
                } else {
                    u32 Z = bv[h].zeros;
                    L = L + Z - l0;
                    R = R + Z - r0;
                }
            } else {
                bool origbit = xbit; // XOR後bit=0へ
                if(!origbit) {
                    L = l0;
                    R = r0;
                } else {
                    u32 Z = bv[h].zeros;
                    L = L + Z - l0;
                    R = R + Z - r0;
                }
            }
        }
        return ret;
    }

    int range_freq(int l, int r, T lower, T upper) const { return range_freq(l, r, upper) - range_freq(l, r, lower); }

    T prev_value(int l, int r, T upper) const {
        int cnt = range_freq(l, r, upper);
        return cnt == 0 ? T(-1) : kth_smallest((u32)l, (u32)r, (u32)(cnt - 1));
    }

    T next_value(int l, int r, T lower) const {
        int cnt = range_freq(l, r, lower);
        return cnt == r - l ? T(-1) : kth_smallest((u32)l, (u32)r, (u32)cnt);
    }
};

// -------------------- SUM + XOR-after-sum Wavelet Matrix --------------------
template <class T, class sum_t = long long, int lg = 31> struct WaveletMatrix_sum_xor : WaveletMatrix<T, lg> {
    using WM = WaveletMatrix<T, lg>;
    using u32 = uint32_t;
    using u64 = uint64_t;

    using WM::a;
    using WM::bv;
    using WM::n;
    using WM::xor_mask;

    // sums[h] : ord_h(=bit h 処理後の並び)での prefix sum（元の値）
    array<vector<sum_t>, lg> sums;

    // low_bv[h][b] : ord_h 上での「bit b」が 1 の位置bitvector（b < h のみ）
    // これでノード区間(=ord_h 上の contiguous interval)内の bit b の1数が O(1)
    array<vector<bit_vector>, lg> low_bv;

    // original order 上の bitvector（区間全体の XOR和に使う）
    array<bit_vector, lg> orig_bv;

    // original order prefix sums（元の値の区間和）
    vector<sum_t> acc;

    WaveletMatrix_sum_xor() = default;

    explicit WaveletMatrix_sum_xor(const vector<T> &_a) {
        n = (int)_a.size();
        a = _a;
        build_all();
    }

    void build() = delete;

    __attribute__((optimize("O3"))) void build_all() {
        // bv / sums / low_bv を作る（ord を回しながら）
        for(int h = 0; h < lg; ++h) {
            bv[h] = bit_vector((u32)n);
            sums[h].assign(n + 1, sum_t());
            low_bv[h].clear();
        }

        vector<T> cur = a, nxt(n);

        for(int h = lg - 1; h >= 0; --h) {
            // bv[h] は ord_{h+1} 上の bit h
            for(int i = 0; i < n; ++i)
                if(((u64)cur[i] >> h) & 1ULL) bv[h].set((u32)i);
            bv[h].build();

            // stable partition -> ord_h
            auto it0 = begin(nxt);
            auto it1 = begin(nxt) + bv[h].zeros;
            for(int i = 0; i < n; ++i) {
                if(bv[h].get((u32)i))
                    *it1++ = cur[i];
                else
                    *it0++ = cur[i];
            }
            swap(cur, nxt); // cur は ord_h

            // sums[h]（元の値）
            for(int i = 0; i < n; ++i) sums[h][i + 1] = sums[h][i] + (sum_t)cur[i];

            // low_bv[h][b] を ord_h 上で構築（b < h のみ）
            if(h > 0) {
                low_bv[h].resize(h);
                for(int b = 0; b < h; ++b) low_bv[h][b].init((u32)n);

                for(int i = 0; i < n; ++i) {
                    u64 v = (u64)cur[i];
                    // 下位 h bit を見て、立ってるbitだけ set
                    // (最悪 O(n*h) だけど lg=31なら許容)
                    for(int b = 0; b < h; ++b)
                        if((v >> b) & 1ULL) low_bv[h][b].set((u32)i);
                }
                for(int b = 0; b < h; ++b) low_bv[h][b].build();
            }
        }

        // original order prefix sum
        acc.assign(n + 1, sum_t());
        for(int i = 0; i < n; ++i) acc[i + 1] = acc[i] + (sum_t)a[i];

        // original order bitvectors
        for(int b = 0; b < lg; ++b) {
            orig_bv[b].init((u32)n);
            for(int i = 0; i < n; ++i)
                if(((u64)a[i] >> b) & 1ULL) orig_bv[b].set((u32)i);
            orig_bv[b].build();
        }
    }

    // Σ_{i in [l,r)} (a[i] ^ xor_mask)
    sum_t range_sum_xor(int l, int r) const {
        if(l >= r) return sum_t();
        sum_t sum = acc[r] - acc[l];
        u64 x = xor_mask;
        sum_t len = (sum_t)(r - l);

        // v^x = v + x - 2*(v&x)
        sum_t res = sum + len * (sum_t)x;

        u64 m = x;
        while(m) {
            int b = __builtin_ctzll(m);
            m &= (m - 1);
            u32 ones = orig_bv[b].rank1((u32)r) - orig_bv[b].rank1((u32)l);
            res -= (sum_t)2 * (sum_t)ones * ((sum_t)1 << b);
        }
        return res;
    }

    // このノード（bit h で分割した直後の ord_h 上の区間 [l,r)）を丸ごと足すときの
    // Σ (value ^ xor_mask) を返す
    //
    // prefix_above : 現在の親ノードの固定上位bit（value >> (h+1)）
    // origbit      : この子ノードの bit h（0/1）
    sum_t node_sum_xor(int h, u32 l, u32 r, u64 prefix_above, bool origbit) const {
        u32 len_u = r - l;
        if(len_u == 0) return sum_t();
        sum_t len = (sum_t)len_u;

        // このノードでは value>>h が一定：high = (prefix_above<<1)|origbit
        u64 high = (prefix_above << 1) | (u64)origbit;

        sum_t sum_orig = sums[h][r] - sums[h][l]; // 元の値の和

        // sum_orig = len*(high<<h) + sum_low
        sum_t base_high = len * (sum_t)(high << h);
        sum_t sum_low = sum_orig - base_high; // 下位 h bit の和（元）

        // XOR後の上位部分： (high ^ (x>>h)) << h が一定
        u64 x = xor_mask;
        u64 high_xored = high ^ (x >> h);
        sum_t res = len * (sum_t)(high_xored << h);

        // 下位 h bit 部分だけ XOR：sum_low_xor = sum_low + Σ_{b<h, x_b=1} (len-2*ones_b)*2^b
        sum_t sum_low_xor = sum_low;

        if(h > 0) {
            u64 x_low = x & ((1ULL << h) - 1ULL); // 下位 h bit
            u64 m = x_low;
            while(m) {
                int b = __builtin_ctzll(m);
                m &= (m - 1);
                // b < h なので low_bv[h][b] がある
                u32 ones = low_bv[h][b].rank1(r) - low_bv[h][b].rank1(l);
                sum_t delta = len - (sum_t)2 * (sum_t)ones;
                sum_low_xor += delta * ((sum_t)1 << b);
            }
        }

        res += sum_low_xor;
        return res;
    }

    // XOR後の「k個最小」の XOR後の和： Σ_{t=0..k-1} kth_smallest(t)   (kは個数)
    sum_t bottom_k_sum(int l, int r, int k) const {
        int len = r - l;
        if(k <= 0) return sum_t();
        if(k >= len) return range_sum_xor(l, r);

        u32 L = (u32)l, R = (u32)r;
        sum_t ans = 0;

        // prefix_above は value>>(h+1) を表す（最初は value>>lg = 0）
        u64 prefix_above = 0;
        u64 x = xor_mask;

        for(int h = lg - 1; h >= 0; --h) {
            // ここでの L,R は ord_{h+1} 上の区間
            u32 l0 = bv[h].rank0(L);
            u32 r0 = bv[h].rank0(R);
            u32 zeros = r0 - l0;
            u32 all = R - L;

            u32 Z = bv[h].zeros;
            u32 l1 = L + Z - l0;
            u32 r1 = R + Z - r0;

            bool xbit = (x >> h) & 1ULL;

            // XOR後 bit=0 側（小さい側）は origbit = xbit
            bool origbit_small = xbit;
            u32 size0 = origbit_small ? (all - zeros) : zeros;

            u32 ls = origbit_small ? l1 : l0;
            u32 rs = origbit_small ? r1 : r0;
            u32 ll = origbit_small ? l0 : l1;
            u32 rl = origbit_small ? r0 : r1;

            if((u32)k < size0) {
                // 小さい側へ降りる
                prefix_above = (prefix_above << 1) | (u64)origbit_small;
                L = ls;
                R = rs;
            } else {
                // 小さい側を丸ごと足す
                ans += node_sum_xor(h, ls, rs, prefix_above, origbit_small);
                k -= (int)size0;

                // 大きい側へ
                bool origbit_large = origbit_small ^ 1;
                prefix_above = (prefix_above << 1) | (u64)origbit_large;
                L = ll;
                R = rl;
            }
        }

        // ここまでで prefix_above は value>>0 = value（葉）
        // 残り k 個は同一値なので k*(value^x)
        ans += (sum_t)k * (sum_t)((u64)prefix_above ^ x);
        return ans;
    }

    // XOR後の「k個最大」の XOR後の和
    sum_t top_k_sum(int l, int r, int k) const {
        int len = r - l;
        if(k <= 0) return sum_t();
        if(k >= len) return range_sum_xor(l, r);
        // 最大k個 = 全体 - 最小(len-k)個
        return range_sum_xor(l, r) - bottom_k_sum(l, r, len - k);
    }
};

