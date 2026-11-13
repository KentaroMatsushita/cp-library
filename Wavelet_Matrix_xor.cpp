//==================== bit_vector ====================
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

//====================================================
// 1) sumなし版：XOR-view Wavelet Matrix
//====================================================
template <class T, int lg = 31> struct WaveletMatrixXor {
    using u32 = uint32_t;
    using u64 = uint64_t;

    static_assert(lg > 0 && lg <= 63, "lg must be in [1,63]");

    int n = 0;
    array<bit_vector, lg> bv;
    u64 xor_mask = 0; // view mask (dataは変えない)

    WaveletMatrixXor() = default;
    explicit WaveletMatrixXor(const vector<T> &a) { build(a); }

    static constexpr u64 allmask() { return (1ULL << lg) - 1ULL; }

    inline void set_xor(u64 x) { xor_mask = x & allmask(); }

    __attribute__((optimize("O3"))) void build(const vector<T> &a) {
        n = (int)a.size();
        for(int h = 0; h < lg; ++h) bv[h] = bit_vector((u32)n);

        vector<u64> cur(n), nxt(n);
        for(int i = 0; i < n; ++i) cur[i] = (u64)a[i] & allmask();

        for(int h = lg - 1; h >= 0; --h) {
            for(int i = 0; i < n; ++i)
                if((cur[i] >> h) & 1ULL) bv[h].set((u32)i);
            bv[h].build();

            auto it0 = nxt.begin();
            auto it1 = nxt.begin() + bv[h].zeros;
            for(int i = 0; i < n; ++i) {
                if(bv[h].get((u32)i))
                    *it1++ = cur[i];
                else
                    *it0++ = cur[i];
            }
            cur.swap(nxt);
        }
        xor_mask = 0;
    }

    // (a[k]^xor_mask)
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
            u32 cnt0 = xbit ? (len - zeros) : zeros; // xor後bit=0の個数

            bool take1 = (k >= cnt0); // xor後bit=1側？
            if(take1) {
                k -= cnt0;
                res |= (1ULL << h);
            }

            // origbit = xorbit ^ xbit
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
        return (T)res; // xor後の値
    }

    T kth_largest(u32 l, u32 r, u32 k) const { return kth_smallest(l, r, (r - l) - k - 1); }

    // count (a[i]^xor_mask) < upper
    int range_freq(int l, int r, u64 upper) const {
        if(l >= r) return 0;
        if(upper >= (1ULL << lg)) return r - l;

        const u64 x = xor_mask;
        int ret = 0;
        u32 L = (u32)l, R = (u32)r;

        for(int h = lg - 1; h >= 0; --h) {
            u32 l0 = bv[h].rank0(L);
            u32 r0 = bv[h].rank0(R);
            u32 zeros = r0 - l0;
            u32 len = R - L;

            bool xbit = (x >> h) & 1ULL;
            bool ubit = (upper >> h) & 1ULL;
            u32 cnt0 = xbit ? (len - zeros) : zeros;

            if(ubit) {
                ret += (int)cnt0;
                bool origbit = xbit ^ 1ULL; // xor後bit=1側へ
                if(!origbit) {
                    L = l0;
                    R = r0;
                } else {
                    u32 Z = bv[h].zeros;
                    L = L + Z - l0;
                    R = R + Z - r0;
                }
            } else {
                bool origbit = xbit; // xor後bit=0側へ
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

    int range_freq(int l, int r, u64 lower, u64 upper) const {
        if(upper <= lower) return 0;
        return range_freq(l, r, upper) - range_freq(l, r, lower);
    }

    T prev_value(int l, int r, u64 upper) const {
        int cnt = range_freq(l, r, upper);
        return cnt == 0 ? T(-1) : kth_smallest((u32)l, (u32)r, (u32)(cnt - 1));
    }

    T next_value(int l, int r, u64 lower) const {
        int cnt = range_freq(l, r, lower);
        return cnt == r - l ? T(-1) : kth_smallest((u32)l, (u32)r, (u32)cnt);
    }
};

//====================================================
// 2) sumあり版：XOR-view + XOR後の和（bottom/top-k）を O(lg) で
//
// 実装方針：
//   各レベル ord[h] 上の区間 [l,r) に対し
//     sum(v^x) = sum(v) + len*x - 2*sum(v&x)
//   を使う。
//   sum(v&x) を “ニブル(4bit)×ブロックprefix” で O(lg/4) で計算。
//   bottom_k_sum は wavelet の分岐で取るノード区間を丸ごと sum(v^x) 加算。
//====================================================
template <class T, class sum_t = long long, int lg = 31, int BLK = 256, int NBITS = 4> struct WaveletMatrixXorSum {
    using u32 = uint32_t;
    using u64 = uint64_t;

    static_assert(lg > 0 && lg <= 63, "lg must be in [1,63]");
    static_assert(NBITS > 0 && NBITS <= 6, "NBITS should be small (e.g. 4)");
    static constexpr int NMASK = 1 << NBITS;
    static constexpr int P = (lg + NBITS - 1) / NBITS;

    int n = 0;
    int nb = 0;

    vector<u64> a;                  // original values (masked)
    array<bit_vector, lg> bv;       // wavelet bits
    vector<vector<u64>> ord;        // ord[h]
    vector<vector<sum_t>> pref_and; // per level, block-prefix for each (p,mask)

    // total sum on original order
    vector<sum_t> acc;
    array<bit_vector, lg> orig_bv;

    u64 xor_mask = 0;

    WaveletMatrixXorSum() = default;
    explicit WaveletMatrixXorSum(const vector<T> &_a) { build(_a); }

    static constexpr u64 allmask() { return (1ULL << lg) - 1ULL; }
    inline void set_xor(u64 x) { xor_mask = x & allmask(); }

    __attribute__((optimize("O3"))) void build(const vector<T> &_a) {
        n = (int)_a.size();
        nb = (n + BLK - 1) / BLK;

        a.resize(n);
        for(int i = 0; i < n; ++i) a[i] = (u64)_a[i] & allmask();

        // original-order prefix sum
        acc.assign(n + 1, 0);
        for(int i = 0; i < n; ++i) acc[i + 1] = acc[i] + (sum_t)a[i];

        // original-order bitvectors
        for(int b = 0; b < lg; ++b) {
            orig_bv[b].init((u32)n);
            for(int i = 0; i < n; ++i)
                if((a[i] >> b) & 1ULL) orig_bv[b].set((u32)i);
            orig_bv[b].build();
        }

        for(int h = 0; h < lg; ++h) bv[h] = bit_vector((u32)n);
        ord.assign(lg, {});
        pref_and.assign(lg, {});

        vector<u64> cur = a, nxt(n);

        static uint8_t and_table[NMASK][NMASK];
        static bool inited = false;
        if(!inited) {
            for(int m = 0; m < NMASK; ++m)
                for(int v = 0; v < NMASK; ++v) and_table[m][v] = (uint8_t)(m & v);
            inited = true;
        }

        for(int h = lg - 1; h >= 0; --h) {
            for(int i = 0; i < n; ++i)
                if((cur[i] >> h) & 1ULL) bv[h].set((u32)i);
            bv[h].build();

            auto it0 = nxt.begin();
            auto it1 = nxt.begin() + bv[h].zeros;
            for(int i = 0; i < n; ++i) {
                if(bv[h].get((u32)i))
                    *it1++ = cur[i];
                else
                    *it0++ = cur[i];
            }
            cur.swap(nxt);
            ord[h] = cur;

            // pref_and[h]
            pref_and[h].assign((size_t)P * NMASK * (nb + 1), 0);

            for(int b = 0; b < nb; ++b) {
                int L = b * BLK;
                int R = min(n, (b + 1) * BLK);

                uint16_t freq[P][NMASK];
                for(int p = 0; p < P; ++p)
                    for(int v = 0; v < NMASK; ++v) freq[p][v] = 0;

                for(int i = L; i < R; ++i) {
                    u64 v = cur[i];
                    for(int p = 0; p < P; ++p) freq[p][(v >> (NBITS * p)) & (NMASK - 1)]++;
                }

                for(int p = 0; p < P; ++p) {
                    for(int m = 0; m < NMASK; ++m) {
                        sum_t s = 0;
                        for(int v = 0; v < NMASK; ++v) s += (sum_t)freq[p][v] * (sum_t)and_table[m][v];
                        s <<= (NBITS * p);

                        size_t base = ((size_t)p * NMASK + (size_t)m) * (nb + 1);
                        pref_and[h][base + (b + 1)] = pref_and[h][base + b] + s;
                    }
                }
            }
        }

        xor_mask = 0;
    }

    // sum(ord[h][i] & mask), i in [l,r)
    inline sum_t range_and_level(int h, u32 l, u32 r, u64 mask) const {
        if(l >= r) return 0;
        mask &= allmask();

        u32 bl = l / BLK;
        u32 br = r / BLK;
        sum_t res = 0;

        if(bl == br) {
            for(u32 i = l; i < r; ++i) res += (sum_t)(ord[h][i] & mask);
            return res;
        }

        u32 Lend = min<u32>(r, (bl + 1) * BLK);
        for(u32 i = l; i < Lend; ++i) res += (sum_t)(ord[h][i] & mask);

        u32 Rstart = max<u32>(l, br * BLK);
        for(u32 i = Rstart; i < r; ++i) res += (sum_t)(ord[h][i] & mask);

        if(bl + 1 < br) {
            for(int p = 0; p < P; ++p) {
                int m = (int)((mask >> (NBITS * p)) & (NMASK - 1));
                size_t base = ((size_t)p * NMASK + (size_t)m) * (nb + 1);
                res += pref_and[h][base + br] - pref_and[h][base + (bl + 1)];
            }
        }
        return res;
    }

    // sum(ord[h][i] ^ xor_mask), i in [l,r)
    inline sum_t range_xor_sum_level(int h, u32 l, u32 r) const {
        if(l >= r) return 0;
        u64 x = xor_mask;
        u32 len = r - l;

        sum_t sum_v = range_and_level(h, l, r, allmask()); // == sum(v)
        sum_t sum_and = range_and_level(h, l, r, x);       // == sum(v&x)

        __int128 ans = (__int128)sum_v + (__int128)len * (__int128)x - (__int128)2 * (__int128)sum_and;
        return (sum_t)ans;
    }

    // Σ_{i in [l,r)} (a[i] ^ xor_mask)  (original order)
    sum_t range_sum_xor(int l, int r) const {
        if(l >= r) return 0;
        u64 x = xor_mask;
        sum_t sum = acc[r] - acc[l];
        sum_t len = (sum_t)(r - l);

        __int128 ans = (__int128)sum + (__int128)len * (__int128)x;
        u64 m = x;
        while(m) {
            int b = __builtin_ctzll(m);
            m &= (m - 1);
            u32 ones = orig_bv[b].rank1((u32)r) - orig_bv[b].rank1((u32)l);
            ans -= (__int128)2 * (__int128)ones * ((__int128)1 << b);
        }
        return (sum_t)ans;
    }

    // kth_smallest (xor view)
    u64 kth_smallest(u32 l, u32 r, u32 k) const {
        u64 res = 0;
        const u64 x = xor_mask;

        for(int h = lg - 1; h >= 0; --h) {
            u32 l0 = bv[h].rank0(l);
            u32 r0 = bv[h].rank0(r);
            u32 zeros = r0 - l0;
            u32 len = r - l;

            bool xbit = (x >> h) & 1ULL;
            u32 cnt0 = xbit ? (len - zeros) : zeros;

            bool take1 = (k >= cnt0);
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
        return res;
    }

    // count (a[i]^xor) < upper
    int range_freq(int l, int r, u64 upper) const {
        if(l >= r) return 0;
        if(upper >= (1ULL << lg)) return r - l;

        const u64 x = xor_mask;
        int ret = 0;
        u32 L = (u32)l, R = (u32)r;

        for(int h = lg - 1; h >= 0; --h) {
            u32 l0 = bv[h].rank0(L);
            u32 r0 = bv[h].rank0(R);
            u32 zeros = r0 - l0;
            u32 len = R - L;

            bool xbit = (x >> h) & 1ULL;
            bool ubit = (upper >> h) & 1ULL;
            u32 cnt0 = xbit ? (len - zeros) : zeros;

            if(ubit) {
                ret += (int)cnt0;
                bool origbit = xbit ^ 1ULL;
                if(!origbit) {
                    L = l0;
                    R = r0;
                } else {
                    u32 Z = bv[h].zeros;
                    L = L + Z - l0;
                    R = R + Z - r0;
                }
            } else {
                bool origbit = xbit;
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

    // XOR-view: sum of k smallest values in a[l,r) (k is count)
    sum_t bottom_k_sum(int l, int r, int k) const {
        int len = r - l;
        if(k <= 0) return 0;
        if(k >= len) return range_sum_xor(l, r);

        u32 L = (u32)l, R = (u32)r;
        sum_t ans = 0;
        u64 x = xor_mask;

        for(int h = lg - 1; h >= 0; --h) {
            u32 l0 = bv[h].rank0(L);
            u32 r0 = bv[h].rank0(R);
            u32 zeros = r0 - l0;
            u32 all = R - L;
            u32 ones = all - zeros;

            u32 Z = bv[h].zeros;
            u32 l1 = L + Z - l0;
            u32 r1 = R + Z - r0;

            bool xbit = (x >> h) & 1ULL;

            // xor後bit=0側（小さい側）の区間
            u32 sl, sr, bl, br, sz0;
            if(!xbit) { // xor後0はorig0
                sz0 = zeros;
                sl = l0;
                sr = r0;
                bl = l1;
                br = r1;
            } else { // xor後0はorig1
                sz0 = ones;
                sl = l1;
                sr = r1;
                bl = l0;
                br = r0;
            }

            if((u32)k < sz0) {
                L = sl;
                R = sr;
            } else {
                ans += range_xor_sum_level(h, sl, sr);
                k -= (int)sz0;
                L = bl;
                R = br;
            }
        }

        ans += range_xor_sum_level(0, L, L + (u32)k);
        return ans;
    }

    // XOR-view: sum of k largest values
    sum_t top_k_sum(int l, int r, int k) const {
        int len = r - l;
        if(k <= 0) return 0;
        if(k >= len) return range_sum_xor(l, r);
        return range_sum_xor(l, r) - bottom_k_sum(l, r, len - k);
    }
};

//====================================================
// 3) mex：永続 implicit segtree（last occurrence の min）+ XOR-view
//
// mex([l,r)) = 最小 y>=0 s.t. y が {a[i]^x} に存在しない
//           <=> v=y^x が [l,r) に存在しない
//           <=> last_pos_r[v] < l
//
// 各値 v の last_pos（prefix r まで）を持つ永続木を作り、
// “min last_pos < l の枝が存在するか” で mex をビットから決める。
//====================================================
template <int lg = 31> struct RangeMexXor {
    using u64 = uint64_t;
    static_assert(lg > 0 && lg <= 63, "lg must be in [1,63]");

    struct Node {
        int lch, rch, mn;
    }; // mn = min last_pos in segment

    static constexpr u64 allmask() { return (1ULL << lg) - 1ULL; }

    vector<Node> st;
    vector<int> root; // root[i] = prefix [0,i)
    u64 xor_mask = 0;

    RangeMexXor() = default;
    explicit RangeMexXor(const vector<u64> &a) { build(a); }

    void set_xor(u64 x) { xor_mask = x & allmask(); }

    int new_node(int lch, int rch, int mn) {
        st.push_back({lch, rch, mn});
        return (int)st.size() - 1;
    }

    int update(int node, int bit, u64 val, int pos) {
        if(bit < 0) return new_node(0, 0, pos);

        int lch = st[node].lch;
        int rch = st[node].rch;
        if((val >> bit) & 1ULL)
            rch = update(rch, bit - 1, val, pos);
        else
            lch = update(lch, bit - 1, val, pos);

        int mn = min(st[lch].mn, st[rch].mn);
        return new_node(lch, rch, mn);
    }

    void build(const vector<u64> &a) {
        int n = (int)a.size();
        st.clear();
        st.reserve((size_t)n * (lg + 2) + 5);

        // null node: mn=-1 (never appeared)
        st.push_back({0, 0, -1});

        root.assign(n + 1, 0);
        for(int i = 0; i < n; ++i) {
            u64 v = a[i] & allmask();
            root[i + 1] = update(root[i], lg - 1, v, i);
        }
        xor_mask = 0;
    }

    // mex of [l,r) under current xor_mask
    u64 mex(int l, int r) const {
        int node = root[r];
        if(st[node].mn >= l) return 1ULL << lg;
        u64 res = 0;
        u64 x = xor_mask;

        for(int bit = lg - 1; bit >= 0; --bit) {
            // y_bit=0 を優先 => v_bit = x_bit
            int prefer_v_bit = (int)((x >> bit) & 1ULL);
            int pref = prefer_v_bit ? st[node].rch : st[node].lch;
            int oth = prefer_v_bit ? st[node].lch : st[node].rch;

            // pref側に “欠けてる値” (last_pos < l) があるならそちらへ
            if(st[pref].mn < l) {
                node = pref; // y_bit=0
            } else {
                node = oth;
                res |= (1ULL << bit); // y_bit=1
            }
        }
        return res;
    }
};
