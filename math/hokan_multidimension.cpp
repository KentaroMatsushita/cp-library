namespace hokan_multidim_impl {
// {rank, det(非正方行列の場合は未定義)} を返す
// 型が double や Rational でも動くはず？(未検証)
//
// pivot 候補 : [0, pivot_end)
template <typename T> std::pair<int, T> GaussElimination(vector<vector<T>> &a, int pivot_end = -1, bool diagonalize = false) {
    int H = a.size(), W = a[0].size(), rank = 0;
    if(pivot_end == -1) pivot_end = W;
    T det = 1;
    for(int j = 0; j < pivot_end; j++) {
        int idx = -1;
        for(int i = rank; i < H; i++) {
            if(a[i][j] != T(0)) {
                idx = i;
                break;
            }
        }
        if(idx == -1) {
            det = 0;
            continue;
        }
        if(rank != idx) det = -det, swap(a[rank], a[idx]);
        det *= a[rank][j];
        if(diagonalize && a[rank][j] != T(1)) {
            T coeff = T(1) / a[rank][j];
            for(int k = j; k < W; k++) a[rank][k] *= coeff;
        }
        int is = diagonalize ? 0 : rank + 1;
        for(int i = is; i < H; i++) {
            if(i == rank) continue;
            if(a[i][j] != T(0)) {
                T coeff = a[i][j] / a[rank][j];
                for(int k = j; k < W; k++) a[i][k] -= a[rank][k] * coeff;
            }
        }
        rank++;
    }
    return make_pair(rank, det);
}

// 解が存在する場合は, 解が v + C_1 w_1 + ... + C_k w_k と表せるとして
// (v, w_1, ..., w_k) を返す
// 解が存在しない場合は空のベクトルを返す
//
// double や Rational でも動くはず？(未検証)
template <typename T> vector<vector<T>> LinearEquation(vector<vector<T>> a, vector<T> b) {
    int H = a.size(), W = a[0].size();
    for(int i = 0; i < H; i++) a[i].push_back(b[i]);
    auto p = GaussElimination(a, W, true);
    int rank = p.first;
    for(int i = rank; i < H; ++i) {
        if(a[i][W] != 0) return vector<vector<T>>{};
    }
    vector<vector<T>> res(1, vector<T>(W));
    vector<int> pivot(W, -1);
    for(int i = 0, j = 0; i < rank; ++i) {
        while(a[i][j] == 0) ++j;
        res[0][j] = a[i][W], pivot[j] = i;
    }
    for(int j = 0; j < W; ++j) {
        if(pivot[j] == -1) {
            vector<T> x(W);
            x[j] = 1;
            for(int k = 0; k < j; ++k) {
                if(pivot[k] != -1) x[k] = -a[pivot[k]][j];
            }
            res.push_back(x);
        }
    }
    return res;
}

namespace ndFOR_impl {
template <int N> struct ndFORarray {
    std::array<int, N> v;
    ndFORarray(std::array<int, N> v_) : v(v_) {}
    struct ndFORitr {
        const std::array<int, N> &v;
        std::array<int, N> tmp;
        bool is_end;
        ndFORitr(const std::array<int, N> &v_) : v(v_), tmp(), is_end(false) {}
        bool operator!=(const ndFORitr &) const { return !is_end; }
        void operator++() {
            int pos = N - 1;
            while(pos != -1) {
                tmp[pos] += 1;
                if(tmp[pos] == v[pos]) {
                    tmp[pos] = 0;
                    pos -= 1;
                } else {
                    break;
                }
            }
            if(pos == -1) { is_end = true; }
        }
        const std::array<int, N> &operator*() const { return tmp; }
    };
    ndFORitr begin() const { return ndFORitr(v); }
    ndFORitr end() const { return ndFORitr(v); }
};

struct ndFORvector {
    std::vector<int> v;
    ndFORvector(std::vector<int> v_) : v(v_) {}
    struct ndFORitr {
        const std::vector<int> &v;
        std::vector<int> tmp;
        bool is_end;
        ndFORitr(const std::vector<int> &v_) : v(v_), tmp(v.size(), 0), is_end(false) {}
        bool operator!=(const ndFORitr &) const { return !is_end; }
        void operator++() {
            int pos = v.size() - 1;
            while(pos != -1) {
                tmp[pos] += 1;
                if(tmp[pos] == v[pos]) {
                    tmp[pos] = 0;
                    pos -= 1;
                } else {
                    break;
                }
            }
            if(pos == -1) { is_end = true; }
        }
        const std::vector<int> &operator*() const { return tmp; }
    };
    ndFORitr begin() const { return ndFORitr(v); }
    ndFORitr end() const { return ndFORitr(v); }
};

// sum, length (if any, -1), limit [l, r)
struct FixedSum_vec {
    int s, n, lim_l, lim_r;
    FixedSum_vec(int s, int n = -1, int lim_l = 0, int lim_r = INT_MAX) : s(s), n(n), lim_l(lim_l), lim_r(lim_r > s ? s + 1 : lim_r) {
        assert(n != -1 or lim_l > 0);
        assert(lim_l < lim_r);
    }
    struct itr {
        vector<int> tmp;
        int s, n, lim_l, lim_r, tsum;
        bool is_end;
        itr(int s, int n, int lim_l, int lim_r) : s(s), n(n), lim_l(lim_l), lim_r(lim_r), tsum(0), is_end(false) {
            if(n == -1)
                tmp.reserve(s / lim_l);
            else { tmp.reserve(n); }
            if(s == 0 and (n == 0 or n == -1)) {
            } else
                next();
        }
        // -1 : short, 0 : ok, 1 : large
        inline int can_fill() {
            if(!empty(tmp) and (tmp.back() == lim_r or tsum > s)) return 1;
            if(n != -1) {
                if((n - tmp.size()) * lim_l > s - tsum) return 1;
                if(s - tsum > (n - tmp.size()) * (lim_r - 1)) return -1;
            }
            if(lim_l) {
                int need = (s - tsum) / lim_l;
                if((s - tsum) > need * (lim_r - 1)) return -1;
            }
            return 0;
        }
        bool operator!=(const itr &) const { return !is_end; }
        void next() {
            while(!empty(tmp)) {
                tmp.back()++;
                tsum++;
                if(int k = can_fill(); k == 1) {
                    tsum -= tmp.back();
                    tmp.pop_back();
                    if(empty(tmp)) {
                        is_end = true;
                        return;
                    }
                } else if(k == 0)
                    break;
            }
            if(can_fill() != 0) {
                is_end = true;
                return;
            }
            while((tsum < s) or ((int)tmp.size() < n)) {
                tmp.emplace_back(lim_l);
                tsum += lim_l;
                if(can_fill() == -1) {
                    if(n != -1) {
                        int need = (s - tsum - (n - tmp.size()) * (lim_r - 1));
                        tmp.back() += need;
                        tsum += need;
                    } else {
                        assert(lim_l);
                        int r = (s - tsum) / lim_l;
                        int need = (s - tsum) - (lim_r - 1) * r;
                        tmp.back() += need;
                        tsum += need;
                    }
                }
            }
        }
        void operator++() {
            if(s == 0 and (n == 0 or n == -1))
                is_end = true;
            else
                next();
        }
        const std::vector<int> &operator*() const { return tmp; }
    };
    itr begin() const { return itr(s, n, lim_l, lim_r); }
    itr end() const { return itr(s, n, lim_l, lim_r); }
};
auto ndFOR(std::vector<int> v) { return ndFORvector(v); }
auto FixedSumVec(int s, int n = -1, int lim_l = 0, int lim_r = INT_MAX) { return FixedSum_vec(s, n, lim_l, lim_r); }
template <class... Ts> auto ndFOR(Ts... v) { return ndFORarray<std::tuple_size<std::tuple<Ts...>>::value>({v...}); }
}; // namespace ndFOR_impl
using ndFOR_impl::FixedSumVec;
using ndFOR_impl::ndFOR;

// v[i] := 列, w[i] := v[i] のときの答え
vector<pair<vector<int>, mint>> hokan_multidim(vector<vector<mint>> v, vector<mint> w, int maxk = 10) {
    // 次数を小さい方から試す
    vector<vi> f;
    rep(k, maxk + 1) {
        fore(w, FixedSumVec(k, si(v[0]))) { f.emplace_back(w); }
        vv(mint, A, si(v), si(f));
        vmint B(si(v));
        rep(i, si(v)) {
            rep(j, si(f)) {
                mint c = 1;
                rep(k, si(v[i])) c *= v[i][k].pow(f[j][k]);
                A[i][j] = c;
            }
            B[i] = w[i];
        }
        auto res = LinearEquation(A, B);
        if(si(res) == 0) continue;
        if(si(res) == 1) {
            vector<pair<vi, mint>> ans;
            rep(i, si(f)) ans.eb(f[i], res[0][i]);
            return ans;
        } else {
            assert(false); // 等式の数不足
        }
    }
    assert(false); // 次数不足または多項式でないまたは愚直ミス
}
} // namespace hokan_multidim_impl

using hokan_multidim_impl::hokan_multidim;