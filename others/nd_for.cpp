namespace ndFOR_impl {
template <int N> struct ndFORarray {
    std::array<int, N> v;
    ndFORarray(std::array<int, N> v_) : v(v_) {}
    struct ndFORitr {
        const std::array<int, N> &v;
        std::array<int, N> tmp;
        bool is_end;
        ndFORitr(const std::array<int, N> &v_) : v(v_), tmp(), is_end(std::any_of(v.begin(),v.end(),[](int x){return x<=0;})) {}
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
        ndFORitr(const std::vector<int> &v_) : v(v_), tmp(v.size(), 0), is_end(std::any_of(v.begin(),v.end(),[](int x){return x<=0;})) {}
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
