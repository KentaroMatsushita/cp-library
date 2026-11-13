struct Subsetsum_restore {
    vector<int> v;
    int n;
    Subsetsum_restore(const vi &v) : v(v) {
        n = 0;
        while((1 << n) < si(v)) n++;
        dp.resize(1 << n + 1);
        build();
    }
    vector<vmint> dp;
    void build() {
        int N = 1 << n;
        rep(i, 1 << n) {
            if(i < si(v)) {
                dp[N + i].resize(v[i] + 1);
                dp[N + i][0] = dp[N + i].back() = rnd(998244353);
            } else {
                dp[N + i] = vmint{1};
            }
        }
        per(i, 1 << n, 1) { dp[i] = dp[i * 2] * dp[i * 2 + 1]; }
    }

    vector<int> sums() {
        vector<int> res;
        rep(i, dp[1].size()) {
            if(dp[1][i].a) res.eb(i);
        }
        return res;
    }

    pair<bool, vector<int>> get(int s) {
        if(s < 0 || s >= int(dp[1].size())) return {false, {}};
        if(dp[1][s].a == 0) return {false, {}};
        vector<int> tmp(1 << n + 1);
        tmp[1] = s;
        rep(i, 1, 1 << n) {
            int l = i * 2, r = i * 2 + 1;
            rep(x, dp[l].size()) {
                int y = tmp[i] - x;
                if(dp[l][x].a and y >= 0 and y < int(dp[r].size()) and dp[r][y].a) {
                    tmp[l] = x, tmp[r] = y;
                    break;
                }
            }
        }
        vector<int> res;
        rep(i, v.size()) if(tmp[i + (1 << n)]) res.emplace_back(i);
        return {true, res};
    }
};

/* less memory, slow
struct Subsetsum_restore {
    vector<int> v;
    int n;
    Subsetsum_restore(const vi &v) : v(v) {}

    vector<int> sums() {
        auto f = calc(0, v.size());
        vector<int> res;
        rep(i, si(f)) if(f[i].a) res.eb(i);
        return res;
    }

    vector<int> get(int s) {
        vector<int> res;
        _get(0, v.size(), s, res);
        return res;
    }

  private:
    void _get(int l, int r, int s, vector<int> &res) {
        if(s == 0) return;
        if(l + 1 == r) {
            assert(v[l] == s);
            res.emplace_back(l);
            return;
        }
        int mid = l + r >> 1;
        auto fl = calc(l, mid, s), fr = calc(mid, r, s);
        rep(i, si(fl)) {
            int j = s - i;
            if(fl[i].a and j < si(fr) and fr[j].a) {
                _get(l, mid, i, res);
                _get(mid, r, j, res);
                return;
            }
        }
    }
    vector<mint> calc(int l, int r, int mx = -1) {
        if(l + 1 == r) {
            if(mx != -1 and v[l] > mx) return vmint{1};
            vmint res(v[l] + 1);
            res[0] = res.back() = 1;
            return res;
        }
        int mid = l + r >> 1;
        auto res = calc(l, mid, mx) * calc(mid, r, mx);
        if(mx != -1 and si(res) > mx + 1) res.resize(mx + 1);
        return res;
    }
};


*/
