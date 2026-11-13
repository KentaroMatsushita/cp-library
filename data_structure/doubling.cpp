template <typename T = int> struct Doubling {
    int n, LOG;
    vector<vector<int>> nxt;
    vector<vector<T>> cost;
    Doubling(int n) : n(n) {
        LOG = n ? 32 - __builtin_clz(n) + 1 : 1;
        nxt.assign(LOG, vector<int>(n));
        cost.assign(LOG, vector<T>(n, T{}));
        for(int i=0;i<LOG;++i) iota(nxt[i].begin(),nxt[i].end(),0);
    }
    Doubling(const vector<int> &_nxt, const vector<T> &_cost, int lg = 0) {
        n = _nxt.size();
        LOG = (lg ? lg : (n ? 32 - __builtin_clz(n) + 1 : 1));
        nxt.assign(LOG, vector<int>(n)), cost.assign(LOG, vector<T>(n));
        nxt[0] = _nxt, cost[0] = _cost;
        rep(i, 1, LOG) rep(j, n) {
            nxt[i][j] = nxt[i - 1][nxt[i - 1][j]];
            cost[i][j] = cost[i - 1][j] + cost[i - 1][nxt[i - 1][j]];
        }
    }
    Doubling(const vector<int> &_nxt) : Doubling(_nxt, vi(si(_nxt), 1)) {}

    pair<int, T> move(int from, T times) {
        if(times < 0) throw std::invalid_argument("negative move count");
        int to = from;
        T res = 0;
        for(int i = 0; times; i++, times >>= 1) {
            if(i >= LOG) throw std::out_of_range("move exceeds precomputed doubling levels");
            if(times & 1) {
                res += cost[i][to];
                to = nxt[i][to];
            }
        }
        return {to, res};
    }

    bool can_move(int l, int r) {
        per(i, LOG, 0) {
            if(nxt[i][l] <= r) { l = nxt[i][l]; }
        }
        return l == r;
    }

    // l から r 以上まで行くのに必要なコスト
    T move_cnt(int l, int r) {
        if(l >= r) return 0;
        T res = 0;
        per(i, LOG, 0) {
            if(nxt[i][l] < r) { res += cost[i][l], l = nxt[i][l]; }
        }
        res += cost[0][l];
        return res;
    }

    // f を満たす限り進むときコストいくつ分進めて、最後はどこか
    // template <typename F> pair<T, int> move_cnt(int l, const F &f) {
    //     assert(f(l));
    //     T res = 0;
    //     per(i, LOG, 0) {
    //         if(f(nxt[i][l])) res += cost[i][l], l = nxt[i][l];
    //     }
    //     return {res, l};
    // }
};
