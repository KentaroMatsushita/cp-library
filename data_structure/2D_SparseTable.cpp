template <typename T = int, const int LLOG = 12, bool ismin = true> struct SparseTable2D {
    vector<vector<T>> dat[LLOG][LLOG];
    int LOG;
    vector<int> lookup;
    int H, W;
    SparseTable2D(int H, int W) : H(H), W(W) {
        LOG = 0;
        while(1 << LOG <= max(H, W)) LOG++;
        if(LOG > LLOG) throw std::invalid_argument("SparseTable2D LLOG too small");
        for(int i = 0; i < LOG; i++)
            for(int j = 0; j < LOG; j++) dat[i][j].assign(H, vector<T>(W));
    }
    SparseTable2D(vector<vector<T>> &a) : SparseTable2D(int(a.size()), a.empty() ? 0 : int(a[0].size())) {
        for(int i = 0; i < a.size(); ++i)
            for(int j = 0; j < W; ++j) dat[0][0][i][j] = a[i][j];
        build();
    }
    void set(int i, int j, T x) { dat[0][0][i][j] = x; }
    void build() {
        lookup.resize(max(H, W) + 1);
        for(int i = 2; i < lookup.size(); ++i) lookup[i] = lookup[i >> 1] + 1;
        for(int k = 0; k < LOG; ++k) {
            for(int l = 0; l < LOG; ++l) {
                if(!k and !l) continue;
                for(int i = 0; i + (1 << k) <= H; ++i)
                    for(int j = 0; j + (1 << l) <= W; ++j) {
                        if(k) {
                            if(ismin)
                                dat[k][l][i][j] = min(dat[k - 1][l][i][j], dat[k - 1][l][i + (1 << (k - 1))][j]);
                            else
                                dat[k][l][i][j] = max(dat[k - 1][l][i][j], dat[k - 1][l][i + (1 << (k - 1))][j]);
                        } else {
                            if(ismin)
                                dat[k][l][i][j] = min(dat[k][l - 1][i][j], dat[k][l - 1][i][j + (1 << (l - 1))]);
                            else
                                dat[k][l][i][j] = max(dat[k][l - 1][i][j], dat[k][l - 1][i][j + (1 << (l - 1))]);
                        }
                    }
            }
        }
    }
    inline T query(int i, int j, int ii, int jj) {
        int b = lookup[ii - i], bb = lookup[jj - j];
        if(ismin)
            return min({dat[b][bb][i][j], dat[b][bb][ii - (1 << b)][j], dat[b][bb][i][jj - (1 << bb)], dat[b][bb][ii - (1 << b)][jj - (1 << bb)]});
        else
            return max({dat[b][bb][i][j], dat[b][bb][ii - (1 << b)][j], dat[b][bb][i][jj - (1 << bb)], dat[b][bb][ii - (1 << b)][jj - (1 << bb)]});
    }
};
