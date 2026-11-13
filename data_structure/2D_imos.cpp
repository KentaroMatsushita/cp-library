template <typename T> struct imos2D {
    vector<vector<T>> v;
    imos2D(int H, int W) : v(H + 1, vector<T>(W + 1)) {}
    void add(int d, int l, int u, int r, T val) { v[d][l] += val, v[d][r] -= val, v[u][l] -= val, v[u][r] += val; }
    void build() {
        for(int i = 0; i < v.size(); i++)
            for(int j = 0; j < (int)v[0].size() - 1; j++) v[i][j + 1] += v[i][j];
        for(int i = 0; i < (int)v.size() - 1; i++)
            for(int j = 0; j < v[0].size(); j++) v[i + 1][j] += v[i][j];
    }
    inline const vector<T> &operator[](int i) const { return v[i]; }
};
