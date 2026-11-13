struct CycleDecomposisiton {
    int n;
    vector<pair<int, int>> v;
    CycleDecomposisiton(Graph &g) : n(g.size()), v(g.size()) {
        rep(i, n) {
            assert(g[i].size() == 2);
            v[i] = pair(g[i][0], g[i][1]);
        }
    }
    CycleDecomposisiton(int n) : n(n), v(n, pair(-1, -1)) {}

    void add(int x, int y) { _add(x, y), _add(y, x); }
    void add_directed(int x, int y) { _add(x, y); }

    vector<vector<int>> build() {
        vector<vector<int>> res;
        vector<int> used(n);
        vector<int> x;
        rep(i, n) {
            if(!used[i]) {
                x.clear();
                x.emplace_back(i);
                if(v[i].first != i) {
                    x.emplace_back(v[i].first);
                    while(true) {
                        if(x[x.size() - 2] == v[x.back()].first)
                            x.emplace_back(v[x.back()].second);
                        else
                            x.emplace_back(v[x.back()].first);
                        if(x.back() == i or x.back() == -1) {
                            x.pop_back();
                            break;
                        }
                    }
                }
                res.eb(x);
                for(auto e : x) used[e] = true;
            }
        }
        return res;
    }

  private:
    void _add(int x, int y) {
        if(v[x].first == -1)
            v[x].first = y;
        else if(v[x].second == -1)
            v[x].second = y;
        else
            assert(false);
    }
};