#include "scc.cpp"

struct TwoSat : SCC<Graph> {
    int n;
    using SCC::SCC;
    TwoSat(int n) : n(n), SCC(n * 2) {}
    // not i = i + n
    inline int rev(int x) { return x >= n ? x - n : x + n; }
    inline int id(int x) { return x < 0 ? n - x - 1 : x - 1; }
    inline void IF(int x, int y) {
        // x => y
        add(x, y), add(rev(y), rev(x));
    }
    void OR(int x, int y) { IF(rev(id(x)), id(y)); }
    vector<int> solve() {
        build();
        vector<int> res(n);
        rep(i, n) {
            if(comp[i] == comp[rev(i)]) return vector<int>();
            res[i] = comp[i] > comp[rev(i)];
        }
        return res;
    }
    bool can() {
        build();
        rep(i, n) {
            if(comp[i] == comp[rev(i)]) return false;
        }
        return true;
    }
};