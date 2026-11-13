// description : [l,r) -> [L,R) に辺を貼ることができる
// 辺の本数は N log N 頂点の個数は + 2 * N
// 元の頂点の番号は sz + N で取得できる
template <typename T = ll> struct SectionGraph {
    int sz;
    WeightedGraph<T> g;
    SectionGraph(int n) {
        sz = 1;
        while(sz < n) sz <<= 1;
        g.resize(sz * 3);
        rep2(i, 1, sz - 1) {
            g[i].eb(i << 1, 0);
            g[i].eb((i << 1) + 1, 0);
            int nxt = i << 1;
            if(nxt < sz) nxt += sz << 1;
            g[nxt].eb(i + (sz << 1), 0);
            g[nxt | 1].eb(i + (sz << 1), 0);
        }
    }
    // [l,r) から [L,R) に辺を貼る
    void add_edge(int l, int r, int L, int R, T w) {
        if(l >= r || L >= R) return;
        T source_weight = (R - L == 1 ? w : T{});
        if(r - l == 1 and R - L == 1) {
            g[sz + l].eb(sz + L, w);
            return;
        }
        int id = g.size();
        if(r - l == 1) {
            id = sz + l;
        } else if(R - L == 1) {
            id = sz + L;
        } else
            g.eb(Edges<T>());
        if(R - L > 1)
            for(int d = R - L; d >= 1; d = R - L) {
                int Lt = L | ((1U << 31) >> __builtin_clz(d));
                int k = __builtin_ctz(Lt);
                g[id].eb(((sz | L) >> k), w);
                L += Lt & (-Lt);
            }
        if(r - l > 1)
            for(int d = r - l; d >= 1; d = r - l) {
                int L = l | ((1U << 31) >> __builtin_clz(d));
                int k = __builtin_ctz(L);
                int node = (sz | l) >> k;
                if(node < sz) node += sz * 2;
                g[node].eb(id, source_weight);
                l += L & (-L);
            }
    }

    void print() {
        rep(i, g.size()) {
            cout << i << " : ";
            for(auto e : g[i]) cout << e << " ";
            cout << endl;
        }
    }
};
