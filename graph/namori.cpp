template <typename G> struct Namori {
    int n;
    G g;
    // 一番近いサイクル上の頂点は何か
    vector<int> root;
    // サイクル
    vector<vector<int>> cycles;
    // cycle 上の辺でないなら -1, そうなら属するサイクルの番号
    vector<int> cycle_num;
    // 森, サイクルから遠ざかる方向の辺だけ
    vector<vector<int>> forest;
    // g から サイクル上の辺を取り除いたもの
    G gg;
    // 親, サイクル上なら -1
    vector<int> par;
    Namori(int n) : n(n), g(n), cycle_num(n, -1), root(n), forest(n), par(n, -1) {}
    Namori(G g) : n(size(g)), g(g), cycle_num(size(g), -1), root(size(g)), forest(size(g)), par(size(g), -1) {}
    void add_edge(int x, int y) {
        g[x].emplace_back(y);
        g[y].emplace_back(x);
    }

    void build() {
        vector<int> cnt(n);
        vector<int> used(n);
        rep(i, n) for(auto e : g[i]) cnt[e]++;
        queue<int> q;
        rep(i, n) if(cnt[i] == 1) {
            q.emplace(i);
            used[i] = true;
        }
        while(!empty(q)) {
            int x = q.front();
            q.pop();
            int p;
            for(auto e : g[x]) {
                if(used[e]) continue;
                p = e;
                forest[p].emplace_back(x);
                break;
            }
            par[x] = p;
            if((--cnt[p]) == 1) {
                used[p] = true;
                q.emplace(p);
            }
        }
        rep(i, n) {
            if(!used[i]) {
                vi cycle;
                int now = i;
                used[i] = true;
                cycle.emplace_back(now);
                cycle_num[now] = cycles.size();
                while(1) {
                    bool ch = false;
                    for(auto e : g[now]) {
                        if(!used[e]) {
                            now = e;
                            ch = true;
                            break;
                        }
                    }
                    if(!ch) break;
                    used[now] = true;
                    cycle.emplace_back(now);
                    cycle_num[now] = cycles.size();
                }
                cycles.emplace_back(cycle);
            }
        }
        for(int r = 0; r < n; r++) {
            if(cycle_num[r] != -1) {
                queue<int> q({r});
                while(!empty(q)) {
                    int x = q.front();
                    root[x] = r;
                    q.pop();
                    for(auto e : forest[x]) { q.emplace(e); }
                }
            }
        }
    }

    const int operator[](int k) const { return root[k]; }

    // サイクル上で i -> i + 1 のコストを並べた vector を返す
    template <class T> vector<vector<T>> get_cycle() {
        vector<vector<T>> RES;
        for(auto &&cycle : cycles) {
            vector<T> res(cycle.size(), -1);
            bool flag = (cycle.size() == 2);
            if(flag) {
                vector<T> v;
                for(auto e : g[cycle[0]]) {
                    if(e == cycle[1]) v.emplace_back(e.cost);
                }
                res[0] = *min_element(v.begin(), v.end());
                res[1] = *max_element(v.begin(), v.end());
            } else {
                for(int i = 0; i < cycle.size(); i++) {
                    int l = cycle[i], r = cycle[(i + 1) % cycle.size()];
                    for(auto e : g[l]) {
                        if(e == r) res[i] = e.cost;
                    }
                }
            }
            RES.emplace_back(res);
        }
        return RES;
    }

    // サイクルから離れる方向だけの重み付きグラフを作る
    template <class T> Wgraph<T> get_weighted_forest() {
        Wgraph<T> res(n);
        queue<int> q;
        vector<int> vis(n);
        rep(i, n) if(cycle_num[i] != -1) q.emplace(i), vis[i] = true;
        while(!empty(q)) {
            int x = q.front();
            q.pop();
            for(auto e : g[x]) {
                if(vis[e]) continue;
                res[x].emplace_back(e, e.cost);
                q.emplace(e);
                vis[e] = true;
            }
        }
        return res;
    }
};