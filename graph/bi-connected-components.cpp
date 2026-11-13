
template <typename G> struct LowLink {
    const G &g;
    vector<int> used, ord, low;
    vector<int> articulation;
    vector<pair<int, int>> bridge;

    LowLink(const G &g) : g(g) {}

    int dfs(int idx, int k, int par) {
        used[idx] = true;
        ord[idx] = k++;
        low[idx] = ord[idx];
        bool is_articulation = false;
        int cnt = 0;
        for(auto &to : g[idx]) {
            if(!used[to]) {
                ++cnt;
                k = dfs(to, k, idx);
                low[idx] = min(low[idx], low[to]);
                is_articulation |= ~par && low[to] >= ord[idx];
                if(ord[idx] < low[to]) bridge.emplace_back(minmax(idx, (int)to));
            } else if(to != par) {
                low[idx] = min(low[idx], ord[to]);
            }
        }
        is_articulation |= par == -1 && cnt > 1;
        if(is_articulation) articulation.push_back(idx);
        return k;
    }

    virtual void build() {
        articulation.clear(); bridge.clear();
        used.assign(g.size(), 0);
        ord.assign(g.size(), 0);
        low.assign(g.size(), 0);
        int k = 0;
        for(int i = 0; i < g.size(); i++) {
            if(!used[i]) k = dfs(i, k, -1);
        }
    }
};
template <typename G> struct BiConnectedComponents : LowLink<G> {
    using LoL = LowLink<G>;

    vector<int> used;
    vector<vector<pair<int, int>>> bc;
    vector<pair<int, int>> tmp;

    BiConnectedComponents(const G &g) : LoL(g) {}

    void dfs(int idx, int par) {
        used[idx] = true;
        for(auto &to : this->g[idx]) {
            if(to == par) continue;
            if(!used[to] || this->ord[to] < this->ord[idx]) { tmp.emplace_back(minmax(idx, to)); }
            if(!used[to]) {
                dfs(to, idx);
                if(this->low[to] >= this->ord[idx]) {
                    bc.emplace_back();
                    for(;;) {
                        auto e = tmp.back();
                        bc.back().emplace_back(e);
                        tmp.pop_back();
                        if(e.first == min(idx, to) && e.second == max(idx, to)) { break; }
                    }
                }
            }
        }
    }

    void build() override {
        bc.clear(); tmp.clear();
        LoL::build();
        used.assign(this->g.size(), 0);
        for(int i = 0; i < used.size(); i++) {
            if(!used[i]) dfs(i, -1);
        }
    }
};
