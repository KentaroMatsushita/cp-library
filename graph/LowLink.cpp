#pragma once
struct LowLink {
  public:
    vector<int> ord, low, articulation;
    vector<pii> bridge;
    const Graph &g;
    virtual void build() {
        articulation.clear();
        bridge.clear();
        used.assign(g.size(), 0);
        ord.assign(g.size(), 0);
        low.assign(g.size(), 0);
        int k = 0;
        for(int i = 0; i < (int)g.size(); i++) {
            if(!used[i]) k = dfs(i, k, -1);
        }
    }

    LowLink(const Graph &g) : g(g) {}

  private:
    vector<int> used;

    int dfs(int idx, int k, int par) {
        used[idx] = true;
        ord[idx] = k++;
        low[idx] = ord[idx];
        bool is_articulation = false, beet = false;
        int cnt = 0;
        for(auto &to : g[idx]) {
            if(to == par && !exchange(beet, true)) { continue; }
            if(!used[to]) {
                ++cnt;
                k = dfs(to, k, idx);
                low[idx] = min(low[idx], low[to]);
                is_articulation |= par >= 0 && low[to] >= ord[idx];
                if(ord[idx] < low[to]) bridge.emplace_back(idx, to);
            } else {
                low[idx] = min(low[idx], ord[to]);
            }
        }
        is_articulation |= par == -1 && cnt > 1;
        if(is_articulation) articulation.push_back(idx);
        return k;
    }
};
