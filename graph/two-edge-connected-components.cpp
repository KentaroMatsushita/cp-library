#include "LowLink.cpp"

struct TwoEdgeConnectedComponents : LowLink {
  public:
    using LowLink::bridge;
    using LowLink::g;
    using LowLink::low;
    using LowLink::LowLink;
    using LowLink::ord;

    vector<int> comp;
    Graph tree;
    vector<vector<int>> group;

    int operator[](const int &k) const { return comp[k]; }

    void build() override {
        LowLink::build();
        comp.assign(g.size(), -1);
        int k = 0;
        for(int i = 0; i < (int)comp.size(); i++) {
            if(comp[i] == -1) dfs(i, -1, k);
        }
        group.assign(k, {});
        for(int i = 0; i < (int)g.size(); i++) { group[comp[i]].emplace_back(i); }
        tree = Graph(k);
        for(auto &e : bridge) {
            tree[comp[e.first]].emplace_back(comp[e.second]);
            tree[comp[e.second]].emplace_back(comp[e.first]);
        }
    }

    TwoEdgeConnectedComponents(const Graph &g) : LowLink(g) { build(); }

  private:
    void dfs(int idx, int par, int &k) {
        if(par >= 0 && ord[par] >= low[idx])
            comp[idx] = comp[par];
        else
            comp[idx] = k++;
        for(auto &to : g[idx]) {
            if(comp[to] == -1) dfs(to, idx, k);
        }
    }
};