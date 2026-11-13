// works on Tree
template <typename G> struct EulerTour {
    const G &g;
    EulerTour(G &g) : g(g) {}
    void dfs(int x, int p, int &num, vector<int> &in, vector<int> &out) {
        in[x] = num++;
        for(auto &e : g[x]) {
            if(e != p) dfs(e, x, num, in, out);
        }
        out[x] = num;
    }
    pair<vector<int>, vector<int>> build(int root = 0) {
        int num = 0;
        vector<int> in(g.size()), out(g.size());
        if(!g.empty()) dfs(root, -1, num, in, out);
        return {in, out};
    }
};
