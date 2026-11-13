// works on Tree
// (u,v) に対して u の方が根に近い時、v を根とする部分木に対応する番号付けを行う
template <typename G> struct EulerTour {
    const G &g;
    EulerTour(G &g) : g(g) {}
    void dfs(int x, int p, int &num, vector<int> &in, vector<int> &out, int edgenum = -1) {
        if(edgenum != -1) in[edgenum] = num++;
        for(auto &e : g[x]) {
            if(e != p) dfs(e, x, num, in, out, e.id);
        }
        if(edgenum != -1) out[edgenum] = num;
    }
    pair<vector<int>, vector<int>> build(int root = 0) {
        int num = 0;
        vector<int> in(g.empty() ? 0 : g.size() - 1), out(in.size());
        if(!g.empty()) dfs(root, -1, num, in, out);
        return {in, out};
    }
};
