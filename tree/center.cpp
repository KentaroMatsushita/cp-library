// Unweighted connected tree; O(N), no recursion. Empty tree => {}.
template <typename G> vector<int> search_center(const G& g) {
    if(g.empty()) return {};
    auto bfs=[&](int s,vector<int>& par) {
        par.assign(g.size(),-1);
        vector<int> d(g.size(),-1);d[s]=0;
        queue<int> q;q.push(s);int far=s;
        while(!q.empty()) {
            int x=q.front();q.pop();
            if(d[x]>d[far]) far=x;
            for(int y:g[x]) if(d[y]<0) {d[y]=d[x]+1;par[y]=x;q.push(y);}
        }
        return far;
    };
    vector<int> par;int a=bfs(0,par),b=bfs(a,par);
    vector<int> path;
    for(int x=b;x!=-1;x=par[x]) path.push_back(x);
    int n=path.size();
    if(n%2) return {path[n/2]};
    return {path[n/2-1],path[n/2]};
}
