template<typename Abel>
struct UnionFind {
    typedef function<Abel(Abel,Abel)> F;
    F plus,inv;
    Abel e;
    vector<int> par;
    vector<int> rank;
    vector<Abel> diff;
    UnionFind(int n=1,Abel s = Abel{}){
        init(n,s,std::plus<Abel>{},std::minus<Abel>{});
    }
    void init(int n,Abel _e,F _plus,F _inv){
        plus=_plus; inv=_inv; e=_e;
        par.resize(n);rank.resize(n);diff.resize(n);
        rep(i,n)par[i]=i,rank[i]=0,diff[i]=_e;
    }
    int root(int x){
        if(par[x]==x)return x;
        int r = root(par[x]);
        diff[x]=plus(diff[x],diff[par[x]]);
        return par[x]=r;
    }
    Abel weight(int x){
        root(x);
        return diff[x];
    }
    bool same(int x,int y){return root(x)==root(y);}
    bool merge(int x,int y,Abel w){
        w = plus(w,weight(x));
        w = inv(w,weight(y));
        x = root(x);y=root(y);
        if(x==y){
            return false;
        }
        if(rank[x]<rank[y]){swap(x,y);w=inv(e,w);}
        if(rank[x]==rank[y])++rank[x];
        par[y]=x;
        diff[y]=w;
        return true;
    }
    Abel val(int x,int y){
        if(!same(x,y))return -1;
        return inv(weight(x),weight(y));
    }
};