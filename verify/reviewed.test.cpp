#include <bits/stdc++.h>
using namespace std;
using Graph=vector<vector<int>>;
using pii=pair<int,int>;
using vi=vector<int>;
#define rep(i,n) for(int i=0;i<int(n);++i)
#define si(a) int((a).size())
#define all(a) (a).begin(),(a).end()
#define UNIQUE(a) do{sort(all(a));(a).erase(unique(all(a)),(a).end());}while(0)
template<class T>bool chmax(T& a,T b){if(a<b){a=b;return true;}return false;}
#include "graph/template.cpp"
template<class T>using Wgraph=WeightedGraph<T>;
#include "graph/dijkstra_restore.cpp"
#include "graph/2-SAT.cpp"
#include "graph/two-edge-connected-components.cpp"
#include "tree/center.cpp"
#include "Strings/LCS.cpp"
#include "dp/convex-convolution.hpp"
#include "data_structure/UnionFind_rollback.cpp"
namespace weighted {
#include "data_structure/UnionFind_weighted.cpp"
}
namespace potential {
#include "data_structure/UnionFind_potential.cpp"
}
void check(bool b){if(!b)throw runtime_error("review test failed");}
int main(){
    check(LCS("axb","ayb")==2);
    mt19937 rng(148712);
    for(int t=0;t<1000;++t){
        string a,b;for(int i=rng()%12;i--;)a+=char('a'+rng()%3);
        for(int i=rng()%12;i--;)b+=char('a'+rng()%3);
        int best=0;for(int mask=0;mask<(1<<a.size());++mask){
            string s;rep(i,a.size())if(mask>>i&1)s+=a[i];
            size_t k=0;for(char c:b)if(k<s.size()&&s[k]==c)++k;
            if(k==s.size())best=max(best,int(k));
        }
        check(LCS(a,b)==best);
        int n=1+rng()%20,m=1+rng()%20;
        vector<int>x(n),y(m);int slope=-10;
        rep(i,n){slope+=rng()%4;x[i]=(i?x[i-1]:0)+slope;}
        for(int& v:y)v=int(rng()%100)-50;
        auto z=convex_min_conv(x,y);rep(k,z.size()){
            int ref=INT_MAX;rep(i,n)if(k>=i&&k-i<m)ref=min(ref,x[i]+y[k-i]);check(z[k]==ref);
        }
        for(int& v:x)v=-v;z=convex_max_conv(x,y);rep(k,z.size()){
            int ref=INT_MIN;rep(i,n)if(k>=i&&k-i<m)ref=max(ref,x[i]+y[k-i]);check(z[k]==ref);
        }
        TwoSat sat(4);vector<pair<int,int>> clauses;
        rep(k,12){int a=1+rng()%4,b=1+rng()%4;if(rng()%2)a=-a;if(rng()%2)b=-b;clauses.emplace_back(a,b);sat.OR(a,b);}
        int valid=0;auto accepts=[&](int mask){for(auto [a,b]:clauses){auto value=[&](int x){return bool(mask>>(abs(x)-1)&1)==(x>0);};if(!value(a)&&!value(b))return false;}return true;};
        rep(mask,16)valid+=accepts(mask);check(sat.can()==(valid>0));auto ans=sat.solve();
        check(ans.empty()==(valid==0));if(!ans.empty()){int mask=0;rep(i,4)mask|=ans[i]<<i;check(accepts(mask));}
    }
    check(monotone_minima(0,0,[](int,int,int){return false;}).empty());
    check(monotone_minima(3,0,[](int,int,int){return false;})==vector<int>(3,-1));
    UnionFind_undo undo(3);undo.unite(0,1);undo.unite(0,1);undo.rollback(2);check(undo.num==3&&!undo.same(0,1));
    auto merge=[](int& a,int b){a+=b;};weighted::UnionFind<int,decltype(merge)> u(4,merge,1);
    u.unite(0,1);u.unite(2,3);u.unite(1,3);check(u.size(0)==4&&u.num==1&&u.get(0)==4);
    potential::UnionFind<long long> p(5);p.merge(0,1,3000000000LL);p.merge(2,3,2);p.merge(4,2,7);
    // Rank-swapped merge must negate the relative potential.
    check(p.val(0,1)==-3000000000LL);check(p.val(4,3)==-9);
    Graph g{{1},{0,2},{1}};TwoEdgeConnectedComponents c(g);check(c.bridge.size()==2);c.build();check(c.bridge.size()==2&&c.group.size()==3);
    check(search_center(g)==vector<int>{1});Graph single(1);check(search_center(single)==vector<int>{0});
    Graph even{{1},{0,2},{1,3},{2}};auto centers=search_center(even);sort(all(centers));check(centers==vector<int>({1,2}));
    Wgraph<int> wg(1);auto [path,d]=dijkstra(wg,0,0);check(d==0&&path==vector<int>{0});
    cout<<"reviewed: LCS/convolution/2SAT differential tests and graph/DSU regressions passed\n";
}
