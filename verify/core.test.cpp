#include <bits/stdc++.h>
using namespace std;
using ll=long long;
using Graph=vector<vector<int>>;
#define rep(i,n) for(int i=0;i<int(n);++i)
#define all(a) (a).begin(),(a).end()
#define UNIQUE(a) do {sort(all(a)); (a).erase(unique(all(a)),(a).end());}while(0)
#include "data_structure/UnionFind.cpp"
#include "data_structure/BIT.cpp"
#include "data_structure/SparseTable.cpp"
#include "data_structure/SegmentTree.cpp"
#include "dp/lis.cpp"
#include "Strings/z_algorithm.cpp"
#include "graph/topological_sort.cpp"
#include "graph/scc.cpp"
#include "graph/template.cpp"
#include "tree/diameter_UnWeighted.cpp"
void check(bool b) {if(!b) throw runtime_error("core test failed");}
int main() {
    auto merge=[](int& a,int b){a+=b;};
    UnionFind_value<int,decltype(merge)> u(vector<int>{2,3,7},merge);
    check(u.unite(0,1)); check(u.size(1)==2 && u.num==2 && u.value(0)==5);
    check(!u.unite(1,0)); check(u.num==2); u.reset(1);check(u.num==3&&u.size(0)==1);
    mt19937 rng(165871);
    for(int t=0;t<1000;++t) {
        vector<int>a(rng()%80);for(int& x:a)x=int(rng()%21)-10;
        BinaryIndexedTree<int> bit(a.size()); RST<int> seg(a);
        SparseTable<int> st(a,[](int a,int b){return min(a,b);},INT_MAX);
        rep(i,a.size())bit.add(i,a[i]);
        rep(k,100) {
            int l=rng()%(a.size()+1),r=rng()%(a.size()+1);if(l>r)swap(l,r);
            check(seg.prod(l,r)==accumulate(a.begin()+l,a.begin()+r,0));
            check(bit.sum(r-1)-bit.sum(l-1)==seg.prod(l,r));
            check(st.query(l,r)==(l==r?INT_MAX:*min_element(a.begin()+l,a.begin()+r)));
        }
        for(bool strict:{false,true}) {
            vector<int> dp(a.size(),1);int ans=0;
            rep(i,a.size()){rep(j,i)if(strict?a[j]<a[i]:a[j]<=a[i])dp[i]=max(dp[i],dp[j]+1);ans=max(ans,dp[i]);}
            auto idx=lis_restore(a,strict);check(int(idx.size())==ans && lis(a,strict)==ans);
            rep(i,idx.size()) if(i){check(idx[i-1]<idx[i]);check(strict?a[idx[i-1]]<a[idx[i]]:a[idx[i-1]]<=a[idx[i]]);}
        }
        auto z=z_algorithm(a);
        rep(i,a.size()){int k=0;while(i+k<int(a.size())&&a[k]==a[i+k])++k;check(z[i]==k);}
        string s;for(int x:a)s+=char('a'+(x+10)%4);
        z=z_algorithm(s);rep(i,s.size()){int k=0;while(i+k<int(s.size())&&s[k]==s[i+k])++k;check(z[i]==k);}
    }
    SCC<> scc(3);scc.add(0,1);scc.build();check(scc.num==3);
    scc.add(1,0);scc.build();check(scc.num==2 && scc[0]==scc[1]);
    Graph tree{{1},{0,2},{1}};Diameter d(tree);check(d.build()==2);check(d.build()==2&&d.path.size()==3);
    Graph empty;Diameter de(empty);check(de.build()==0&&de.path.empty());
    Graph dag{{1,2},{2},{}};auto ord=topological_sort(dag);check(ord==vector<int>({0,1,2}));
    istringstream input("1 2 5 2 3 7");auto* old=cin.rdbuf(input.rdbuf());
    auto weighted=make2<int>(3,2);cin.rdbuf(old);check(weighted.size()==3&&weighted[1].size()==2);
    cout<<"core: 1000 randomized arrays and regressions passed\n";
}
