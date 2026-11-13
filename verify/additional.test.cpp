#include <bits/stdc++.h>
using namespace std;
using ll=long long;
using vi=vector<int>;
using Graph=vector<vector<int>>;
#define REP2(i,n) for(int i=0;i<int(n);++i)
#define REP3(i,l,r) for(int i=int(l);i<int(r);++i)
#define SELECT_REP(_1,_2,_3,NAME,...) NAME
#define rep(...) SELECT_REP(__VA_ARGS__,REP3,REP2)(__VA_ARGS__)
#define per(i,r,l) for(int i=int(r)-1;i>=int(l);--i)
#define all(a) (a).begin(),(a).end()
#define eb emplace_back
#define si(a) int((a).size())
#define rep2(i,l,r) for(int i=int(l);i<=int(r);++i)
#define fi first
#define se second
int rnd(int n){static mt19937 r(1727);return r()%n;}
#include "graph/template.cpp"
template<class T>using Wgraph=WeightedGraph<T>;
#include "graph/SectionGraph.cpp"
#include "graph/dijkstra.cpp"
#include "graph/get_cycle.cpp"
#include "graph/MinCostFlow.cpp"
#include "data_structure/doubling.cpp"
#include "data_structure/2D_SparseTable.cpp"
#include "data_structure/dynamic_li_chao_tree.cpp"
#include "data_structure/treap.cpp"
#include "others/fastset.cpp"
#include "others/adjacent_manager.cpp"
#include "others/interval_manager.cpp"
#include "others/nd_for.cpp"
#include "math/kth_root.cpp"
#include "math/extgcd.cpp"
#include "math/prime_sieve.cpp"
#include "math/fft_aribitary_mod.cpp"
#include "linear_algebra/gauss_jordan.cpp"
#include "tree/Rerooting.cpp"
namespace f2 {
#include "linear_algebra/F2Matrix.cpp"
}
namespace fixed_matrix {
#include "math/matrix_square.cpp"
}
namespace random_impl {
#include "others/random.cpp"
}
namespace color_undo {
#include "data_structure/UnionFind_undo_coloring.cpp"
}
void check(bool b){if(!b)throw runtime_error("additional regression test");}
int main(){
    mt19937 rng(893471);
    for(int t=0;t<1000;++t){
        uint64_t a=(uint64_t(rng())<<32)|rng(),k=1+rng()%70,r=floor_root(a,k);
        auto fits=[&](uint64_t x){__uint128_t p=1;for(uint64_t j=0;j<k;++j){if(x&&p>a/x)return false;p*=x;}return p<=a;};
        check(fits(r));if(r!=UINT64_MAX)check(!fits(r+1));
    }
    check(floor_root(UINT64_MAX,2)==4294967295ULL);
    for(ll a=-5;a<=5;++a)for(ll b=-5;b<=5;++b)if(a||b)for(ll v=-20;v<=20;++v){
        EXTgcd<ll> eg(a,b);ll cnt=0,minx=100;
        for(ll x=-5;x<=5;++x)for(ll y=-5;y<=5;++y)if(a*x+b*y==v){++cnt;minx=min(minx,x);}
        check(eg.count(v,-5,-5,5,5)==cnt);auto [ok,xy]=eg.get(v,-5,-5,5,5);
        check(ok==bool(cnt));if(ok)check(xy.first==minx&&a*xy.first+b*xy.second==v&&xy.second>=-5&&xy.second<=5);
    }
    for(int n=0;n<1000;++n){auto primes=prime_sieve(n);vector<int> ref;for(int i=2;i<=n;++i){bool p=true;for(int j=2;j*j<=i;++j)if(i%j==0)p=false;if(p)ref.push_back(i);}check(primes==ref);}
    FastSet<2> fs(1000);set<int> ref;
    for(int t=0;t<10000;++t){int x=rng()%1000;if(rng()%2){fs.insert(x);ref.insert(x);}else{fs.erase(x);ref.erase(x);}
        int q=int(rng()%1200)-100;auto it=ref.lower_bound(q);check(fs.next(q)==(it==ref.end()?1000:*it));it=ref.upper_bound(q);check(fs.prev(q)==(it==ref.begin()?-1:*prev(it)));check(fs.empty()==ref.empty());}
    FastSet<1> no_elements(0);check(no_elements.empty()&&no_elements.next(0)==0&&no_elements.prev(0)==-1);
    auto adj=get_adj_manager<int,ll>([](int a,int b){return ll(b-a);},0LL);adj.add(2);adj.add(5);adj.add(5);check(adj()==3);adj.erase(99);check(adj()==3);
    IntervalManager<int> intervals(0);vector<int> values(30);
    for(int t=0;t<2000;++t){int l=rng()%31,r=rng()%31,x=rng()%4;if(l>r)swap(l,r);intervals.update(l,r,x);fill(values.begin()+l,values.begin()+r,x);const auto& c=intervals;rep(i,30)check(c[i]==values[i]);}
    IntervalManager<string> strings("");strings.update(1,4,"a");check(strings[2]=="a"&&strings[0]=="");
    int total=0;for(auto ij:ndFOR(2,3)){check(ij[0]<2&&ij[1]<3);++total;}check(total==6);
    for(auto ij:ndFOR(2,0)){(void)ij;check(false);}
    for(int t=0;t<300;++t){
        int h=1+rng()%5,w=1+rng()%5;BitMatrix a(h,w);vector<int>b(h),solution;rep(i,h){b[i]=rng()%2;rep(j,w)a[i][j]=rng()%2;}
        bool possible=false;rep(mask,1<<w){bool ok=true;rep(i,h){int v=0;rep(j,w)if(mask>>j&1)v^=a[i][j];if(v!=b[i])ok=false;}possible|=ok;}
        int rank=linear_equation(a,b,solution);check((rank>=0)==possible);if(possible)rep(i,h){int v=0;rep(j,w)v^=int(a[i][j])&solution[j];check(v==b[i]);}
        f2::Matrix<8,8> m=f2::Matrix<8,8>::I(5);rep(i,30){int a=rng()%5,b=rng()%5;if(a!=b)m[a]^=m[b];}auto inv=m.inv();
        rep(i,5)rep(j,5){int v=0;rep(k,5)v^=m[i][k]&inv[k][j];check(v==(i==j));}
        vector<vector<int>> grid(5,vector<int>(7));for(auto& row:grid)for(int& v:row)v=rng()%100;
        SparseTable2D<int,4> st(grid);rep(k,50){int l=rng()%5,r=l+1+rng()%(5-l),a=rng()%7,b=a+1+rng()%(7-a),mi=INT_MAX;rep(i,l,r)rep(j,a,b)mi=min(mi,grid[i][j]);check(st.query(l,a,r,b)==mi);}
    }
    fixed_matrix::Matrix<ll,3> matrix(2);matrix[0][0]=1;matrix[0][1]=2;matrix[1][0]=3;matrix[1][1]=4;
    check((matrix*matrix)[0][0]==7);check(vector<ll>({5,6})*matrix==vector<ll>({23,34}));
    random_impl::Random random(1);vector<int> empty;random.shuffle(empty);check(random.prob(1)&&!random.prob(0));set<vector<int>> permutations;
    rep(t,100){vector<int> v{0,1,2};random.shuffle(v);permutations.insert(v);}check(permutations.size()==6);check(random_impl::hash64(LLONG_MIN)!=random_impl::hash64(1));
    Doubling<ll> doubling(vector<int>{1,2,0},vector<ll>{2,3,4},10);check(doubling.move(0,5)==make_pair(2,14LL));Doubling<int> identity(3);check(identity.move(1,2)==make_pair(1,0));
    color_undo::UnionFind uf(3);check(uf.unite(0,1));check(uf.unite(1,2));check(!uf.unite(2,0));uf.rollback(3);check(uf.num==3);
    SectionGraph<ll> sg(6);Wgraph<ll> brute(6);
    rep(t,30){int l=rng()%6,r=l+1+rng()%(6-l),a=rng()%6,b=a+1+rng()%(6-a);ll cost=rng()%10;sg.add_edge(l,r,a,b,cost);rep(i,l,r)rep(j,a,b)brute[i].emplace_back(j,cost);}
    rep(s,6){auto a=dijkstra(sg.g,sg.sz+s,1LL<<60),b=dijkstra(brute,s,1LL<<60);rep(i,6)check(a[sg.sz+i]==b[i]);}
    DynamicLiChaoTree<ll,-50,50,(1LL<<60)> lc;vector<tuple<int,int,ll,ll>> lines;
    rep(t,200){ll a=int(rng()%31)-15,b=int(rng()%101)-50;int l=int(rng()%101)-50,r=l+1+rng()%(51-l);lc.add_segment(l,r,a,b);lines.emplace_back(l,r,a,b);for(ll x=-50;x<=50;++x){ll ref=1LL<<60;for(auto [l,r,a,b]:lines)if(l<=x&&x<r)ref=min(ref,a*x+b);check(lc.query(x)==ref);}}
    lc.clear();check(lc.query(0)==(1LL<<60));
    MinCostFlow<ll,ll> flow(4);flow.add_edge(0,1,2,-3);flow.add_edge(1,2,2,5);check(flow.flow(0,2,2)==make_pair(2LL,4LL));
    Graph tree{{1},{0,2},{1}};ReRooting<int,int> reroot(tree,[](int x,auto){return x+1;},[](int a,int b){return max(a,b);},0);check(reroot.build()==vector<int>({2,1,2}));check(reroot.build()==vector<int>({2,1,2}));
    Graph self{{1},{1}};DirectedCycle dc(self);check(dc.build()==vector<int>{1});check(dc.build()==vector<int>{1});
    Graph parallel{{1,1},{0,0}};UndirectedCycle uc(parallel);check(uc.build().size()==2);
    Treap<int,int> tr(20,[](int a,int b){return a+b;},[](int a,int b){return a+b;},-1,0);for(int v:{1,3,5})tr.insert(v);check(tr.erase(2,5));check(!tr.count(3)&&tr.count(1)&&tr.count(5));
    ArbitraryModConvolution<ll> conv;
    rep(t,200){int n=1+rng()%20,m=1+rng()%20;vector<ll>a(n),b(m);for(auto& x:a)x=int(rng()%200)-100;for(auto& x:b)x=int(rng()%200)-100;auto c=conv.multiply(a,b,1000000007);vector<ll> ref(n+m-1);rep(i,n)rep(j,m)ref[i+j]=(ref[i+j]+a[i]*b[j])%1000000007;for(auto& x:ref)if(x<0)x+=1000000007;check(c==ref);auto shortc=conv.multiply(a,b,1000000007,1);check(shortc[0]==ref[0]);}
    cout<<"additional: roots/Diophantine/GF2/intervals/FastSet/LiChao/FFT/flow/tree regressions passed\n";
}
