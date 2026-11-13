#include <bits/stdc++.h>
#include <ext/pb_ds/assoc_container.hpp>
using namespace std;
using ll=long long;
using vi=vector<int>;
using pii=pair<int,int>;
using Graph=vector<vector<int>>;
#define rep(i,n) for(int i=0;i<int(n);++i)
#define fi first
#define se second
#define si(a) int((a).size())
template<class T>bool chmin(T& a,T b){if(b<a){a=b;return true;}return false;}
const pii dx4[]={{1,0},{0,1},{-1,0},{0,-1}};
bool inc(int x,int l,int r){return l<=x&&x<r;}
#include "others/Grid.cpp"
#include "others/2D_calculations.cpp"
#include "data_structure/BFS_restore.cpp"
#include "tree/centroid.cpp"
#include "tree/EulerTour.cpp"
#include "data_structure/hash_map_gnu.cpp"
#include "math/binom.cpp"
#include "math/fft.cpp"
#include "graph/TSP.cpp"
#include "graph/template.cpp"
template<class T>using Wgraph=WeightedGraph<T>;
#include "math/mobius_function_table.cpp"
#include "math/fast_gcd.cpp"
#include "data_structure/BIT01.cpp"
#include "data_structure/inversion_number.cpp"
namespace conv_and {
#include "math/and_convolution.cpp"
}
namespace conv_or {
#include "math/or_convolution.cpp"
}
namespace conv_xor {
#include "math/xor_convolution.cpp"
}
namespace edge_tour {
#include "tree/EulerTour_Edgebase.cpp"
}
namespace weighted_diameter {
#include "tree/diameter_weighted.cpp"
}
void check(bool b){if(!b)throw runtime_error("short source regression");}
int main(){
    auto sum=pii{1,2}+pii{3,4};check(sum==pii(4,6));
    Grid grid(3,4);for(int i=0;i<3;++i)for(int j=0;j<4;++j){int count=0;for(auto [x,y]:grid(i,j)){check(inc(x,0,3)&&inc(y,0,4));++count;}check(count==4-(i==0)-(i==2)-(j==0)-(j==3));}
    Grid tiny(1,1);for(auto p:tiny(0,0)){(void)p;check(false);}
    Graph g{{1},{0,2},{1}};
    check(BFS(g,2,0)==vector<int>({2,1,0}));check(BFS(g,1,1)==vector<int>{1});
    EulerTour tour(g);auto [in,out]=tour.build(2);check(in[2]==0&&out[2]==3&&out[0]-in[0]==1);
    struct E {int to,id;operator int()const{return to;}};
    vector<vector<E>> edges{{{1,0}},{{0,0},{2,1}},{{1,1}}};
    edge_tour::EulerTour et(edges);auto [ei,eo]=et.build(2);check(ei[1]==0&&eo[1]==2&&eo[0]-ei[0]==1);
    Wgraph<int> wg(3);wg[0].emplace_back(1,4);wg[1].emplace_back(0,4);wg[1].emplace_back(2,7);wg[2].emplace_back(1,7);
    weighted_diameter::Diameter<int> wd(wg);check(wd.build()==11);check(wd.build()==11&&wd.path.size()==3);
    Graph empty;check(centroid(empty).empty());check(centroid(g)==vector<int>{1});
    hash_map<int,int> h;h[42]=3;check(h[42]==3);
    check(C(5,2)==10&&C(5,8)==0&&C(-1,0)==0);
    bool overflow=false;try{C(67,33);}catch(const overflow_error&){overflow=true;}check(overflow);
    mt19937 rng(516238);
    FFT fft;
    for(int n:{1,2,4,8,16,32}){
        vector<complex<double>> a(n);for(auto& v:a)v={double(rng()%21)-10,double(rng()%21)-10};auto original=a;
        fft.fft(a);for(int k=0;k<n;++k){complex<double> ref=0;for(int j=0;j<n;++j)ref+=original[j]*polar(1.0,-2*acos(-1.0)*j*k/n);check(abs(a[k]-ref)<1e-8);}
        fft.fft(a,true);rep(i,n)check(abs(a[i]-original[i])<1e-8);
    }
    vector<complex<double>> bad(3);bool rejected=false;try{fft.fft(bad);}catch(const invalid_argument&){rejected=true;}check(rejected);
    vector<vector<long long>> d{{0,3000000000LL},{3000000000LL,0}};check(TSP(d)==6000000000LL);check(TSP(vector<vector<int>>{})==0);
    for(int t=0;t<1000;++t){
        int n=1+rng()%64;BIT01 bits(n);vector<int>ref(n);vector<int>a(n);
        rep(k,200){int x=rng()%n,y=rng()%2;bits.set(x,y);ref[x]=y;int l=rng()%(n+1),r=rng()%(n+1);if(l>r)swap(l,r);check(bits.query(l,r)==accumulate(ref.begin()+l,ref.begin()+r,0));}
        for(int& v:a)v=int(rng()%50)-25;long long inv=0;rep(i,n)rep(j,i)inv+=a[j]>a[i];check(inversion(a,false)==inv);
        check(fgcd(rng()%100000,rng()%100000)>=0);
        int m=1<<(rng()%5);vector<long long>x(m),y(m);for(auto& v:x)v=rng()%20;for(auto& v:y)v=rng()%20;
        vector<long long> ca(m),co(m),cx(m);rep(i,m)rep(j,m){ca[i&j]+=x[i]*y[j];co[i|j]+=x[i]*y[j];cx[i^j]+=x[i]*y[j];}
        check(conv_and::and_conv(x,y)==ca);check(conv_or::or_conv(x,y)==co);check(conv_xor::xor_conv(x,y)==cx);
    }
    auto mu=mobius(100);for(int n=1;n<=100;++n){int sum=0;for(int k=1;k<=n;++k)if(n%k==0)sum+=mu[k];check(sum==(n==1));}
    cout<<"short sources: FFT/bitwise convolution/BIT/inversion differential tests and regressions passed\n";
}
