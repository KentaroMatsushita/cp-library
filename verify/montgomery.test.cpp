#include "my_template.hpp"
#include "modular/mod1000000007.cpp"
#include "modular/mod998244353.cpp"
template<class FPS> void verify_fps(){
    using M=typename FPS::value_type;
    mt19937 gen(146178);
    for(int n:{0,1,7,140}){
        FPS a(n),b(n);for(auto& x:a)x=gen()%1000;for(auto& x:b)x=gen()%1000;
        auto c=a*b;FPS reference(n?2*n-1:0);rep(i,n)rep(j,n)reference[i+j]+=a[i]*b[j];assert(c==reference);
        if(n){FPS f{1,1};auto inverse=f.inv(n);rep(i,n)assert(inverse[i]==M(i%2?-1:1));}
    }
    FPS zero(4);auto one=zero.pow(0,7);assert(si(one)==7&&one[0]==M(1));rep(i,1,7)assert(one[i]==M(0));
}
int main(){
    // Large transforms in these bundled implementations require AVX2.
    if(!__builtin_cpu_supports("avx2"))throw runtime_error("AVX2 required for bundled NTT tests");
    verify_fps<Modular107::fps>();verify_fps<Modular998::fps>();
    assert(Modular107::C(10,3).get()==120);assert(Modular998::C(10,3).get()==120);
    cout<<"Montgomery NTT/FPS regressions passed for both fixed moduli\n";
}
