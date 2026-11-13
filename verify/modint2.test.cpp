#include "my_template.hpp"
#include "modular/modint2.cpp"
int main(){
    mt19937 gen(481612);
    for(int n:{0,1,3,31,64,128}) {
        Poly a(n),b(n);for(auto& x:a)x=gen()%MOD;for(auto& x:b)x=gen()%MOD;
        auto c=a*b;Poly reference(n?2*n-1:0);rep(i,n)rep(j,n)reference[i+j]+=a[i]*b[j];assert(c==reference);
    }
    for(int n:{0,1,2,7,16,65}) {
        Poly f{0,1};auto e=exp(f,n);assert(si(e)==n);rep(i,n)assert(e[i]==ifact(i));
        auto zero=exp(Poly{},n);assert(si(zero)==n);rep(i,n)assert(zero[i]==mint(i==0));
        if(n){auto l=log(e,n);rep(i,n)assert(l[i]==mint(i==1));}
        auto one=pow(Poly(4),0,n);assert(si(one)==n);rep(i,n)assert(one[i]==mint(i==0));
        auto invf=inv(Poly{1,-1},n);rep(i,n)assert(invf[i]==mint(1));
    }
    cout<<"modint2: NTT and FPS differential identities passed\n";
}
