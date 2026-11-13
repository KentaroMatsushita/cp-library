#include "my_template.hpp"
#include "modular/modint1.cpp"
int main(){
    modular::modint<17> custom(16);++custom;assert(custom.a==0);custom--;assert(custom.a==16);
    assert((modular::modint<17>(5)/modular::modint<17>(3)).a==13);
    mt19937 gen(1287);
    for(int n:{0,1,3,99,100,160}) {
        Poly a(n),b(n);for(auto& x:a)x=gen()%MOD;for(auto& x:b)x=gen()%MOD;
        auto c=a*b;Poly reference(n?2*n-1:0);
        rep(i,n)rep(j,n)reference[i+j]+=a[i]*b[j];assert(c==reference);
    }
    assert(C(5,2)==mint(10));
    cout<<"modint1: nondefault modulus and wide CRT convolution passed\n";
}
