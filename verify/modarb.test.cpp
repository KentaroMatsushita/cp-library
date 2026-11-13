#include "my_template.hpp"
#include "modular/mod_arbitary.cpp"
int main(){
    mint::set_mod(1000000007);
    Poly a{1,2,3},b{3,4};auto c=a*b;assert(c==Poly({3,10,17,12}));
    auto inverse=Poly{1,1}.inv(20);rep(i,20)assert(inverse[i]==mint(i%2?-1:1));
    auto one=Poly{}.pow(0,10);assert(one[0]==mint(1));
    cout<<"runtime modulus FPS regressions passed\n";
}
