#include <bits/stdc++.h>
using namespace std;
using ll=long long;
#define rep(i,n) for(int i=0;i<int(n);++i)
#define si(x) int((x).size())
#include ".vscode/geometry/Points.cpp"
#include ".vscode/geometry/geometry.hpp"
int main() {
    cp_integer_geometry::Point<ll> zero{0,0};assert(!(zero<zero));
    vector<cp_integer_geometry::Point<ll>> v{{0,0},{1,0},{0,0},{0,1},{-1,0},{0,-1}};
    for(auto a:v)for(auto b:v) {assert(!(a<b&&b<a));for(auto c:v)if(a<b&&b<c)assert(a<c);}
    const Geometry::Point a{2,3};auto b=a*4LL;
    assert(a.x==2&&a.y==3&&b.x==8&&b.y==12);
    cout<<"legacy geometry regressions passed\n";
}
