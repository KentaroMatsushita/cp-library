#include "geometry/delaunay.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace cp_geometry;
void require(bool condition,int line) {
    if(!condition) throw std::runtime_error("Delaunay verification failed at line "+std::to_string(line));
}
#define CHECK_DELAUNAY(condition) require(bool(condition),__LINE__)
using Wide=__int128_t;
Wide orientation(Point a,Point b,Point c) {
    Wide ax=static_cast<long long>(a.x)-static_cast<long long>(c.x);
    Wide ay=static_cast<long long>(a.y)-static_cast<long long>(c.y);
    Wide bx=static_cast<long long>(b.x)-static_cast<long long>(c.x);
    Wide by=static_cast<long long>(b.y)-static_cast<long long>(c.y);
    return ax*by-ay*bx;
}
Wide circle(Point a,Point b,Point c,Point d) {
    Wide ax=static_cast<long long>(a.x)-static_cast<long long>(d.x);
    Wide ay=static_cast<long long>(a.y)-static_cast<long long>(d.y);
    Wide bx=static_cast<long long>(b.x)-static_cast<long long>(d.x);
    Wide by=static_cast<long long>(b.y)-static_cast<long long>(d.y);
    Wide cx=static_cast<long long>(c.x)-static_cast<long long>(d.x);
    Wide cy=static_cast<long long>(c.y)-static_cast<long long>(d.y);
    return (ax*ax+ay*ay)*(bx*cy-by*cx)+(bx*bx+by*by)*(cx*ay-cy*ax)+
           (cx*cx+cy*cy)*(ax*by-ay*bx);
}
int wide_sign(Wide value) { return (value>0)-(value<0); }
void verify_integer(const std::vector<Point>& input) {
    auto result=delaunay_triangulation(input);
    CHECK_DELAUNAY(result.original_to_site.size()==input.size());
    CHECK_DELAUNAY(result.representative.size()==result.sites.size());
    for(std::size_t i=0;i<input.size();++i) {
        auto site=result.original_to_site[i];
        CHECK_DELAUNAY(site<result.sites.size());
        CHECK_DELAUNAY(result.sites[site].x==input[i].x && result.sites[site].y==input[i].y);
        CHECK_DELAUNAY(result.representative[site]<=i);
    }
    const auto count=result.sites.size();
    std::set<std::pair<std::size_t,std::size_t>> edges(result.edges.begin(),result.edges.end());
    CHECK_DELAUNAY(edges.size()==result.edges.size());
    std::vector<std::size_t> component(count);
    std::iota(component.begin(),component.end(),std::size_t(0));
    auto root=[&](auto&& self,std::size_t x)->std::size_t {
        return component[x]==x?x:component[x]=self(self,component[x]);
    };
    for(auto edge:result.edges) {
        CHECK_DELAUNAY(edge.first<edge.second && edge.second<count);
        component[root(root,edge.first)]=root(root,edge.second);
        for(std::size_t i=0;i<count;++i) if(i!=edge.first&&i!=edge.second) {
            Point a=result.sites[edge.first],b=result.sites[edge.second],p=result.sites[i];
            if(orientation(a,b,p)==0) CHECK_DELAUNAY(dot(p-a,p-b)>0);
        }
    }
    for(std::size_t i=1;i<count;++i) CHECK_DELAUNAY(root(root,i)==root(root,0));
    // No two edges cross in their interiors.
    for(std::size_t i=0;i<result.edges.size();++i) for(std::size_t j=0;j<i;++j) {
        auto e=result.edges[i],f=result.edges[j];
        if(e.first==f.first||e.first==f.second||e.second==f.first||e.second==f.second) continue;
        Point a=result.sites[e.first],b=result.sites[e.second],c=result.sites[f.first],d=result.sites[f.second];
        CHECK_DELAUNAY(!(wide_sign(orientation(a,b,c))*wide_sign(orientation(a,b,d))<0 &&
                        wide_sign(orientation(c,d,a))*wide_sign(orientation(c,d,b))<0));
    }
    bool noncollinear=false;
    for(std::size_t i=2;i<count;++i) if(orientation(result.sites[0],result.sites[1],result.sites[i])!=0)
        noncollinear=true;
    if(!noncollinear) {
        CHECK_DELAUNAY(result.triangles.empty());
        CHECK_DELAUNAY(result.edges.size()==(count?count-1:0));
        CHECK_DELAUNAY(result.hull.size()==count);
        return;
    }
    CHECK_DELAUNAY(result.triangles.size()==2*count-2-result.hull.size());
    CHECK_DELAUNAY(result.edges.size()==3*count-3-result.hull.size());
    std::set<std::array<std::size_t,3>> triangles;
    Wide triangle_area=0;
    for(auto triangle:result.triangles) {
        CHECK_DELAUNAY(triangle[0]<triangle[1] && triangle[0]<triangle[2]);
        CHECK_DELAUNAY(triangles.insert(triangle).second);
        Point a=result.sites[triangle[0]],b=result.sites[triangle[1]],c=result.sites[triangle[2]];
        Wide area=orientation(a,b,c);
        CHECK_DELAUNAY(area>0); triangle_area+=area;
        for(std::size_t i=0;i<count;++i) CHECK_DELAUNAY(circle(a,b,c,result.sites[i])<=0);
        for(int k=0;k<3;++k) {
            auto x=triangle[k],y=triangle[(k+1)%3];
            if(x>y) std::swap(x,y);
            CHECK_DELAUNAY(edges.count({x,y})==1);
        }
    }
    Wide hull_area=0;
    for(std::size_t i=1;i+1<result.hull.size();++i)
        hull_area+=orientation(result.sites[result.hull[0]],result.sites[result.hull[i]],result.sites[result.hull[i+1]]);
    CHECK_DELAUNAY(triangle_area==hull_area);
    // Construction depends only on sites, not duplicate ordering.
    auto again=delaunay_triangulation(result.sites);
    CHECK_DELAUNAY(again.edges==result.edges && again.triangles==result.triangles);
}
}

int main(int argc,char** argv) {
    if(argc>1 && std::string(argv[1])=="--benchmark") {
        std::mt19937_64 random(718234);
        std::vector<Point> points;
        points.reserve(100000);
        for(int i=0;i<100000;++i) points.push_back({Real(random()%100000000),Real(random()%100000000)});
        const auto started=std::chrono::steady_clock::now();
        auto triangulation=delaunay_triangulation(points);
        CHECK_DELAUNAY(triangulation.triangles.size()==2*triangulation.sites.size()-2-triangulation.hull.size());
        CHECK_DELAUNAY(triangulation.edges.size()==3*triangulation.sites.size()-3-triangulation.hull.size());
        std::cout<<"Delaunay benchmark: "<<triangulation.sites.size()<<" sites, "<<triangulation.triangles.size()
                 <<" triangles, "<<std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count()<<" seconds\n";
        return 0;
    }
    for(auto points:std::vector<std::vector<Point>>{{},{{1,2}},{{1,2},{1,2}},{{0,0},{3,4}},
        {{0,0},{1,0},{2,0},{3,0}},{{0,3},{0,0},{0,2},{0,1}},
        {{0,0},{2,0},{0,2}},{{0,0},{2,0},{2,2},{0,2}},
        {{0,0},{1,0},{2,0},{2,1},{2,2},{1,2},{0,2},{0,1},{1,1}}}) verify_integer(points);
    std::vector<Point> circle_points;
    for(int x=-5;x<=5;++x) for(int y=-5;y<=5;++y) if(x*x+y*y==25) circle_points.push_back({Real(x),Real(y)});
    verify_integer(circle_points);
    for(int width=2;width<=12;++width) {
        std::vector<Point> grid;
        for(int x=0;x<width;++x) for(int y=0;y<width;++y) grid.push_back({Real(x),Real(y)});
        verify_integer(grid);
    }
    std::mt19937_64 random(328754128);
    for(int trial=0;trial<1500;++trial) {
        std::vector<Point> points;
        for(int i=0,count=int(random()%35);i<count;++i)
            points.push_back({Real(int(random()%61)-30),Real(int(random()%61)-30)});
        verify_integer(points);
    }
    for(int trial=0;trial<100;++trial) {
        std::vector<Point> points;
        for(int i=0;i<15;++i)
            points.push_back({Real(int(random()%200000001)-100000000),Real(int(random()%200000001)-100000000)});
        verify_integer(points);
    }
    const Real tiny_delta=std::ldexp(Real(1),4-std::numeric_limits<Real>::digits);
    for(Real delta:{tiny_delta,-tiny_delta}) {
        auto result=delaunay_triangulation({{0,0},{1,0},{1,1},{0,1+delta}});
        auto diagonal=delta>0?std::make_pair(std::size_t(0),std::size_t(3)):
                              std::make_pair(std::size_t(1),std::size_t(2));
        CHECK_DELAUNAY(std::binary_search(result.edges.begin(),result.edges.end(),diagonal));
    }
    for(Real scale:{1e-2000L,1.L,1e2000L}) {
        auto result=delaunay_triangulation({{0,0},{scale,0},{scale,scale},{0,scale}});
        CHECK_DELAUNAY(result.triangles.size()==2 && result.edges.size()==5);
    }
    {
        // Consecutive Fibonacci vectors have determinant exactly 1, although
        // products around 1e31 cancel. Plain long double loses this orientation.
        const Point a{2111485077978050.L,3416454622906707.L};
        const Point b{3416454622906707.L,5527939700884757.L};
        CHECK_DELAUNAY(orientation({0,0},a,b)==1);
        auto result=delaunay_triangulation({{0,0},a,b,a+b});
        CHECK_DELAUNAY(result.triangles.size()==2 && result.edges.size()==5);
        CHECK_DELAUNAY(std::binary_search(result.edges.begin(),result.edges.end(),
                                         std::make_pair(std::size_t(1),std::size_t(2))));
    }
    {
        Real offset=1e15L;
        auto result=delaunay_triangulation({{offset,offset},{offset+1,offset},{offset,offset+1},{offset+1,offset+1}});
        CHECK_DELAUNAY(result.triangles.size()==2 && result.edges.size()==5);
    }
    bool invalid=false;
    try { (void)delaunay_triangulation({{0,0},{std::numeric_limits<Real>::infinity(),0}}); }
    catch(const std::invalid_argument&) { invalid=true; }
    CHECK_DELAUNAY(invalid);
    std::cout<<"Delaunay: 1600 integer empty-circle/topology differential cases, grids, cocircular/collinear/duplicate and near-cocircular cases passed\n";
}
