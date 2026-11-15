#include <bits/stdc++.h>
#include "geometry/geometry.hpp"
using namespace std;
using namespace cp_geometry;
#define check(b) do {if(!(b)) throw runtime_error(string("geometry failure at line ")+to_string(__LINE__));}while(0)
bool close(Real a,Real b) {return abs(a-b)<1e-7L;}
vector<Point> brute_intersections(const vector<Point>& p,const Line& l) {
    vector<Point> q;
    auto add=[&](Point x){for(Point y:q) if(near(x,y)) return; q.push_back(x);};
    for(size_t i=0;i<p.size();++i) {
        Point a=p[i],b=p[(i+1)%p.size()];
        Real x=cross(l.direction(),a-l.a),y=cross(l.direction(),b-l.a);
        if(sign(x)==0) add(a);
        if(sign(x)*sign(y)<0) add(a+(b-a)*(x/(x-y)));
    }
    sort(q.begin(),q.end(),[&](Point a,Point b){return dot(l.direction(),a)<dot(l.direction(),b);});
    return q;
}
int main() {
    constexpr Point origin;
    constexpr Point constructed(2,3);
    static_assert(origin.x==0 && origin.y==0);
    static_assert(constructed.x==2 && constructed.y==3);
    std::vector<Point> constructed_points;
    constructed_points.emplace_back(4,5);
    check(constructed_points[0].x==4 && constructed_points[0].y==5);
    mt19937 rng(928431);
    auto random_point=[&](){return Point{Real(int(rng()%201)-100),Real(int(rng()%201)-100)};};
    for(int t=0;t<1000;++t) {
        vector<Point> cloud;
        for(int k=0;k<int(rng()%70);++k) cloud.push_back(random_point());
        ConvexPolygon poly(cloud); const auto& p=poly.vertices();
        Real diameter=0;
        for(Point a:p) for(Point b:p) diameter=max(diameter,norm2(a-b));
        check(close(diameter,convex_diameter2(p)));
        for(int k=0;k<100;++k) {
            Point a=random_point(),b=random_point();
            check(poly.contains(a)==cp_geometry::contains(p,a));
            if(!p.empty()) {
                Real best=dot(a,p[poly.support(a)]);
                for(Point v:p) check(best+1e-7L>=dot(a,v));
            }
            if(near(a,b)) continue;
            Line l(a,b);
            auto x=poly.line_intersections(l),y=brute_intersections(p,l);
            check(x.size()==y.size());
            for(size_t i=0;i<x.size();++i) check(norm2(x[i]-y[i])<1e-12L);
        }
        // All vertices, edge-coincidence and tangencies.
        for(int i=0;i<int(p.size());++i) {
            check(poly.contains(p[i])==Location::boundary);
            if(p.size()>1) {
                Line l(p[i],p[(i+1)%p.size()]);
                auto x=poly.line_intersections(l),y=brute_intersections(p,l);
                check(x.size()==y.size());
                for(size_t j=0;j<x.size();++j) check(norm2(x[j]-y[j])<1e-12L);
            }
        }
    }
    Point low{-100,-100},high{100,100};
    vector<Point> box{low,{100,-100},high,{-100,100}};
    for(int t=0;t<2000;++t) {
        vector<Line> lines;
        vector<Point> reference=box;
        for(int k=0;k<int(rng()%30);++k) {
            Point a=random_point(),b=random_point(); if(near(a,b)) continue;
            lines.emplace_back(a,b);
            reference=convex_cut(reference,lines.back());
        }
        auto result=half_plane_intersection(lines,low,high);
        check(close(area2(result),area2(reference)));
        auto bounded=lines;
        for(int i=0;i<4;++i) bounded.emplace_back(box[i],box[(i+1)%4]);
        auto general=half_plane_intersection(bounded);
        check(close(area2(general.polygon),area2(reference)));
        if(general.kind!=RegionKind::EMPTY) for(const auto& l:bounded)
            check(cross(l.direction(),general.witness-l.a)>-1e-7L);
        for(Point p:result) for(const auto& l:lines) check(cross(l.direction(),p-l.a)>-1e-7L);
    }
    // Parallel/repeated/opposed constraints, empty and lower-dimensional intersections.
    for(auto lines:vector<vector<Line>>{
        {},{Line({0,0},{1,0}),Line({0,1},{2,1})},
        {Line({0,0},{1,0}),Line({1,0},{0,0})},
        {Line({0,1},{1,1}),Line({1,0},{0,0})},
        {Line({0,0},{1,0}),Line({0,0},{1,0})}}) {
        auto ref=box; for(const auto& l:lines) ref=convex_cut(ref,l);
        check(close(area2(half_plane_intersection(lines,low,high)),area2(ref)));
    }
    check(segments_intersect({0,0},{0,0},{-1,0},{1,0}));
    check(!segments_intersect({0,0},{0,0},{1,0},{2,0}));
    check(close(segment_distance({0,0},{0,0},{3,4}),5));
    check(near(projection(Line({0,0},{2,0}),{1,3}),{1,0}));
    check(half_plane_intersection(std::vector<Line>{}).kind==RegionKind::UNBOUNDED);
    check(half_plane_intersection(vector<Line>{Line({0,0},{1,0})}).kind==RegionKind::UNBOUNDED);
    check(half_plane_intersection(vector<Line>{Line({0,1},{1,1}),Line({1,0},{0,0})}).kind==RegionKind::EMPTY);
    check(half_plane_intersection(vector<Line>{Line({0,0},{1,0}),Line({1,0},{0,0})}).kind==RegionKind::UNBOUNDED);
    vector<Line> point{Line({0,0},{1,0}),Line({0,0},{0,1}),Line({0,0},{-1,0}),Line({0,0},{0,-1})};
    check(half_plane_intersection(point).kind==RegionKind::DEGENERATE);
    cout<<"geometry: 1000 hulls, 100000 queries, 2000 half-plane differential cases passed\n";
}
