#include "geometry/halfplane.hpp"
#include <cassert>
#include <iostream>
#include <random>
using namespace cp_geometry;
bool close(Real a,Real b) { return std::abs(a-b)<1e-7L*(1+std::abs(b)); }
int main() {
    std::vector<HalfPlane> box{{-1,0,0},{0,-1,0},{1,0,4},{0,1,3}};
    auto region=half_plane_intersection(box);
    assert(region.kind==RegionKind::POLYGON && close(polygon_area(region.polygon),12));
    assert(half_plane_intersection(std::vector<HalfPlane>{}).kind==RegionKind::UNBOUNDED);
    assert(half_plane_intersection(std::vector<HalfPlane>{{0,0,-1}}).kind==RegionKind::EMPTY);
    assert(half_plane_intersection(std::vector<HalfPlane>{{0,0,0},{0,0,1}}).kind==RegionKind::UNBOUNDED);
    assert(half_plane_intersection(std::vector<HalfPlane>{{1,0,0},{-1,0,-1}}).kind==RegionKind::EMPTY);
    assert(half_plane_intersection(std::vector<HalfPlane>{{1,0,0},{-1,0,0}}).kind==RegionKind::UNBOUNDED);
    auto point=half_plane_intersection(std::vector<HalfPlane>{{1,0,0},{-1,0,0},{0,1,0},{0,-1,0}});
    assert(point.kind==RegionKind::DEGENERATE && near(point.witness,{0,0}));
    // Point segments must use the same distance tolerance as a point hull.
    assert(!on_segment({0,0},{0,0},{1e-8L,0}));
    auto hit=line_intersection(Line({0,0},{1e-8L,0}),Line({1,1},{1,1+1e-8L}));
    assert(near(hit,{1,0}));
    std::mt19937 rng(1925367);
    for(int trial=0;trial<3000;++trial) {
        std::vector<HalfPlane> constraints=box;
        std::vector<Point> clipped{{0,0},{4,0},{4,3},{0,3}};
        for(int i=0;i<int(rng()%25);++i) {
            int a=int(rng()%21)-10,b=int(rng()%21)-10,c=int(rng()%61)-20;
            if(a==0&&b==0) continue;
            constraints.push_back({Real(a),Real(b),Real(c)});
            Point normal{Real(a),Real(b)};
            Point origin=normal*(Real(c)/norm2(normal));
            clipped=convex_cut(clipped,Line(origin,origin+perpendicular(normal)));
        }
        // Multiplying any inequality by a positive scalar must not alter it.
        for(auto& h:constraints) {
            Real scale=std::pow(10.L,int(rng()%17)-8);
            h.a*=scale; h.b*=scale; h.c*=scale;
        }
        auto answer=half_plane_intersection(constraints);
        assert(close(polygon_area(answer.polygon),polygon_area(clipped)));
        assert(close(polygon_area(half_plane_intersection(constraints,{-10,-10},{10,10})),polygon_area(clipped)));
        if(answer.kind!=RegionKind::EMPTY) for(const auto& h:constraints) {
            Real scale=std::hypot(h.a,h.b);
            assert((h.a*answer.witness.x+h.b*answer.witness.y-h.c)/scale<1e-7L);
            for(auto vertex:answer.polygon)
                assert((h.a*vertex.x+h.b*vertex.y-h.c)/scale<1e-7L);
        }
    }
    std::vector<Point> hull{{0,0},{5,0},{4,3},{1,4}};
    for(int trial=0;trial<8;++trial) {
        if(trial==4) std::reverse(hull.begin(),hull.end());
        std::rotate(hull.begin(),hull.begin()+1,hull.end());
        auto p=ConvexPolygon::from_convex_hull(hull);
        for(int x=-1;x<=6;++x) for(int y=-1;y<=5;++y)
            assert(p.contains({Real(x),Real(y)})==contains(hull,{Real(x),Real(y)}));
        for(int x=-8;x<=8;++x) for(int y=-8;y<=8;++y) {
            Point direction{Real(x),Real(y)};
            Real actual=dot(p.vertices()[p.support(direction)],direction);
            for(auto v:hull) assert(actual+1e-10L>=dot(v,direction));
        }
    }
    std::cout<<"halfplane: 3000 coefficient/scaling cases, ordered hulls and tiny lines passed\n";
}
