#include "geometry/halfplane.hpp"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using namespace cp_geometry;
void check(bool result, int line) {
    if (!result) throw std::runtime_error("geometry independent review at line " + std::to_string(line));
}
#define CHECK_REVIEW(condition) check(bool(condition), __LINE__)
bool close(Real a, Real b) { return std::abs(a-b) < 1e-7L * (1 + std::abs(a) + std::abs(b)); }
std::vector<Point> intersections(const std::vector<Point>& polygon, const Line& line) {
    std::vector<Point> result;
    auto add = [&](Point point) {
        for (Point previous : result) if (near(previous, point)) return;
        result.push_back(point);
    };
    Point direction = line.direction();
    for (std::size_t i = 0; i < polygon.size(); ++i) {
        Point a = polygon[i], b = polygon[(i+1)%polygon.size()];
        Real x = cross(direction,a-line.a), y = cross(direction,b-line.a);
        if (sign(x) == 0) add(a);
        if (sign(x)*sign(y) < 0) add(a+(b-a)*(x/(x-y)));
    }
    std::sort(result.begin(),result.end(),[&](Point a,Point b){return dot(direction,a)<dot(direction,b);});
    return result;
}
}

int main() {
    using namespace cp_geometry;
    // Rectangle edges and caller-supplied lines must use the same angle scale.
    // This used to miss the parallel group, then throw while intersecting it.
    for (Real slope : {0.L, 5e-13L, -5e-13L}) {
        std::vector<HalfPlane> constraints{{slope,-1,0}};
        auto polygon = half_plane_intersection(constraints,{-1000,-1000},{1000,1000});
        CHECK_REVIEW(close(polygon_area(polygon), 2000000));
    }
    for (Real intercept : {-1.L,1.L}) {
        std::vector<HalfPlane> constraints{{-1,0,0},{1,0,1000},{0,-1,0},{0,1,1000},
                                           {-5e-13L,-1,-intercept}};
        auto general = half_plane_intersection(constraints);
        CHECK_REVIEW(general.kind == RegionKind::POLYGON);
        CHECK_REVIEW(close(polygon_area(general.polygon),intercept<0 ? 1000000 : 999000));
        auto bounded = half_plane_intersection(constraints,{-1000,-1000},{1000,1000});
        CHECK_REVIEW(close(polygon_area(bounded),polygon_area(general.polygon)));
    }
    // Even extreme but finite coefficient scaling represents the same triangle.
    for (Real scale : {1e-3000L, 1.L, 1e3000L}) {
        std::vector<HalfPlane> constraints{{-scale,0,-scale},{0,-scale,2*scale},{scale,scale,5*scale}};
        auto region = half_plane_intersection(constraints);
        CHECK_REVIEW(region.kind == RegionKind::POLYGON);
        CHECK_REVIEW(close(polygon_area(region.polygon),18));
        CHECK_REVIEW(region.witness.x >= 1-1e-10L && region.witness.y >= -2-1e-10L);
        CHECK_REVIEW(region.witness.x+region.witness.y <= 5+1e-10L);
    }
    std::mt19937 random(4927354);
    {
        Real offset = 1e12L;
        std::vector<Point> clockwise{{offset,offset},{offset,offset+1},
                                    {offset+1,offset+1},{offset+1,offset}};
        CHECK_REVIEW(close(area2(clockwise),-2));
        auto polygon = ConvexPolygon::from_convex_hull(clockwise);
        CHECK_REVIEW(polygon.contains({offset+0.5L,offset+0.5L}) == Location::inside);
        CHECK_REVIEW(close(polygon_area(polygon.vertices()),1));
    }
    auto point = [&]{return Point{Real(int(random()%101)-50),Real(int(random()%101)-50)};};
    for (int trial = 0; trial < 400; ++trial) {
        std::vector<Point> cloud;
        for (int i = 0; i < int(random()%45)+3; ++i) cloud.push_back(point());
        auto hull = convex_hull(cloud);
        if (trial%2) std::reverse(hull.begin(),hull.end());
        std::rotate(hull.begin(),hull.begin()+random()%hull.size(),hull.end());
        auto polygon = ConvexPolygon::from_convex_hull(hull);
        for (int query = 0; query < 40; ++query) {
            Point a = point(), b = point();
            CHECK_REVIEW(polygon.contains(a) == contains(hull,a));
            int at = polygon.support(a);
            for (Point vertex : hull) CHECK_REVIEW(dot(a,polygon.vertices()[at])+1e-9L >= dot(a,vertex));
            if (near(a,b)) continue;
            Line line(a,b);
            auto expected = intersections(hull,line), actual = polygon.line_intersections(line);
            CHECK_REVIEW(expected.size() == actual.size());
            for (std::size_t i = 0; i < actual.size(); ++i) CHECK_REVIEW(near(actual[i],expected[i]));
        }
    }
    std::cout << "geometry review: near-parallel rectangle clipping, extreme coefficient scaling, 16000 ordered-hull queries passed\n";
}
