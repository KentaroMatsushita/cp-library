#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>
#include "geometry/algorithms.hpp"

using cp_geometry::Point;
using cp_geometry::Real;

void require(bool condition, int line) {
    if (!condition) throw std::runtime_error("geometry algorithms failure at line " + std::to_string(line));
}
#define CHECK(condition) require(bool(condition), __LINE__)

bool close(Real a, Real b, Real tolerance = 1e-9L) {
    return std::abs(a-b) <= tolerance * (1 + std::abs(a) + std::abs(b));
}
bool equal_hulls(const std::vector<Point>& a, const std::vector<Point>& b) {
    if (a.size() != b.size()) return false;
    for (Point point : a) {
        bool found = false;
        for (Point other : b) if (cp_geometry::near(point, other)) found = true;
        if (!found) return false;
    }
    return true;
}

int main() {
    std::mt19937_64 random(93876153);
    auto random_point = [&] {
        return Point{Real(int(random()%81)-40), Real(int(random()%81)-40)};
    };
    CHECK(!cp_geometry::closest_pair({}));
    CHECK(!cp_geometry::closest_pair({{3,4}}));
    bool invalid_point_threw = false;
    try { (void)cp_geometry::closest_pair({{std::numeric_limits<Real>::infinity(),0}}); }
    catch (const std::invalid_argument&) { invalid_point_threw = true; }
    CHECK(invalid_point_threw);
    for (int iteration = 0; iteration < 3000; ++iteration) {
        std::vector<Point> points;
        for (int i = 0, count = int(random()%100)+2; i < count; ++i) points.push_back(random_point());
        auto result = cp_geometry::closest_pair(points);
        Real brute = std::numeric_limits<Real>::infinity();
        for (std::size_t i = 0; i < points.size(); ++i)
            for (std::size_t j = 0; j < i; ++j)
                brute = (std::min)(brute, cp_geometry::norm2(points[i] - points[j]));
        CHECK(result && result->first < points.size() && result->second < points.size());
        CHECK(result->first != result->second);
        CHECK(close(result->squared_distance, brute));
        CHECK(close(cp_geometry::norm2(points[result->first] - points[result->second]), brute));
    }
    {
        // Same x coordinates, dense strips, tiny nonzero separation and large translation.
        std::vector<Point> points;
        for (int i = 0; i < 50000; ++i) points.push_back({Real(i%2)*1e-8L, Real(i)});
        points.push_back({0, 1e-8L});
        CHECK(close(cp_geometry::closest_pair(points)->squared_distance, 1e-16L, 1e-24L));
        points = {{-1e9L,1e9L},{1e9L,-1e9L}};
        CHECK(close(cp_geometry::closest_pair(points)->squared_distance, 8e18L));
        points = {{1,1},{2,2},{1,1}};
        CHECK(cp_geometry::closest_pair(points)->squared_distance == 0);
    }
    for (int iteration = 0; iteration < 3000; ++iteration) {
        std::vector<Point> a, b;
        for (int i = 0, count = int(random()%20); i < count; ++i) a.push_back(random_point());
        for (int i = 0, count = int(random()%20); i < count; ++i) b.push_back(random_point());
        auto hull_a = cp_geometry::convex_hull(a), hull_b = cp_geometry::convex_hull(b);
        if (!hull_a.empty()) std::rotate(hull_a.begin(), hull_a.begin() + random()%hull_a.size(), hull_a.end());
        if (!hull_b.empty()) std::rotate(hull_b.begin(), hull_b.begin() + random()%hull_b.size(), hull_b.end());
        std::vector<Point> sums;
        for (Point x : a) for (Point y : b) sums.push_back(x+y);
        auto expected = cp_geometry::convex_hull(sums);
        auto actual = cp_geometry::minkowski_sum_convex(hull_a, hull_b);
        CHECK(equal_hulls(actual, expected));
        CHECK(equal_hulls(cp_geometry::minkowski_sum(a,b), expected));
        CHECK(close(cp_geometry::area2(actual), cp_geometry::area2(expected)));
        if (actual.size() > 2) for (std::size_t i = 0; i < actual.size(); ++i)
            CHECK(cp_geometry::sign(cp_geometry::cross(actual[(i+1)%actual.size()]-actual[i],
                  actual[(i+2)%actual.size()]-actual[(i+1)%actual.size()])) > 0);
    }
    for (auto a : std::vector<std::vector<Point>>{{},{{1,2}},{{-3,0},{7,0}},{{2,-3},{2,5}},
                  {{-2,0},{0,-2},{2,0},{0,2}}})
        for (auto b : std::vector<std::vector<Point>>{{},{{0,0}},{{-4,0},{8,0}},{{1,2},{2,1}},
                      {{-1,-1},{1,-1},{1,1},{-1,1}}}) {
            a = cp_geometry::convex_hull(a); b = cp_geometry::convex_hull(b);
            std::vector<Point> sums;
            for (Point x : a) for (Point y : b) sums.push_back(x+y);
            CHECK(equal_hulls(cp_geometry::minkowski_sum_convex(a,b), cp_geometry::convex_hull(sums)));
        }
    for (int iteration = 0; iteration < 1000; ++iteration) {
        Point a = random_point(), b = random_point(), c = random_point();
        if (cp_geometry::sign(cp_geometry::cross(b-a,c-a)) == 0) continue;
        Point center = cp_geometry::polygon_centroid({a,b,c});
        CHECK(cp_geometry::near(center, (a+b+c)/3));
        CHECK(cp_geometry::near(cp_geometry::polygon_centroid({c,b,a}), center));
        Point shift{1e12L,-1e12L};
        CHECK(cp_geometry::norm2(cp_geometry::polygon_centroid({a+shift,b+shift,c+shift})-(center+shift)) < 1e-12L);
    }
    // Concave L: signed fan triangles must give the area centroid, not the vertex average.
    auto l_center = cp_geometry::polygon_centroid({{0,0},{3,0},{3,1},{1,1},{1,3},{0,3}});
    CHECK(cp_geometry::near(l_center, {1.1L,1.1L}));
    bool threw = false;
    try { (void)cp_geometry::polygon_centroid({{0,0},{1,0},{2,0}}); }
    catch (const std::invalid_argument&) { threw = true; }
    CHECK(threw);

    CHECK(close(cp_geometry::segment_segment_distance({0,0},{2,0},{1,-1},{1,1}), 0));
    CHECK(close(cp_geometry::segment_segment_distance({0,0},{1,0},{2,1},{2,3}), std::sqrt(Real(2))));
    CHECK(close(cp_geometry::segment_segment_distance({0,0},{0,0},{3,4},{3,4}), 5));
    for (int iteration = 0; iteration < 1000; ++iteration) {
        Point a = random_point(), b = random_point(), p = random_point();
        if (cp_geometry::near(a,b)) continue;
        cp_geometry::Line line(a,b);
        Point reflected = cp_geometry::reflection(line,p);
        CHECK(cp_geometry::near(cp_geometry::reflection(line,reflected),p));
        CHECK(cp_geometry::near(cp_geometry::projection(line,p),(p+reflected)/2));
        Point c = random_point(), d = random_point();
        Real distance = cp_geometry::segment_segment_distance(a,b,c,d);
        CHECK(close(distance,cp_geometry::segment_segment_distance(d,c,b,a)));
        CHECK(close(distance,cp_geometry::segment_segment_distance(a,b,d,c)));
    }
    std::cout << "geometry algorithms: 3000 closest-pair and 3000 Minkowski brute comparisons, degeneracies, centroid/reflection/distance passed\n";
}
