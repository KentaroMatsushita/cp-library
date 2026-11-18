#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

#define INT(value) int value
#define IN(...) caller_input(__VA_ARGS__)
#define OUT(...) caller_output(__VA_ARGS__)
#define ON 54321
#define min(first, second) ((first) < (second) ? (first) : (second))
#define max(first, second) ((first) > (second) ? (first) : (second))
#define all(value) (value).begin(), (value).end()
#define rep(index, count) for (int index = 0; index < (count); ++index)
#include "geometry/power_diagram.hpp"
static_assert(min(3, 4) == 3 && max(3, 4) == 4 && ON == 54321);
#ifndef INT
#error INT must be preserved
#endif
#ifndef IN
#error IN must be preserved
#endif
#ifndef OUT
#error OUT must be preserved
#endif
#undef min
#undef max
#undef ON

using namespace cp_geometry;

static void close(Real actual, Real expected, Real tolerance = 1e-8L) {
    if (std::abs(actual-expected) > tolerance) {
        std::cerr << "actual=" << double(actual) << " expected=" << double(expected) << '\n';
        std::abort();
    }
}

static Real stable_area(const std::vector<Point>& polygon) {
    if (polygon.size() < 3) return 0;
    Real result = 0;
    for (std::size_t index = 1; index+1 < polygon.size(); ++index)
        result += cross(polygon[index]-polygon[0], polygon[index+1]-polygon[0])/2;
    return std::abs(result);
}

// Independent reference: sequential clipping with unnormalized coefficients.
static std::vector<Point> brute_cell(const std::vector<PowerSite>& sites, std::size_t selected,
                                     Point lower, Point upper) {
    std::vector<Point> polygon{lower, {upper.x, lower.y}, upper, {lower.x, upper.y}};
    for (std::size_t other = 0; other < sites.size(); ++other) {
        if (other == selected) continue;
        const auto first = sites[selected], second = sites[other];
        const Point normal = (second.point-first.point)*2;
        const Real constant = norm2(second.point)-norm2(first.point)+first.weight-second.weight;
        if (normal.x == 0 && normal.y == 0) {
            if (constant < 0 || (constant == 0 && other < selected)) return {};
            continue;
        }
        auto value = [&](Point point) { return dot(normal, point)-constant; };
        std::vector<Point> clipped;
        for (std::size_t index = 0; index < polygon.size(); ++index) {
            const Point start = polygon[index], finish = polygon[(index+1)%polygon.size()];
            const Real start_value = value(start), finish_value = value(finish);
            if (start_value <= 0) clipped.push_back(start);
            if ((start_value < 0 && finish_value > 0) || (start_value > 0 && finish_value < 0))
                clipped.push_back(start+(finish-start)*(start_value/(start_value-finish_value)));
        }
        polygon = std::move(clipped);
    }
    return polygon;
}

// Independent numerical integration of a UNION of disks clipped by the box.
static Real numeric_union_area(const std::vector<Circle>& circles, Point lower, Point upper) {
    constexpr int slices = 24000;
    const Real width = (upper.x-lower.x)/slices;
    Real result = 0;
    for (int slice = 0; slice < slices; ++slice) {
        const Real x = lower.x+(slice+Real(0.5))*width;
        std::vector<std::pair<Real, Real>> intervals;
        for (const auto& circle : circles) {
            const Real dx = std::abs(x-circle.center.x);
            if (dx >= circle.radius) continue;
            const Real dy = std::sqrt((circle.radius-dx)*(circle.radius+dx));
            const Real low = std::fmax(lower.y, circle.center.y-dy);
            const Real high = std::fmin(upper.y, circle.center.y+dy);
            if (low < high) intervals.push_back({low, high});
        }
        std::sort(intervals.begin(), intervals.end());
        Real high = lower.y;
        for (const auto& interval : intervals) {
            if (interval.second > high) {
                result += (interval.second-std::fmax(interval.first, high))*width;
                high = interval.second;
            }
        }
    }
    return result;
}

int main() {
    const Point lower{-6, -5}, upper{6, 5};
    const Real pi = std::acos(Real(-1));
    assert(power_diagram({}, lower, upper).empty());
    auto region = power_cell({{{2, 3}, -7}}, 0);
    assert(region.kind == RegionKind::UNBOUNDED);
    close(region.witness.x, 2); close(region.witness.y, 3);
    close(power_value({{2, 3}, -7}, {5, 7}), 32);
    close(stable_area(power_cell({{{2, 3}, -7}}, 0, lower, upper)), 120);
    auto cells = power_diagram({{{-1, 0}, 4}, {{1, 0}, 0}}, {-2, -2}, {2, 2});
    close(stable_area(cells[0]), 12); close(stable_area(cells[1]), 4);
    auto halfplane = power_bisector({{-1, 0}, 4}, {{1, 0}, 0});
    close(halfplane.a, 1); close(halfplane.b, 0); close(halfplane.c, 1);
    halfplane = power_bisector({{0, 0}, 3}, {{0, 0}, 4});
    assert(halfplane.a == 0 && halfplane.b == 0 && halfplane.c < 0);
    cells = power_diagram({{{0, 0}, 3}, {{0, 0}, 3}, {{0, 0}, 1}}, lower, upper);
    close(stable_area(cells[0]), 120); assert(cells[1].empty() && cells[2].empty());
    cells = power_diagram({{{0, 0}, 1}, {{0, 0}, 3}, {{0, 0}, 3}}, lower, upper);
    assert(cells[0].empty() && cells[2].empty()); close(stable_area(cells[1]), 120);
    assert(power_cell({{{0, 0}, 1}, {{0, 0}, 3}}, 0).kind == RegionKind::EMPTY);
    std::vector<PowerSite> point_cell{{{0, 0}, 0}, {{2, 0}, 4}, {{-2, 0}, 4},
                                      {{0, 2}, 4}, {{0, -2}, 4}};
    region = power_cell(point_cell, 0);
    assert(region.kind == RegionKind::DEGENERATE);
    close(region.witness.x, 0); close(region.witness.y, 0);
    assert(power_cell(point_cell, 0, lower, upper).empty());
    point_cell[0].weight = -1;
    assert(power_cell(point_cell, 0).kind == RegionKind::EMPTY);
    assert(power_cell(point_cell, 0, lower, upper).empty());
    region = power_cell({{{0, 0}, 0}, {{2, 0}, 4}, {{-2, 0}, 4}}, 0);
    assert(region.kind == RegionKind::UNBOUNDED); close(region.witness.x, 0);
    std::vector<PowerSite> bounded{{{0, 0}, 0}, {{2, 0}, 0}, {{-2, 0}, 0},
                                   {{0, 2}, 0}, {{0, -2}, 0}};
    region = power_cell(bounded, 0);
    assert(region.kind == RegionKind::POLYGON); close(stable_area(region.polygon), 4);
    close(stable_area(power_cell(bounded, 0, lower, upper)), 4);

    bool rejected = false;
    try { (void)power_cell({}, 0, lower, upper); }
    catch (const std::out_of_range&) { rejected = true; }
    assert(rejected);
    rejected = false;
    try { (void)power_diagram({}, {0, 0}, {0, 1}); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
    rejected = false;
    try { (void)power_diagram({}, lower, {std::numeric_limits<Real>::infinity(), 1}); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
    rejected = false;
    try { (void)power_cell({{{0, 0}, std::numeric_limits<Real>::quiet_NaN()}}, 0); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);

    const std::vector<PowerSite> translated_sites{{{1e12L-1, -1e12L}, 4}, {{1e12L+1, -1e12L}, 0}};
    auto translated_cell = power_cell(translated_sites, 0, {1e12L-2, -1e12L-2}, {1e12L+2, -1e12L+2});
    close(stable_area(translated_cell), 12);
    halfplane = power_bisector(translated_sites[0], translated_sites[1]);
    close(halfplane.c-1e12L, 1);

    std::mt19937_64 random(82035719);
    auto coordinate = [&] { return Real(int(random()%81)-40)/5; };
    for (int iteration = 0; iteration < 800; ++iteration) {
        std::vector<PowerSite> sites;
        const int count = 1+int(random()%13);
        for (int index = 0; index < count; ++index) {
            if (index && random()%7 == 0) sites.push_back(sites[random()%index]);
            else sites.push_back({{coordinate(), coordinate()}, Real(int(random()%101)-50)/3});
        }
        cells = power_diagram(sites, lower, upper);
        assert(cells.size() == sites.size());
        Real total_area = 0;
        for (std::size_t index = 0; index < sites.size(); ++index) {
            const Real area = stable_area(cells[index]);
            total_area += area;
            close(area, stable_area(brute_cell(sites, index, lower, upper)), 1e-7L);
            for (Point point : cells[index]) {
                assert(point.x >= lower.x-1e-8L && point.x <= upper.x+1e-8L);
                assert(point.y >= lower.y-1e-8L && point.y <= upper.y+1e-8L);
                for (const auto& competitor : sites)
                    assert(power_value(sites[index], point) <= power_value(competitor, point)+1e-7L);
            }
            for (int query = 0; query < 5; ++query) {
                Point point{coordinate(), coordinate()};
                for (std::size_t other = 0; other < sites.size(); ++other) {
                    auto inequality = power_bisector(sites[index], sites[other]);
                    const Real residual = inequality.a*point.x+inequality.b*point.y-inequality.c;
                    const Real difference = power_value(sites[index], point)-power_value(sites[other], point);
                    if (std::abs(residual) > 1e-8L && std::abs(difference) > 1e-8L)
                        assert((residual > 0) == (difference > 0));
                }
            }
        }
        close(total_area, 120, 1e-7L);
    }

    const Point galaxy_lower{-10000, -10000}, galaxy_upper{10000, 10000};
    close(circle_territory_area({Circle({1, 1}, 2), Circle({4, 5}, 2),
          Circle({-2, 0}, 1), Circle({5, -7}, 3)}, 0, galaxy_lower, galaxy_upper),
          Real(718205)/16128-4*pi, 1e-7L);
    close(circle_territory_area({Circle({0, 0}, 1000), Circle({0, 4000}, 1000)},
          0, galaxy_lower, galaxy_upper), 236858407.346410206761537L, 1e-6L);
    close(circle_territory_area({Circle({1, -1}, 1), Circle({2, -1}, 2),
          Circle({0, 0}, 2), Circle({-3, -2}, 4)}, 0, galaxy_lower, galaxy_upper), 0, 1e-8L);
    close(circle_territory_area({Circle({15, -3}, 7), Circle({13, -4}, 2), Circle({14, -8}, 3),
          Circle({22, 1}, 4), Circle({42, 7}, 1), Circle({-1, 0}, 2), Circle({-3, 0}, 1),
          Circle({0, -15}, 5)}, 0, galaxy_lower, galaxy_upper), 163.354738705171682L, 1e-7L);
    close(circle_territory_area({Circle({0, 0}, 1)}, 0, {-2, -2}, {2, 2}), 16-pi);
    close(circle_territory_area({Circle({0, 0}, 0)}, 0, {-2, -2}, {2, 2}), 16);
    assert(circle_territory_area({Circle({0, 0}, 1), Circle({0, 0}, 1)}, 1, lower, upper) == 0);
    for (int iteration = 0; iteration < 20; ++iteration) {
        std::vector<Circle> circles;
        for (int index = 0; index < 5; ++index)
            circles.emplace_back(Point{coordinate()/2, coordinate()/2}, Real(random()%500)/100);
        Real uncovered = 0;
        for (std::size_t index = 0; index < circles.size(); ++index)
            uncovered += circle_territory_area(circles, index, lower, upper);
        close(uncovered, 120-numeric_union_area(circles, lower, upper), 3e-4L);
    }
    std::cout << "power diagram: 800 clipping/partition comparisons, four Galaxy War samples, 20 union integrations passed\n";
}
