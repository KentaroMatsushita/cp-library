#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

#define INT(value) int value
#define IN(...) input_macro_must_be_preserved(__VA_ARGS__)
#define OUT(...) output_macro_must_be_preserved(__VA_ARGS__)
#define ON 123456
#define min(first, second) ((first) < (second) ? (first) : (second))
#define max(first, second) ((first) > (second) ? (first) : (second))
#define all(value) (value).begin(), (value).end()
#define rep(index, count) for (int index = 0; index < (count); ++index)
#include "geometry/circle.hpp"
static_assert(min(7, 8) == 7 && max(7, 8) == 8 && ON == 123456);
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

static void close(Real actual, Real expected, Real tolerance = 1e-9L) {
    if (std::abs(actual-expected) > tolerance) {
        std::cerr << "actual " << double(actual) << ", expected " << double(expected) << '\n';
        std::abort();
    }
}

static void on_circle(Point point, const Circle& circle) {
    close(std::hypot(point.x-circle.center.x, point.y-circle.center.y), circle.radius, 1e-8L);
}

// Independent vertical integration: polygon edge crossings are paired into intervals.
static Real numeric_polygon_area(const Circle& circle, const std::vector<Point>& polygon) {
    constexpr int slices = 18000;
    const Real left = circle.center.x-circle.radius;
    const Real width = 2*circle.radius/slices;
    Real result = 0;
    std::vector<Real> crossings;
    for (int slice = 0; slice < slices; ++slice) {
        const Real x = left+(slice+Real(0.5))*width;
        const Real offset = std::sqrt(std::fmax(Real(0),
            circle.radius*circle.radius-(x-circle.center.x)*(x-circle.center.x)));
        const Real circle_low = circle.center.y-offset, circle_high = circle.center.y+offset;
        crossings.clear();
        for (std::size_t index = 0; index < polygon.size(); ++index) {
            Point first = polygon[index], second = polygon[(index+1)%polygon.size()];
            if (first.x > second.x) std::swap(first, second);
            if (first.x <= x && x < second.x)
                crossings.push_back(first.y+(second.y-first.y)*(x-first.x)/(second.x-first.x));
        }
        std::sort(crossings.begin(), crossings.end());
        assert(crossings.size()%2 == 0);
        for (std::size_t index = 1; index < crossings.size(); index += 2)
            result += std::fmax(Real(0), std::fmin(circle_high, crossings[index])-
                                             std::fmax(circle_low, crossings[index-1]))*width;
    }
    return result;
}

static Real numeric_circle_area(const Circle& first, const Circle& second) {
    constexpr int slices = 18000;
    const Real left = std::fmax(first.center.x-first.radius, second.center.x-second.radius);
    const Real right = std::fmin(first.center.x+first.radius, second.center.x+second.radius);
    if (right <= left) return 0;
    const Real width = (right-left)/slices;
    Real result = 0;
    for (int slice = 0; slice < slices; ++slice) {
        const Real x = left+(slice+Real(0.5))*width;
        const Real first_offset = std::sqrt(std::fmax(Real(0),
            first.radius*first.radius-(x-first.center.x)*(x-first.center.x)));
        const Real second_offset = std::sqrt(std::fmax(Real(0),
            second.radius*second.radius-(x-second.center.x)*(x-second.center.x)));
        result += width*std::fmax(Real(0),
            std::fmin(first.center.y+first_offset, second.center.y+second_offset)-
            std::fmax(first.center.y-first_offset, second.center.y-second_offset));
    }
    return result;
}

int main() {
    const Real pi = std::acos(Real(-1));
    bool rejected = false;
    try { Circle invalid({0, 0}, -1); } catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
    rejected = false;
    try { circumcircle({0, 0}, {1, 1}, {2, 2}); }
    catch (const std::invalid_argument&) { rejected = true; }
    assert(rejected);
    const Circle unit({0, 0}, 1);
    auto diameter = circle_from_diameter({-2, 1}, {4, 1});
    close(diameter.center.x, 1); close(diameter.center.y, 1); close(diameter.radius, 3);
    auto outer = circumcircle({0, 0}, {4, 0}, {0, 3});
    close(outer.center.x, 2); close(outer.center.y, 1.5L); close(outer.radius, 2.5L);
    auto inner = incircle({0, 0}, {4, 0}, {0, 3});
    close(inner.center.x, 1); close(inner.center.y, 1); close(inner.radius, 1);
    auto reversed_inner = incircle({0, 3}, {4, 0}, {0, 0});
    close(reversed_inner.center.x, 1); close(reversed_inner.center.y, 1);

    const Line horizontal({-3, 0}, {3, 0});
    auto points = line_circle_intersections(horizontal, unit);
    assert(points.size() == 2); close(points[0].x, -1); close(points[1].x, 1);
    assert(line_circle_intersections(Line({0, 2}, {1, 2}), unit).empty());
    points = line_circle_intersections(Line({0, 1}, {1, 1}), unit);
    assert(points.size() == 1); close(points[0].x, 0); close(points[0].y, 1);
    assert(segment_circle_intersections({-0.5L, 0}, {0.5L, 0}, unit).empty());
    assert(segment_circle_intersections({-1, 0}, {1, 0}, unit).size() == 2);
    assert(segment_circle_intersections({1, 0}, {1, 0}, unit).size() == 1);
    assert(segment_circle_intersections({0, 0}, {0, 0}, unit).empty());
    assert(segment_circle_intersections({2, 0}, {-2, 0}, unit).front().x == 1);
    assert(circle_intersections(unit, unit).coincident);
    assert(circle_intersections(Circle({0, 0}, 0), Circle({0, 0}, 0)).points.size() == 1);
    assert(!circle_intersections(Circle({0, 0}, 0), Circle({0, 0}, 0)).coincident);
    assert(circle_intersections(Circle({1, 0}, 0), unit).points.size() == 1);
    assert(circle_intersections(unit, Circle({2, 0}, 1)).points.size() == 1);
    assert(circle_intersections(unit, Circle({1, 0}, 2)).points.size() == 1);
    assert(circle_intersections(unit, Circle({3, 0}, 1)).points.empty());
    assert(circle_intersections(unit, Circle({0, 0}, 2)).points.empty());
    close(circle_overlap_area(unit, unit), pi);
    close(circle_overlap_area(unit, Circle({2, 0}, 1)), 0);
    close(circle_overlap_area(unit, Circle({1, 0}, 1)), 2*pi/3-std::sqrt(Real(3))/2);
    close(circle_overlap_area(unit, Circle({2, 0}, 2)), 1.40306643968573875104L);
    close(circle_overlap_area(Circle({1, 0}, 1), Circle({0, 0}, 3)), pi);
    close(circle_overlap_area(Circle({0, 0}, 0), unit), 0);

    std::vector<Point> enclosing{{-2, -2}, {2, -2}, {2, 2}, {-2, 2}};
    close(circle_polygon_intersection_area(unit, enclosing), pi);
    std::vector<Point> small{{0, 0}, {0.5L, 0}, {0.5L, 0.5L}, {0, 0.5L}};
    close(circle_polygon_intersection_area(unit, small), 0.25L);
    std::vector<Point> quadrant{{0, 0}, {3, 0}, {0, 3}};
    close(circle_polygon_intersection_area(unit, quadrant), pi/4);
    std::vector<Point> concave{{-2, -2}, {2, -2}, {2, 0}, {0, 0}, {0, 2}, {-2, 2}};
    close(circle_polygon_intersection_area(unit, concave), 3*pi/4);
    std::reverse(concave.begin(), concave.end());
    close(circle_polygon_intersection_area(unit, concave), 3*pi/4);
    close(circle_polygon_intersection_area(Circle({0, 0}, 5), {{1, 1}, {4, 1}, {5, 5}}),
          4.639858417607L, 1e-10L);
    close(circle_polygon_intersection_area(Circle({0, 0}, 5),
          {{0, 0}, {-3, -6}, {1, -3}, {5, -4}}), 11.787686807576L, 1e-10L);
    close(circle_polygon_intersection_area(unit, {{4, 4}, {5, 4}, {5, 5}, {4, 5}}), 0);
    close(circle_polygon_intersection_area(unit, {}), 0);
    close(circle_polygon_intersection_area(unit, {{0, 0}, {1, 0}}), 0);
    std::vector<Point> touching{{-2, 1}, {2, 1}, {2, 2}, {-2, 2}};
    close(circle_polygon_intersection_area(unit, touching), 0);
    assert(tangent_points(unit, {0, 0}).empty());
    assert(tangent_points(Circle({0, 0}, 1e-15L), {0, 0}).empty());
    assert(tangent_points(unit, {1, 0}).size() == 1);
    points = tangent_points(unit, {2, 0});
    assert(points.size() == 2);
    for (Point point : points) { on_circle(point, unit); close(dot(point, point-Point{2, 0}), 0); }
    assert(common_tangents(unit, unit).coincident);
    assert(common_tangents(Circle({0, 0}, 0), Circle({0, 0}, 0)).coincident);
    assert(common_tangents(unit, Circle({3, 0}, 1)).tangents.size() == 4);
    assert(common_tangents(unit, Circle({2, 0}, 1)).tangents.size() == 3);
    assert(common_tangents(unit, Circle({1, 0}, 1)).tangents.size() == 2);
    assert(common_tangents(unit, Circle({1, 0}, 2)).tangents.size() == 1);
    assert(common_tangents(unit, Circle({0, 0}, 2)).tangents.empty());
    assert(common_tangents(Circle({0, 0}, 0), Circle({3, 0}, 0)).tangents.size() == 1);
    assert(common_tangents(Circle({0, 0}, 0), Circle({3, 0}, 1)).tangents.size() == 2);

    std::mt19937_64 random(972541);
    auto coordinate = [&] { return Real(int(random()%2001)-1000)/200; };
    for (int iteration = 0; iteration < 1500; ++iteration) {
        Circle first({coordinate(), coordinate()}, Real(random()%800+1)/200);
        Circle second({coordinate(), coordinate()}, Real(random()%800+1)/200);
        auto intersection = circle_intersections(first, second);
        assert(!intersection.coincident);
        for (Point point : intersection.points) { on_circle(point, first); on_circle(point, second); }
        const Real distance = std::hypot(first.center.x-second.center.x, first.center.y-second.center.y);
        const std::size_t expected = distance > first.radius+second.radius ||
                                    distance < std::abs(first.radius-second.radius) ? 0 : 2;
        if (std::abs(distance-first.radius-second.radius) > 1e-8L &&
            std::abs(distance-std::abs(first.radius-second.radius)) > 1e-8L)
            assert(intersection.points.size() == expected);
        close(circle_overlap_area(first, second), circle_overlap_area(second, first), 1e-8L);
        const auto tangents = common_tangents(first, second);
        for (const auto& tangent : tangents.tangents) {
            on_circle(tangent.first_contact, first); on_circle(tangent.second_contact, second);
            close(cross(tangent.line.direction(), tangent.first_contact-tangent.line.a), 0, 1e-8L);
            close(cross(tangent.line.direction(), tangent.second_contact-tangent.line.a), 0, 1e-8L);
            close(dot(tangent.line.direction(), tangent.first_contact-first.center), 0, 1e-8L);
            close(dot(tangent.line.direction(), tangent.second_contact-second.center), 0, 1e-8L);
        }
        if (iteration < 80) close(circle_overlap_area(first, second),
                                  numeric_circle_area(first, second), 2e-4L);
        Point first_endpoint{coordinate(), coordinate()}, second_endpoint{coordinate(), coordinate()};
        for (Point point : segment_circle_intersections(first_endpoint, second_endpoint, first)) {
            on_circle(point, first);
            close(cross(second_endpoint-first_endpoint, point-first_endpoint), 0, 1e-8L);
            assert(dot(point-first_endpoint, point-second_endpoint) <= 1e-8L);
        }
    }
    // Alternating radii produce simple concave polygons, with arbitrary circle centers.
    for (int iteration = 0; iteration < 55; ++iteration) {
        const int count = 6+int(random()%9)*2;
        std::vector<Point> polygon;
        for (int index = 0; index < count; ++index) {
            const Real angle = 2*pi*index/count;
            const Real radius = (index%2 ? Real(0.5) : Real(3)) + Real(random()%100)/100;
            polygon.push_back({radius*std::cos(angle), radius*std::sin(angle)});
        }
        Circle circle({coordinate()/2, coordinate()/2}, Real(random()%600+1)/200);
        const Real area = circle_polygon_intersection_area(circle, polygon);
        close(area, numeric_polygon_area(circle, polygon), 3e-4L);
        assert(area <= circle_area(circle)+1e-9L);
        assert(area <= std::abs(area2(polygon))/2+1e-9L);
        const Point shift{1234, -918};
        for (auto& point : polygon) point = point+shift;
        close(area, circle_polygon_intersection_area(Circle(circle.center+shift, circle.radius), polygon), 1e-8L);
        std::reverse(polygon.begin(), polygon.end());
        close(area, circle_polygon_intersection_area(Circle(circle.center+shift, circle.radius), polygon), 1e-8L);
    }
}
