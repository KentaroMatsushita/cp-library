#pragma once
#include "detail/macro_push.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <vector>
#include "base.hpp"

namespace cp_geometry {

struct Circle {
    Point center;
    Real radius;
    Circle(Point circle_center, Real circle_radius)
        : center(circle_center), radius(circle_radius) {
        if (!std::isfinite(center.x) || !std::isfinite(center.y) ||
            !std::isfinite(radius) || radius < 0)
            throw std::invalid_argument("circle needs a finite center and nonnegative radius");
    }
};

inline Real circle_area(const Circle& circle) {
    return std::acos(Real(-1)) * circle.radius * circle.radius;
}

inline Circle circle_from_diameter(Point first, Point second) {
    return Circle((first + second) / 2, std::hypot(first.x-second.x, first.y-second.y) / 2);
}

// Three noncollinear points; the computation is translated to the first point.
inline Circle circumcircle(Point first, Point second, Point third) {
    const Point second_delta = second-first, third_delta = third-first;
    const Real determinant = 2 * cross(second_delta, third_delta);
    if (determinant == 0)
        throw std::invalid_argument("circumcircle needs three noncollinear points");
    const Point offset = (perpendicular(second_delta)*norm2(third_delta) -
                          perpendicular(third_delta)*norm2(second_delta)) / determinant;
    const Point center = first + offset;
    return Circle(center, std::hypot(offset.x, offset.y));
}

inline Circle incircle(Point first, Point second, Point third) {
    const Real first_weight = std::hypot(second.x-third.x, second.y-third.y);
    const Real second_weight = std::hypot(first.x-third.x, first.y-third.y);
    const Real third_weight = std::hypot(first.x-second.x, first.y-second.y);
    const Real perimeter = first_weight+second_weight+third_weight;
    const Real twice_area = std::abs(cross(second-first, third-first));
    if (twice_area == 0 || perimeter == 0)
        throw std::invalid_argument("incircle needs three noncollinear points");
    return Circle((first*first_weight+second*second_weight+third*third_weight)/perimeter,
                  twice_area/perimeter);
}

// Boundary intersections, ordered along line.direction(); tangent points occur once.
inline std::vector<Point> line_circle_intersections(const Line& line, const Circle& circle) {
    const Point delta = line.direction();
    const Real length = std::hypot(delta.x, delta.y);
    const Point unit = delta/length;
    const Point center_delta = circle.center-line.a;
    const Real signed_distance = cross(unit, center_delta);
    const Real distance = std::abs(signed_distance);
    if (distance > circle.radius+eps) return {};
    const Point foot = line.a+unit*dot(center_delta, unit);
    if (std::abs(distance-circle.radius) <= eps) return {foot};
    const Real offset = std::sqrt(std::fmax(Real(0),
        (circle.radius-distance)*(circle.radius+distance)));
    return {foot-unit*offset, foot+unit*offset};
}

// Segment endpoints may coincide. Points are ordered from first toward second.
inline std::vector<Point> segment_circle_intersections(Point first, Point second,
                                                      const Circle& circle) {
    const Point delta = second-first;
    const Real length = std::hypot(delta.x, delta.y);
    if (length == 0) {
        const Real distance = std::hypot(first.x-circle.center.x, first.y-circle.center.y);
        return std::abs(distance-circle.radius) <= eps ? std::vector<Point>{first}
                                                     : std::vector<Point>{};
    }
    const Point unit = delta/length;
    const Point center_delta = circle.center-first;
    const Real distance = std::abs(cross(unit, center_delta));
    if (distance > circle.radius+eps) return {};
    const Real along = dot(center_delta, unit);
    const Real offset = std::abs(distance-circle.radius) <= eps ? 0 :
        std::sqrt(std::fmax(Real(0), (circle.radius-distance)*(circle.radius+distance)));
    std::vector<Point> result;
    for (Real coordinate : {along-offset, along+offset}) {
        if (coordinate < -eps || coordinate > length+eps) continue;
        coordinate = std::fmax(Real(0), std::fmin(length, coordinate));
        const Point point = first+unit*coordinate;
        if (result.empty() || !near(result.back(), point)) result.push_back(point);
    }
    return result;
}

struct CircleIntersection {
    // Exactly equal circles with positive radius have infinitely many intersections.
    bool coincident = false;
    std::vector<Point> points;
};

inline CircleIntersection circle_intersections(const Circle& first, const Circle& second) {
    const Point delta = second.center-first.center;
    const Real distance = std::hypot(delta.x, delta.y);
    if (distance == 0) {
        if (first.radius != second.radius) return {};
        if (first.radius == 0) return {false, {first.center}};
        return {true, {}};
    }
    if (first.radius == 0) {
        return std::abs(distance-second.radius) <= eps
            ? CircleIntersection{false, {first.center}} : CircleIntersection{};
    }
    if (second.radius == 0) {
        return std::abs(distance-first.radius) <= eps
            ? CircleIntersection{false, {second.center}} : CircleIntersection{};
    }
    const Real radius_sum = first.radius+second.radius;
    const Real radius_difference = first.radius-second.radius;
    if (distance > radius_sum+eps || distance < std::abs(radius_difference)-eps) return {};
    const Point unit = delta/distance;
    const Real along = (distance + radius_difference*radius_sum/distance)/2;
    const Point foot = first.center+unit*along;
    if (std::abs(distance-radius_sum) <= eps ||
        std::abs(distance-std::abs(radius_difference)) <= eps) return {false, {foot}};
    const Real height = std::sqrt(std::fmax(Real(0),
        (first.radius-along)*(first.radius+along)));
    const Point normal = perpendicular(unit)*height;
    return {false, {foot-normal, foot+normal}};
}

// Area of the intersection of the two closed disks. O(1).
inline Real circle_overlap_area(const Circle& first, const Circle& second) {
    const Real radius_min = std::fmin(first.radius, second.radius);
    if (radius_min == 0) return 0;
    const Real distance = std::hypot(second.center.x-first.center.x,
                                   second.center.y-first.center.y);
    const Real radius_sum = first.radius+second.radius;
    if (distance >= radius_sum) return 0;
    if (distance <= std::abs(first.radius-second.radius))
        return std::acos(Real(-1))*radius_min*radius_min;
    const Real along_first = (distance +
        (first.radius-second.radius)*radius_sum/distance)/2;
    const Real along_second = distance-along_first;
    const Real height = std::sqrt(std::fmax(Real(0),
        (first.radius-along_first)*(first.radius+along_first)));
    const Real first_angle = std::atan2(height, along_first);
    const Real second_angle = std::atan2(height, along_second);
    const Real result = first.radius*first.radius*
                            (first_angle-std::sin(2*first_angle)/2) +
                        second.radius*second.radius*
                            (second_angle-std::sin(2*second_angle)/2);
    return std::fmax(Real(0), std::fmin(result, std::acos(Real(-1))*radius_min*radius_min));
}

namespace circle_detail {
// Signed contribution of a directed polygon edge translated to the circle center.
// Splitting at circle intersections leaves either a triangle or a circular sector.
inline Real edge_disk_area(Point first, Point second, Real radius) {
    const Point delta = second-first;
    const Real length = std::hypot(delta.x, delta.y);
    if (length == 0) return 0;
    const Point unit = delta/length;
    const Real along = -dot(first, unit);
    const Real distance = std::abs(cross(unit, first));
    std::array<Real, 4> cuts{0, length, 0, 0};
    std::size_t cut_count = 2;
    if (distance < radius) {
        const Real offset = std::sqrt((radius-distance)*(radius+distance));
        for (Real coordinate : {along-offset, along+offset})
            if (0 < coordinate && coordinate < length) cuts[cut_count++] = coordinate;
    }
    std::sort(cuts.begin(), cuts.begin()+cut_count);
    Real result = 0;
    for (std::size_t index = 1; index < cut_count; ++index) {
        const Point start = first+unit*cuts[index-1];
        const Point finish = first+unit*cuts[index];
        const Point middle = first+unit*((cuts[index-1]+cuts[index])/2);
        if (std::hypot(middle.x, middle.y) < radius)
            result += cross(start, finish)/2;
        else
            result += radius*radius*std::atan2(cross(start, finish), dot(start, finish))/2;
    }
    return result;
}
} // namespace circle_detail

// A simple polygon, either orientation, including concave polygons. O(N).
// Empty polygons and fewer than three vertices have area zero.
inline Real circle_polygon_intersection_area(const Circle& circle,
                                             const std::vector<Point>& polygon) {
    if (circle.radius == 0 || polygon.size() < 3) return 0;
    Real result = 0, correction = 0;
    for (std::size_t index = 0; index < polygon.size(); ++index) {
        const Real contribution = circle_detail::edge_disk_area(
            polygon[index]-circle.center, polygon[(index+1)%polygon.size()]-circle.center,
            circle.radius);
        const Real adjusted = contribution-correction;
        const Real updated = result+adjusted;
        correction = (updated-result)-adjusted;
        result = updated;
    }
    return std::abs(result);
}

// Contact points on a circle of tangent lines through point. O(1).
// For a radius-zero circle, the only contact point is its center.
inline std::vector<Point> tangent_points(const Circle& circle, Point point) {
    if (circle.radius == 0) return {circle.center};
    const Point delta = point-circle.center;
    const Real distance = std::hypot(delta.x, delta.y);
    if (distance == 0) return {};
    if (distance < circle.radius-eps) return {};
    const Point unit = delta/distance;
    if (std::abs(distance-circle.radius) <= eps) return {circle.center+unit*circle.radius};
    const Real cosine = circle.radius/distance;
    const Real sine = std::sqrt(std::fmax(Real(0), (1-cosine)*(1+cosine)));
    const Point foot = circle.center+unit*(circle.radius*cosine);
    const Point offset = perpendicular(unit)*(circle.radius*sine);
    return {foot-offset, foot+offset};
}

struct CommonTangent {
    Point first_contact, second_contact;
    Line line;
};
struct CircleTangents {
    // Equal concentric circles, including two equal point circles: infinitely many lines.
    bool coincident = false;
    std::vector<CommonTangent> tangents;
};

// All distinct inner and outer common tangents, including touching/point circles. O(1).
inline CircleTangents common_tangents(const Circle& first, const Circle& second) {
    const Point delta = second.center-first.center;
    const Real distance = std::hypot(delta.x, delta.y);
    if (distance == 0) return {first.radius == second.radius, {}};
    const Point unit = delta/distance;
    CircleTangents result;
    for (Real second_sign : {Real(-1), Real(1)}) {
        const Real radius_difference = first.radius-second_sign*second.radius;
        if (std::abs(radius_difference) > distance+eps) continue;
        const Real cosine = std::fmax(Real(-1), std::fmin(Real(1), radius_difference/distance));
        const Real sine = std::sqrt(std::fmax(Real(0), (1-cosine)*(1+cosine)));
        for (Real turn : {Real(-1), Real(1)}) {
            const Point normal = unit*cosine+perpendicular(unit)*(sine*turn);
            const Point first_contact = first.center+normal*first.radius;
            const Point second_contact = second.center+normal*(second.radius*second_sign);
            const Point direction = perpendicular(normal);
            bool duplicate = false;
            for (const auto& tangent : result.tangents) {
                if (std::abs(cross(direction, tangent.line.direction())) <= eps &&
                    std::abs(cross(direction, tangent.line.a-first_contact)) <= eps) {
                    duplicate = true;
                    break;
                }
            }
            if (!duplicate) result.tangents.push_back(
                {first_contact, second_contact, Line(first_contact, first_contact+direction)});
        }
    }
    return result;
}

} // namespace cp_geometry
#include "detail/macro_pop.hpp"
