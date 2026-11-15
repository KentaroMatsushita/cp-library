#pragma once
#include "detail/macro_push.hpp"
#include "convex.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <deque>
#include <limits>
#include <optional>
#include <stdexcept>
#include <vector>

namespace cp_geometry {

struct ClosestPairResult {
    std::size_t first, second; // Indices in the original input, first != second.
    Real squared_distance;
};

// O(N log N) time, O(N) extra memory. Empty optional for fewer than two points.
// Ties are arbitrary. Coordinates and squared distances must fit in Real.
inline std::optional<ClosestPairResult> closest_pair(const std::vector<Point>& points) {
    struct IndexedPoint { Point point; std::size_t index; };
    std::vector<IndexedPoint> sorted;
    sorted.reserve(points.size());
    for (std::size_t i = 0; i < points.size(); ++i) {
        if (!std::isfinite(points[i].x) || !std::isfinite(points[i].y))
            throw std::invalid_argument("closest_pair needs finite coordinates");
        sorted.push_back({points[i], i});
    }
    if (points.size() < 2) return std::nullopt;
    std::sort(sorted.begin(), sorted.end(), [](const auto& a, const auto& b) {
        return a.point.x != b.point.x ? a.point.x < b.point.x : a.point.y < b.point.y;
    });
    for (std::size_t i = 1; i < sorted.size(); ++i) {
        if (sorted[i-1].point.x == sorted[i].point.x && sorted[i-1].point.y == sorted[i].point.y)
            return ClosestPairResult{sorted[i-1].index, sorted[i].index, 0};
    }
    ClosestPairResult best{0, 1, norm2(points[0] - points[1])};
    auto consider = [&](const IndexedPoint& a, const IndexedPoint& b) {
        Real distance = norm2(a.point - b.point);
        if (distance < best.squared_distance) best = {a.index, b.index, distance};
    };
    auto by_y = [](const IndexedPoint& a, const IndexedPoint& b) {
        return a.point.y != b.point.y ? a.point.y < b.point.y : a.point.x < b.point.x;
    };
    std::vector<IndexedPoint> buffer(sorted.size());
    auto solve = [&](auto&& self, std::size_t left, std::size_t right) -> void {
        if (right - left <= 3) {
            for (std::size_t i = left; i < right; ++i)
                for (std::size_t j = left; j < i; ++j) consider(sorted[i], sorted[j]);
            std::sort(sorted.begin() + left, sorted.begin() + right, by_y);
            return;
        }
        std::size_t middle = left + (right - left) / 2;
        Real split_x = sorted[middle].point.x;
        self(self, left, middle);
        self(self, middle, right);
        std::merge(sorted.begin() + left, sorted.begin() + middle,
                   sorted.begin() + middle, sorted.begin() + right,
                   buffer.begin() + left, by_y);
        std::copy(buffer.begin() + left, buffer.begin() + right, sorted.begin() + left);
        // The same buffer now holds the y-sorted strip. Packing bounds the
        // number of earlier points within sqrt(best.squared_distance).
        std::size_t count = 0;
        for (std::size_t i = left; i < right; ++i) {
            Real dx = sorted[i].point.x - split_x;
            if (dx * dx >= best.squared_distance) continue;
            for (std::size_t j = count; j > 0; --j) {
                Real dy = sorted[i].point.y - buffer[left + j - 1].point.y;
                if (dy * dy >= best.squared_distance) break;
                consider(sorted[i], buffer[left + j - 1]);
            }
            buffer[left + count++] = sorted[i];
        }
    };
    solve(solve, 0, sorted.size());
    return best;
}

namespace algorithms_detail {
inline void rotate_bottommost(std::vector<Point>& polygon) {
    if (polygon.empty()) return;
    auto first = std::min_element(polygon.begin(), polygon.end(), [](Point a, Point b) {
        return a.y != b.y ? a.y < b.y : a.x < b.x;
    });
    std::rotate(polygon.begin(), first, polygon.end());
}

// The input is already a convex boundary in CCW order. Remove redundant
// collinear vertices in linear time, including the closing two corners.
inline std::vector<Point> strict_boundary(const std::vector<Point>& polygon) {
    std::deque<Point> boundary;
    for (Point point : polygon) {
        if (!boundary.empty() && near(boundary.back(), point)) continue;
        while (boundary.size() > 1 &&
               sign(cross(boundary.back() - boundary[boundary.size()-2], point - boundary.back())) <= 0)
            boundary.pop_back();
        boundary.push_back(point);
    }
    if (boundary.size() > 1 && near(boundary.front(), boundary.back())) boundary.pop_back();
    while (boundary.size() > 2) {
        if (sign(cross(boundary.back() - boundary[boundary.size()-2], boundary.front() - boundary.back())) <= 0) {
            boundary.pop_back();
        } else if (sign(cross(boundary.front() - boundary.back(), boundary[1] - boundary.front())) <= 0) {
            boundary.pop_front();
        } else break;
    }
    std::vector<Point> result(boundary.begin(), boundary.end());
    rotate_bottommost(result);
    return result;
}
} // namespace algorithms_detail

// O(N+M), including degeneracies. Inputs must be strict CCW convex boundaries
// without a repeated closing vertex, or an empty/one-point/two-point hull.
// Result is strict CCW, starting at its minimum (y,x) vertex.
inline std::vector<Point> minkowski_sum_convex(std::vector<Point> a, std::vector<Point> b) {
    if (a.empty() || b.empty()) return {};
    algorithms_detail::rotate_bottommost(a);
    algorithms_detail::rotate_bottommost(b);
    if (a.size() == 1) {
        for (auto& point : b) point = point + a[0];
        return b;
    }
    if (b.size() == 1) {
        for (auto& point : a) point = point + b[0];
        return a;
    }
    std::vector<Point> result;
    result.reserve(a.size() + b.size());
    std::size_t i = 0, j = 0;
    while (i < a.size() || j < b.size()) {
        result.push_back(a[i % a.size()] + b[j % b.size()]);
        if (i == a.size()) { ++j; continue; }
        if (j == b.size()) { ++i; continue; }
        Point u = a[(i+1) % a.size()] - a[i];
        Point v = b[(j+1) % b.size()] - b[j];
        if (angle_less(u, v)) ++i;
        else if (angle_less(v, u)) ++j;
        else { ++i; ++j; }
    }
    return algorithms_detail::strict_boundary(result);
}

// Arbitrary point clouds are first replaced by their convex hulls.
// O(N log N + M log M); this is the sum of convex hulls, not a nonconvex union.
inline std::vector<Point> minkowski_sum(std::vector<Point> a, std::vector<Point> b) {
    return minkowski_sum_convex(convex_hull(std::move(a)), convex_hull(std::move(b)));
}

// Area centroid of a simple polygon, in either orientation. O(N).
// Throws for fewer than three points or absolute doubled area <= eps.
inline Point polygon_centroid(const std::vector<Point>& polygon) {
    if (polygon.size() < 3) throw std::invalid_argument("polygon_centroid needs positive area");
    Point origin = polygon[0], moment{};
    Real twice_area = 0;
    for (std::size_t i = 1; i + 1 < polygon.size(); ++i) {
        Point a = polygon[i] - origin, b = polygon[i+1] - origin;
        Real triangle_area = cross(a, b);
        twice_area += triangle_area;
        moment = moment + (a + b) * triangle_area;
    }
    if (std::abs(twice_area) <= eps) throw std::invalid_argument("polygon_centroid needs positive area");
    return origin + moment / (3 * twice_area);
}

// Closed segments, including point segments. O(1).
inline Real segment_segment_distance(Point a, Point b, Point c, Point d) {
    if (segments_intersect(a, b, c, d)) return 0;
    return (std::min)({segment_distance(a, b, c), segment_distance(a, b, d),
                       segment_distance(c, d, a), segment_distance(c, d, b)});
}

// Reflection of a point across an infinite line. O(1).
inline Point reflection(const Line& line, Point point) {
    return projection(line, point) * 2 - point;
}

} // namespace cp_geometry
#include "detail/macro_pop.hpp"
