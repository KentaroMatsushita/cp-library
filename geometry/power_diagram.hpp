#pragma once
#include "detail/macro_push.hpp"
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include "halfplane.hpp"
#include "circle.hpp"

namespace cp_geometry {

// The power distance to a site is |query-point|^2-weight.
// Larger weights enlarge the cell. Negative weights are allowed.
struct PowerSite {
    Point point;
    Real weight = 0;
};

namespace power_detail {
inline void validate_site(const PowerSite& site) {
    if (!std::isfinite(site.point.x) || !std::isfinite(site.point.y) || !std::isfinite(site.weight))
        throw std::invalid_argument("power sites need finite coordinates and weights");
}
inline void validate_sites(const std::vector<PowerSite>& sites) {
    for (const auto& site : sites) validate_site(site);
}
inline void validate_bounds(Point lower, Point upper) {
    if (!std::isfinite(lower.x) || !std::isfinite(lower.y) ||
        !std::isfinite(upper.x) || !std::isfinite(upper.y) ||
        !(lower.x < upper.x && lower.y < upper.y))
        throw std::invalid_argument("power diagram needs a finite positive bounding rectangle");
}
inline Point translated(Point point, Point origin) {
    const Point result = point-origin;
    if (!std::isfinite(result.x) || !std::isfinite(result.y))
        throw std::overflow_error("translated power coordinates are out of range");
    return result;
}

// The query is expressed relative to first.point. No large squared coordinates
// are subtracted. Unit normals also avoid squaring the center separation.
inline HalfPlane local_bisector(const PowerSite& first, const PowerSite& second) {
    const Point delta = translated(second.point, first.point);
    const Real separation = std::hypot(delta.x, delta.y);
    if (separation == 0) return {0, 0, first.weight >= second.weight ? Real(0) : Real(-1)};
    if (!std::isfinite(separation))
        throw std::overflow_error("power site separation is out of range");
    const Point unit = delta/separation;
    const Real offset = separation/2+(first.weight/2-second.weight/2)/separation;
    if (!std::isfinite(offset))
        throw std::overflow_error("power bisector is out of range");
    return {unit.x, unit.y, offset};
}

inline std::vector<HalfPlane> constraints(const std::vector<PowerSite>& sites, std::size_t index) {
    if (index >= sites.size()) throw std::out_of_range("power cell site index");
    std::vector<HalfPlane> result;
    result.reserve(sites.size()-1);
    for (std::size_t other = 0; other < sites.size(); ++other) {
        if (other == index) continue;
        const auto& selected = sites[index];
        const auto& competing = sites[other];
        if (selected.point.x == competing.point.x && selected.point.y == competing.point.y &&
            selected.weight == competing.weight) {
            // Identical sites would otherwise share a positive-area cell.
            // The smallest original input index owns that cell.
            if (other < index) return {{0, 0, -1}};
            continue;
        }
        result.push_back(local_bisector(selected, competing));
    }
    return result;
}

inline std::vector<Point> bounded_cell(const std::vector<PowerSite>& sites, std::size_t index,
                                       Point lower, Point upper) {
    const auto halfplanes = constraints(sites, index);
    const Point origin = sites[index].point;
    auto polygon = half_plane_intersection(halfplanes, translated(lower, origin), translated(upper, origin));
    for (auto& point : polygon) point = point+origin;
    return polygon;
}

inline Real polygon_area(const std::vector<Point>& polygon) {
    if (polygon.size() < 3) return 0;
    const Point origin = polygon.front();
    Real twice_area = 0;
    for (std::size_t index = 1; index+1 < polygon.size(); ++index)
        twice_area += cross(polygon[index]-origin, polygon[index+1]-origin);
    return std::abs(twice_area)/2;
}
} // namespace power_detail

inline Real power_value(const PowerSite& site, Point query) {
    power_detail::validate_site(site);
    if (!std::isfinite(query.x) || !std::isfinite(query.y))
        throw std::invalid_argument("power query must be finite");
    const Real result = norm2(power_detail::translated(query, site.point))-site.weight;
    if (!std::isfinite(result)) throw std::overflow_error("power value is out of range");
    return result;
}

// Closed half-plane where first has no greater power distance than second.
// Same-position sites yield the entire plane or the impossible half-plane.
// The individual function keeps equality; index tie-breaking is applied by cells.
inline HalfPlane power_bisector(const PowerSite& first, const PowerSite& second) {
    power_detail::validate_site(first);
    power_detail::validate_site(second);
    HalfPlane result = power_detail::local_bisector(first, second);
    result.c += result.a*first.point.x+result.b*first.point.y;
    if (!std::isfinite(result.c)) throw std::overflow_error("power bisector is out of range");
    return result;
}

// Intersect the selected cell with the explicit rectangle. O(N log N), O(N) memory.
// Empty, point, and segment cells return {}. Positive-area results are CCW.
inline std::vector<Point> power_cell(const std::vector<PowerSite>& sites, std::size_t index,
                                    Point lower, Point upper) {
    power_detail::validate_sites(sites);
    power_detail::validate_bounds(lower, upper);
    return power_detail::bounded_cell(sites, index, lower, upper);
}

// No artificial rectangle. Expected O(N log N), worst O(N^2), inherited from HPI.
// DEGENERATE reports one witness for a bounded point/segment, not its endpoints.
inline HalfPlaneRegion power_cell(const std::vector<PowerSite>& sites, std::size_t index) {
    power_detail::validate_sites(sites);
    auto region = half_plane_intersection(power_detail::constraints(sites, index));
    if (region.kind != RegionKind::EMPTY) {
        const Point origin = sites[index].point;
        region.witness = region.witness+origin;
        for (auto& point : region.polygon) point = point+origin;
    }
    return region;
}

// All cells in the same bounding rectangle, in original site order.
// O(N^2 log N) time, O(N^2) total output memory in this independent-cell method.
// This does not construct the adjacency graph of a power diagram.
inline std::vector<std::vector<Point>> power_diagram(const std::vector<PowerSite>& sites,
                                                    Point lower, Point upper) {
    power_detail::validate_sites(sites);
    power_detail::validate_bounds(lower, upper);
    std::vector<std::vector<Point>> result;
    result.reserve(sites.size());
    for (std::size_t index = 0; index < sites.size(); ++index)
        result.push_back(power_detail::bounded_cell(sites, index, lower, upper));
    return result;
}

// The power of a circle is |query-center|^2-radius^2.
inline std::vector<PowerSite> power_sites_from_circles(const std::vector<Circle>& circles) {
    std::vector<PowerSite> sites;
    sites.reserve(circles.size());
    for (const auto& circle : circles) {
        if (!std::isfinite(circle.radius) || circle.radius < 0)
            throw std::invalid_argument("circle power needs a finite nonnegative radius");
        PowerSite site{circle.center, circle.radius*circle.radius};
        power_detail::validate_site(site);
        sites.push_back(site);
    }
    return sites;
}

// Area assigned to this circle by minimum power, outside every input disk.
// Inside its power cell, being outside its OWN disk implies being outside all
// disks, so only one circle-polygon intersection is subtracted. O(N log N).
inline Real circle_territory_area(const std::vector<Circle>& circles, std::size_t index,
                                   Point lower, Point upper) {
    const auto polygon = power_cell(power_sites_from_circles(circles), index, lower, upper);
    if (polygon.empty()) return 0;
    const Real area = power_detail::polygon_area(polygon);
    const Real excluded = circle_polygon_intersection_area(circles[index], polygon);
    return std::fmax(Real(0), area-excluded);
}

} // namespace cp_geometry
#include "detail/macro_pop.hpp"
