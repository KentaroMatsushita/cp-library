#pragma once
#include "detail/macro_push.hpp"
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <vector>
#include <utility>

namespace cp_geometry {
using Real = long double;
inline constexpr Real eps = 1e-12L;
inline int sign(Real x) { return (x > eps) - (x < -eps); }
struct Point {
    Real x = 0, y = 0;
    constexpr Point() = default;
    constexpr Point(Real x, Real y) : x(x), y(y) {}
    Point operator+(Point p) const { return {x+p.x,y+p.y}; }
    Point operator-(Point p) const { return {x-p.x,y-p.y}; }
    Point operator*(Real k) const { return {x*k,y*k}; }
    Point operator/(Real k) const { return {x/k,y/k}; }
};
inline Real dot(Point a, Point b) { return a.x*b.x+a.y*b.y; }
inline Real cross(Point a, Point b) { return a.x*b.y-a.y*b.x; }
inline Real norm2(Point p) { return dot(p,p); }
inline Point perpendicular(Point p) { return {-p.y,p.x}; }
inline bool near(Point a, Point b) { return norm2(a-b)<=eps*eps; }
inline bool on_segment(Point a, Point b, Point p) {
    if (a.x==b.x && a.y==b.y) return near(a,p);
    const Real length=std::hypot(b.x-a.x,b.y-a.y);
    return std::abs(cross((b-a)/length,p-a))<=eps && dot(p-a,p-b)<=eps*length;
}
struct Line {
    Point a,b;
    Line(Point a, Point b):a(a),b(b) {
        if(!std::isfinite(a.x)||!std::isfinite(a.y)||!std::isfinite(b.x)||!std::isfinite(b.y)||near(a,b))
            throw std::invalid_argument("line needs finite distinct points");
    }
    Point direction() const { return b-a; }
};
// Requires nonparallel lines.
inline Point line_intersection(const Line& l, const Line& m) {
    const Point u=l.direction()/std::hypot(l.direction().x,l.direction().y);
    const Point v=m.direction()/std::hypot(m.direction().x,m.direction().y);
    const Real d=cross(u,v);
    if(sign(d)==0) throw std::invalid_argument("parallel lines");
    return l.a+u*(cross(m.a-l.a,v)/d);
}
inline Point projection(const Line& l, Point p) {
    Point v=l.direction(); return l.a+v*(dot(p-l.a,v)/norm2(v));
}
inline Real segment_distance(Point a, Point b, Point p) {
    if(near(a,b)) return std::sqrt(norm2(p-a));
    Real t=std::clamp(dot(p-a,b-a)/norm2(b-a),Real(0),Real(1));
    return std::sqrt(norm2(p-(a+(b-a)*t)));
}
inline bool segments_intersect(Point a, Point b, Point c, Point d) {
    if(on_segment(a,b,c)||on_segment(a,b,d)||on_segment(c,d,a)||on_segment(c,d,b)) return true;
    return sign(cross(b-a,c-a))*sign(cross(b-a,d-a))<0 &&
           sign(cross(d-c,a-c))*sign(cross(d-c,b-c))<0;
}
inline Real area2(const std::vector<Point>& p) {
    // Translate before multiplying: global-coordinate shoelace terms may cancel
    // completely for a small polygon far from the origin. Kahan accumulation
    // also reduces cancellation between positive/negative concave fan triangles.
    Real s=0,correction=0;
    for(std::size_t i=1;i+1<p.size();++i) {
        const Real contribution=cross(p[i]-p[0],p[i+1]-p[0])-correction;
        const Real updated=s+contribution;
        correction=(updated-s)-contribution;
        s=updated;
    }
    return s;
}
enum class Location { outside=0, boundary=1, inside=2 };
inline Location contains(const std::vector<Point>& poly,Point p) {
    bool in=false;
    for(std::size_t i=0;i<poly.size();++i) {
        Point a=poly[i],b=poly[(i+1)%poly.size()];
        if(on_segment(a,b,p)) return Location::boundary;
        if(a.y>b.y) std::swap(a,b);
        if(a.y<=p.y && p.y<b.y && cross(b-a,p-a)>0) in=!in;
    }
    return in?Location::inside:Location::outside;
}
inline Real length(Point p) { return std::hypot(p.x,p.y); }
inline Real distance(Point a,Point b) { return length(a-b); }
inline Point rotate(Point p,Real radians) {
    const Real c=std::cos(radians),s=std::sin(radians);
    return {p.x*c-p.y*s,p.x*s+p.y*c};
}
inline Real polygon_area(const std::vector<Point>& polygon) { return std::abs(area2(polygon))/2; }

} // namespace cp_geometry
#include "detail/macro_pop.hpp"
