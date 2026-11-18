#pragma once
#include "detail/macro_push.hpp"
#include "convex.hpp"
#include <deque>
#include <limits>
#include <random>

namespace cp_geometry {
namespace detail {
// The angular order is cyclic: nearly equal directions can lie on opposite
// sides of angle zero. Keep the stricter line without breaking that order.
inline void merge_cyclic_halfplanes(std::vector<Line>& lines) {
    std::size_t first=0;
    while(lines.size()-first>1 &&
          sign(cross(lines[first].direction(),lines.back().direction()))==0 &&
          dot(lines[first].direction(),lines.back().direction())>0) {
        if(cross(lines[first].direction(),lines.back().a-lines[first].a)>0) ++first;
        else lines.pop_back();
    }
    lines.erase(lines.begin(),lines.begin()+first);
}
} // namespace detail
// Intersect left half-planes inside an EXPLICIT finite bounding rectangle.
// Returns only positive-area intersections; empty/point/segment => {}.
// O(N log N), including sorting. The rectangle is part of the constraints.
inline std::vector<Point> half_plane_intersection(std::vector<Line> lines,
                                                Point lower,Point upper) {
    if(!(lower.x<upper.x && lower.y<upper.y)) throw std::invalid_argument("invalid bounding box");
    Point a=lower,b{upper.x,lower.y},c=upper,d{lower.x,upper.y};
    lines.insert(lines.end(),{Line(a,b),Line(b,c),Line(c,d),Line(d,a)});
    for(auto& line:lines) {
        const Point direction=line.direction()/length(line.direction());
        line.b=line.a+direction;
    }
    std::sort(lines.begin(),lines.end(),[](const Line& l,const Line& r){return angle_less(l.direction(),r.direction());});
    std::vector<Line> unique;
    for(const auto& l:lines) {
        if(!unique.empty() && sign(cross(unique.back().direction(),l.direction()))==0 && dot(unique.back().direction(),l.direction())>0) {
            if(cross(unique.back().direction(),l.a-unique.back().a)>0) unique.back()=l;
        } else unique.push_back(l);
    }
    detail::merge_cyclic_halfplanes(unique);
    std::deque<Line> q;
    auto outside=[](const Line& l,Point p){return sign(cross(l.direction(),p-l.a))<0;};
    auto meet=[](const Line& l,const Line& r){return line_intersection(l,r);};
    for(const auto& l:unique) {
        while(q.size()>1 && outside(l,meet(q[q.size()-2],q.back()))) q.pop_back();
        while(q.size()>1 && outside(l,meet(q[0],q[1]))) q.pop_front();
        if(!q.empty() && sign(cross(q.back().direction(),l.direction()))==0) return {};
        q.push_back(l);
    }
    while(q.size()>2 && outside(q.front(),meet(q[q.size()-2],q.back()))) q.pop_back();
    while(q.size()>2 && outside(q.back(),meet(q[0],q[1]))) q.pop_front();
    if(q.size()<3 || sign(cross(q.front().direction(),q.back().direction()))==0) return {};
    std::vector<Point> result;
    for(std::size_t i=0;i<q.size();++i) result.push_back(meet(q[i],q[(i+1)%q.size()]));
    // eps is used for distance/angle predicates, not as an area cutoff.
    if(!(area2(result)>0)) return {};
    return result;
}
enum class RegionKind { EMPTY, UNBOUNDED, POLYGON, DEGENERATE };
struct HalfPlaneRegion {
    RegionKind kind;
    Point witness{}; // A feasible point for every nonempty result.
    std::vector<Point> polygon; // Only populated for POLYGON.
};
// General half-plane intersection, without an artificial bounding box.
// Randomized incremental minimum-norm feasibility: expected O(N), worst O(N^2).
// Sorting/deque: O(N log N). Total expected O(N log N), worst O(N^2).
// DEGENERATE reports a bounded point/segment via a feasible witness, not endpoints.
inline HalfPlaneRegion half_plane_intersection(std::vector<Line> lines) {
    for(auto& l:lines) l.b=l.a+l.direction()/length(l.direction());
    std::mt19937 random_engine(712367);
    std::shuffle(lines.begin(),lines.end(),random_engine);
    Point feasible{};
    for(std::size_t i=0;i<lines.size();++i) {
        const auto& l=lines[i];
        if(sign(cross(l.direction(),feasible-l.a))>=0) continue;
        Real left=-std::numeric_limits<Real>::infinity(),right=-left;
        Point v=l.direction();
        for(std::size_t j=0;j<i;++j) {
            Real a=cross(lines[j].direction(),v);
            Real b=cross(lines[j].direction(),l.a-lines[j].a);
            if(sign(a)==0) {if(sign(b)<0) return {RegionKind::EMPTY,{},{}};}
            else if(a>0) left=std::max(left,-b/a);
            else right=std::min(right,-b/a);
        }
        if(left>right+eps) return {RegionKind::EMPTY,{},{}};
        Real t=std::max(left,std::min(right,-dot(l.a,v)));
        feasible=l.a+v*t;
    }
    if(lines.size()<3) return {RegionKind::UNBOUNDED,feasible,{}};
    std::sort(lines.begin(),lines.end(),[](const Line& a,const Line& b){return angle_less(a.direction(),b.direction());});
    // A nonzero recession direction exists iff constraint directions fit in a semicircle.
    const Real pi=std::acos(Real(-1));
    for(std::size_t i=0;i<lines.size();++i) {
        Point a=lines[i].direction(),b=lines[(i+1)%lines.size()].direction();
        Real aa=std::atan2(a.y,a.x),bb=std::atan2(b.y,b.x);
        if(aa<0) aa+=2*pi;
        if(bb<0) bb+=2*pi;
        Real gap=bb-aa;
        if(i+1==lines.size()) gap+=2*pi;
        if(gap>=pi-eps) return {RegionKind::UNBOUNDED,feasible,{}};
    }
    std::vector<Line> unique;
    for(const auto& l:lines) {
        if(!unique.empty() && sign(cross(unique.back().direction(),l.direction()))==0 && dot(unique.back().direction(),l.direction())>0) {
            if(cross(unique.back().direction(),l.a-unique.back().a)>0) unique.back()=l;
        } else unique.push_back(l);
    }
    detail::merge_cyclic_halfplanes(unique);
    std::deque<Line> q;
    auto outside=[](const Line& l,Point p){return sign(cross(l.direction(),p-l.a))<0;};
    auto meet=[](const Line& l,const Line& r){return line_intersection(l,r);};
    for(const auto& l:unique) {
        while(q.size()>1 && outside(l,meet(q[q.size()-2],q.back()))) q.pop_back();
        while(q.size()>1 && outside(l,meet(q[0],q[1]))) q.pop_front();
        if(!q.empty() && sign(cross(q.back().direction(),l.direction()))==0)
            return {RegionKind::DEGENERATE,feasible,{}};
        q.push_back(l);
    }
    while(q.size()>2 && outside(q.front(),meet(q[q.size()-2],q.back()))) q.pop_back();
    while(q.size()>2 && outside(q.back(),meet(q[0],q[1]))) q.pop_front();
    if(q.size()<3 || sign(cross(q.front().direction(),q.back().direction()))==0)
        return {RegionKind::DEGENERATE,feasible,{}};
    std::vector<Point> result;
    for(std::size_t i=0;i<q.size();++i) result.push_back(meet(q[i],q[(i+1)%q.size()]));
    if(!(area2(result)>0)) return {RegionKind::DEGENERATE,feasible,{}};
    return {RegionKind::POLYGON,feasible,std::move(result)};
}

// Closed linear inequality a*x+b*y <= c. Zero normals represent either the
// entire plane (c>=0) or an impossible constraint (c<0).
struct HalfPlane { Real a=0,b=0,c=0; };
namespace detail {
inline bool coefficient_lines(const std::vector<HalfPlane>& constraints,std::vector<Line>& lines) {
    lines.reserve(constraints.size());
    bool feasible=true;
    for(const HalfPlane& h:constraints) {
        if(!std::isfinite(h.a)||!std::isfinite(h.b)||!std::isfinite(h.c))
            throw std::invalid_argument("half-plane coefficients must be finite");
        const Real scale=std::max(std::abs(h.a),std::abs(h.b));
        if(scale==0) { if(h.c<0) feasible=false; continue; }
        const Point normal{h.a/scale,h.b/scale};
        const Real norm=length(normal);
        const Point unit=normal/norm;
        const Real offset=(h.c/scale)/norm;
        if(!std::isfinite(offset)) throw std::overflow_error("half-plane boundary is out of range");
        const Point anchor=unit*offset;
        // cross(perpendicular(unit), p-anchor) = offset-dot(unit,p).
        lines.emplace_back(anchor,anchor+perpendicular(unit));
    }
    return feasible;
}
} // namespace detail
inline HalfPlaneRegion half_plane_intersection(const std::vector<HalfPlane>& constraints) {
    std::vector<Line> lines;
    if(!detail::coefficient_lines(constraints,lines)) return {RegionKind::EMPTY,{},{}};
    return half_plane_intersection(std::move(lines));
}
inline std::vector<Point> half_plane_intersection(const std::vector<HalfPlane>& constraints,
                                                Point lower,Point upper) {
    if(!(lower.x<upper.x && lower.y<upper.y)) throw std::invalid_argument("invalid bounding box");
    std::vector<Line> lines;
    if(!detail::coefficient_lines(constraints,lines)) return {};
    return half_plane_intersection(std::move(lines),lower,upper);
}

} // namespace cp_geometry
#include "detail/macro_pop.hpp"
