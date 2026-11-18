#pragma once
#include "detail/macro_push.hpp"
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("no-fast-math")
#endif
#include "base.hpp"
#include "delaunay.hpp"
#include "halfplane.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cp_geometry {

enum class VoronoiEdgeKind { segment, ray, line };

// origin + direction*t, with t in [0,1], [0,+inf), or (-inf,+inf).
// Site indices always refer to the ORIGINAL input, using its duplicate owner.
struct VoronoiEdge {
    std::size_t first_site, second_site;
    VoronoiEdgeKind kind;
    Point origin, direction;

    bool contains_parameter(Real parameter) const {
        if(kind==VoronoiEdgeKind::line) return true;
        if(parameter < -eps) return false;
        return kind==VoronoiEdgeKind::ray || parameter<=1+eps;
    }

    // Clip any kind of edge to a closed finite rectangle. A point contact is
    // returned twice. No intersection => nullopt. O(1).
    std::optional<std::pair<Point,Point>> clip(Point lower,Point upper) const {
        if(!std::isfinite(lower.x)||!std::isfinite(lower.y)||
           !std::isfinite(upper.x)||!std::isfinite(upper.y)||
           !(lower.x<upper.x && lower.y<upper.y))
            throw std::invalid_argument("invalid Voronoi clipping box");
        Real left=kind==VoronoiEdgeKind::line ? -std::numeric_limits<Real>::infinity() : 0;
        Real right=kind==VoronoiEdgeKind::segment ? 1 : std::numeric_limits<Real>::infinity();
        auto restrict_axis=[&](Real coordinate,Real speed,Real minimum,Real maximum) {
            if(speed==0) return coordinate>=minimum-eps && coordinate<=maximum+eps;
            Real a=(minimum-coordinate)/speed,b=(maximum-coordinate)/speed;
            if(a>b) std::swap(a,b);
            left=std::max(left,a); right=std::min(right,b);
            return left<=right;
        };
        if(!restrict_axis(origin.x,direction.x,lower.x,upper.x) ||
           !restrict_axis(origin.y,direction.y,lower.y,upper.y)) return std::nullopt;
        return std::pair<Point,Point>{origin+direction*left,origin+direction*right};
    }

    // Intersection with a CLOSED segment. A coinciding overlap returns its
    // endpoints, not an infinite list. O(1), ordered from first toward second.
    std::vector<Point> segment_intersections(Point first,Point second) const {
        const Real squared_length=norm2(direction);
        if(squared_length==0) return {};
        const Point side=second-first,delta=first-origin;
        if(first.x==second.x && first.y==second.y) {
            Real parameter=dot(delta,direction)/squared_length;
            if(contains_parameter(parameter) &&
               std::abs(cross(direction/length(direction),delta))<=eps) return {first};
            return {};
        }
        const Real determinant=cross(direction,side);
        if(determinant!=0) {
            const Real t=cross(delta,side)/determinant;
            const Real u=cross(delta,direction)/determinant;
            if(!contains_parameter(t)||u < -eps||u>1+eps) return {};
            return {first+side*std::clamp(u,Real(0),Real(1))};
        }
        if(std::abs(cross(direction/length(direction),delta))>eps) return {};
        Real a=dot(first-origin,direction)/squared_length;
        Real b=dot(second-origin,direction)/squared_length;
        Real left=std::min(a,b),right=std::max(a,b);
        if(kind!=VoronoiEdgeKind::line) left=std::max(left,Real(0));
        if(kind==VoronoiEdgeKind::segment) right=std::min(right,Real(1));
        if(left>right) return {};
        Point p=origin+direction*left,q=origin+direction*right;
        if(a>b) std::swap(p,q);
        if(near(p,q)) return {p};
        return {p,q};
    }
};

struct VoronoiVertex {
    Point point;
    std::size_t site; // One incident site, in the original input.
};

namespace voronoi_detail {
inline void validate_box(Point lower,Point upper) {
    if(!std::isfinite(lower.x)||!std::isfinite(lower.y)||
       !std::isfinite(upper.x)||!std::isfinite(upper.y)||
       !(lower.x<upper.x && lower.y<upper.y))
        throw std::invalid_argument("invalid Voronoi clipping box");
}
inline Point triangle_center(Point a,Point b,Point c) {
    const Point u=b-a,v=c-a;
    const Real positive=u.x*v.y,negative=u.y*v.x;
    const Real determinant=positive-negative;
    Point center;
    if(std::isfinite(determinant) && std::abs(determinant)>
       64*std::numeric_limits<Real>::epsilon()*(std::abs(positive)+std::abs(negative))) {
        center=a+(perpendicular(u)*norm2(v)-perpendicular(v)*norm2(u))/(2*determinant);
    } else {
        // A certified noncollinear triangle can still have cross(u,v)==0 after
        // rounded products cancel. Evaluate the translated circumcenter's
        // determinant AND numerator with error-free expansions in scaled units.
        const Real magnitude=(std::max)({std::abs(a.x),std::abs(a.y),std::abs(b.x),
                                        std::abs(b.y),std::abs(c.x),std::abs(c.y)});
        int exponent=0;
        if(magnitude!=0) std::frexp(magnitude,&exponent);
        auto scaled=[&](Point p) {return Point(std::ldexp(p.x,-exponent),std::ldexp(p.y,-exponent));};
        const Point sa=scaled(a),sb=scaled(b),sc=scaled(c);
        const auto ux=delaunay_detail::difference(sb.x,sa.x),uy=delaunay_detail::difference(sb.y,sa.y);
        const auto vx=delaunay_detail::difference(sc.x,sa.x),vy=delaunay_detail::difference(sc.y,sa.y);
        const auto uu=delaunay_detail::add(delaunay_detail::multiply(ux,ux),delaunay_detail::multiply(uy,uy));
        const auto vv=delaunay_detail::add(delaunay_detail::multiply(vx,vx),delaunay_detail::multiply(vy,vy));
        const auto den=delaunay_detail::determinant(ux,uy,vx,vy);
        const auto x=delaunay_detail::add(delaunay_detail::multiply(uu,vy),delaunay_detail::multiply(vv,uy),-1);
        const auto y=delaunay_detail::add(delaunay_detail::multiply(vv,ux),delaunay_detail::multiply(uu,vx),-1);
        auto estimate=[](const delaunay_detail::Expansion& expansion) {
            Real result=0; for(Real part:expansion) result+=part; return result;
        };
        const Real denominator=2*estimate(den);
        if(denominator==0) throw std::runtime_error("degenerate Delaunay triangle");
        center={std::ldexp(sa.x+estimate(x)/denominator,exponent),
                std::ldexp(sa.y+estimate(y)/denominator,exponent)};
    }
    if(!std::isfinite(center.x)||!std::isfinite(center.y))
        throw std::overflow_error("Voronoi vertex out of range");
    return center;
}
} // namespace voronoi_detail

// Ordinary Euclidean Voronoi diagram via the Delaunay dual. O(N log N)
// construction and O(N) storage. See delaunay.hpp for predicate precision.
// Exact duplicate sites are owned by their smallest ORIGINAL index.
class VoronoiDiagram {
    std::vector<Point> input_;
    DelaunayTriangulation triangulation_;
    std::vector<std::vector<std::size_t>> neighbors_;
    std::vector<VoronoiEdge> edges_;
    std::vector<VoronoiVertex> vertices_;

    void build_dual() {
        const auto& sites=triangulation_.sites;
        const auto& triangles=triangulation_.triangles;
        neighbors_.resize(sites.size());
        using Edge=std::pair<std::size_t,std::size_t>;
        std::map<Edge,std::vector<std::pair<std::size_t,std::size_t>>> faces;
        std::vector<Point> centers;
        centers.reserve(triangles.size());
        for(std::size_t t=0;t<triangles.size();++t) {
            const auto& triangle=triangles[t];
            centers.push_back(voronoi_detail::triangle_center(sites[triangle[0]],sites[triangle[1]],sites[triangle[2]]));
            for(int k=0;k<3;++k) {
                auto edge=std::minmax(triangle[k],triangle[(k+1)%3]);
                faces[edge].push_back({t,triangle[(k+2)%3]});
            }
        }
        std::vector<std::size_t> parent(triangles.size());
        std::iota(parent.begin(),parent.end(),0);
        auto root=[&](std::size_t v) {
            while(parent[v]!=v) { parent[v]=parent[parent[v]]; v=parent[v]; }
            return v;
        };
        // A triangulation may choose arbitrary diagonals on cocircular sites.
        // Merge their coincident circumcenters and omit the zero-length duals.
        for(const auto& item:faces) if(item.second.size()==2) {
            const auto a=item.second[0].first,b=item.second[1].first;
            if(near(centers[a],centers[b])) parent[root(a)]=root(b);
        }
        for(std::size_t i=0;i<parent.size();++i) parent[i]=root(i);
        for(std::size_t i=0;i<triangles.size();++i) if(parent[i]==i)
            vertices_.push_back({centers[i],triangulation_.representative[triangles[i][0]]});

        for(const auto& pair:triangulation_.edges) {
            const auto u=pair.first,v=pair.second;
            neighbors_[u].push_back(v); neighbors_[v].push_back(u);
            const Point delta=sites[v]-sites[u];
            const Point normal=perpendicular(delta)/length(delta);
            const auto found=faces.find(pair);
            const auto original_u=triangulation_.representative[u];
            const auto original_v=triangulation_.representative[v];
            if(found==faces.end()) {
                edges_.push_back({original_u,original_v,VoronoiEdgeKind::line,
                                  sites[u]+delta/2,normal});
            } else if(found->second.size()==1) {
                const auto face=found->second[0];
                // Use certified CCW topology instead of recomputing a rounded
                // orientation, which may be zero for nearly collinear sites.
                bool forward=false;
                const auto& triangle=triangles[face.first];
                for(int k=0;k<3;++k)
                    if(triangle[k]==u && triangle[(k+1)%3]==v) forward=true;
                const Point outward=forward ? normal*-1 : normal;
                edges_.push_back({original_u,original_v,VoronoiEdgeKind::ray,
                                  centers[parent[face.first]],outward});
            } else {
                const auto a=parent[found->second[0].first],b=parent[found->second[1].first];
                if(a==b) continue;
                edges_.push_back({original_u,original_v,VoronoiEdgeKind::segment,
                                  centers[a],centers[b]-centers[a]});
            }
        }
    }

public:
    explicit VoronoiDiagram(std::vector<Point> sites)
        : input_(std::move(sites)),triangulation_(delaunay_triangulation(input_)) { build_dual(); }

    const std::vector<Point>& sites() const { return input_; }
    const DelaunayTriangulation& triangulation() const { return triangulation_; }
    const std::vector<VoronoiEdge>& edges() const { return edges_; }
    const std::vector<VoronoiVertex>& vertices() const { return vertices_; }
    std::size_t site_owner(std::size_t index) const {
        if(index>=input_.size()) throw std::out_of_range("Voronoi site index");
        return triangulation_.representative[triangulation_.original_to_site[index]];
    }

    // O(d log d) for d Delaunay neighbors. The original index must be used.
    // A duplicate non-owner or zero-area/empty clipped cell returns {}.
    std::vector<Point> cell(std::size_t index,Point lower,Point upper) const {
        voronoi_detail::validate_box(lower,upper);
        if(site_owner(index)!=index) return {};
        const std::size_t unique=triangulation_.original_to_site[index];
        const Point origin=triangulation_.sites[unique];
        std::vector<HalfPlane> constraints;
        constraints.reserve(neighbors_[unique].size());
        for(const auto other:neighbors_[unique]) {
            const Point delta=triangulation_.sites[other]-origin;
            constraints.push_back({delta.x,delta.y,norm2(delta)/2});
        }
        auto polygon=half_plane_intersection(constraints,lower-origin,upper-origin);
        for(auto& point:polygon) point=point+origin;
        return polygon;
    }

    // All assigned cells, indexed like the original input. O(N log N) total,
    // O(N) total output for a RECTANGULAR domain (not N independent all-pair HPI).
    std::vector<std::vector<Point>> cells(Point lower,Point upper) const {
        voronoi_detail::validate_box(lower,upper);
        std::vector<std::vector<Point>> result(input_.size());
        for(const auto original:triangulation_.representative) result[original]=cell(original,lower,upper);
        return result;
    }
};

inline VoronoiDiagram voronoi_diagram(std::vector<Point> sites) {
    return VoronoiDiagram(std::move(sites));
}
inline std::vector<std::vector<Point>> voronoi_cells(std::vector<Point> sites,Point lower,Point upper) {
    return VoronoiDiagram(std::move(sites)).cells(lower,upper);
}
inline std::vector<Point> voronoi_cell(std::vector<Point> sites,std::size_t index,Point lower,Point upper) {
    return VoronoiDiagram(std::move(sites)).cell(index,lower,upper);
}

struct EmptyCircleResult {
    Point center;
    Real squared_radius;
    std::size_t nearest_site;
};

// Largest empty circle whose CENTER lies in the box. The disk itself need not
// fit inside the box. Input points are forbidden in its open interior.
// O(N log N), nullopt if there are no sites. This is NOT a disk-inscribed-in-box API.
inline std::optional<EmptyCircleResult> largest_empty_circle_in_box(const VoronoiDiagram& diagram,
                                                                  Point lower,Point upper) {
    auto cells=diagram.cells(lower,upper);
    std::optional<EmptyCircleResult> best;
    for(std::size_t i=0;i<cells.size();++i) for(const auto vertex:cells[i]) {
        const Real squared_radius=norm2(vertex-diagram.sites()[i]);
        if(!best || squared_radius>best->squared_radius) best=EmptyCircleResult{vertex,squared_radius,i};
    }
    return best;
}

// A SIMPLE polygon, in either orientation, may be concave. Only the CENTER is
// constrained to it, not the entire disk. O(N*M + N log N) with M polygon vertices.
// The diagram's finite vertices, boundary intersections and polygon vertices
// contain a maximizer. Nullopt for no sites or a polygon with fewer than 3 vertices.
inline std::optional<EmptyCircleResult> largest_empty_circle_in_polygon(const VoronoiDiagram& diagram,
                                                                      const std::vector<Point>& polygon) {
    if(diagram.sites().empty() || polygon.size()<3) return std::nullopt;
    std::optional<EmptyCircleResult> best;
    auto consider=[&](Point point,std::size_t site) {
        const Real squared_radius=norm2(point-diagram.sites()[site]);
        if(!best || squared_radius>best->squared_radius) best=EmptyCircleResult{point,squared_radius,site};
    };
    for(const auto& vertex:diagram.vertices())
        if(contains(polygon,vertex.point)!=Location::outside) consider(vertex.point,vertex.site);
    for(const auto& edge:diagram.edges()) for(std::size_t i=0;i<polygon.size();++i)
        for(const auto intersection:edge.segment_intersections(polygon[i],polygon[(i+1)%polygon.size()]))
            consider(intersection,edge.first_site);
    for(const auto point:polygon) {
        std::size_t closest=0;
        for(std::size_t i=1;i<diagram.sites().size();++i)
            if(norm2(point-diagram.sites()[i])<norm2(point-diagram.sites()[closest])) closest=i;
        consider(point,diagram.site_owner(closest));
    }
    return best;
}

} // namespace cp_geometry
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif
#include "detail/macro_pop.hpp"
