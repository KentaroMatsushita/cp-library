#include "geometry/voronoi.hpp"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <limits>
#include <random>
#include <vector>

namespace review {
using cp_geometry::Point;
using cp_geometry::Real;
using cp_geometry::VoronoiDiagram;
using cp_geometry::VoronoiEdge;
using cp_geometry::VoronoiEdgeKind;
using cp_geometry::norm2;
using cp_geometry::dot;
using cp_geometry::cross;
using cp_geometry::polygon_area;

bool close(Real a, Real b, Real tolerance=1e-8L) {
    return std::abs(a-b)<=tolerance*(1+std::abs(a)+std::abs(b));
}
bool close_point(Point a, Point b, Real tolerance=1e-8L) {
    return std::hypot(a.x-b.x,a.y-b.y)<=tolerance*(1+std::hypot(a.x,a.y)+std::hypot(b.x,b.y));
}
bool same_point(Point a, Point b) { return a.x==b.x && a.y==b.y; }

// Independent oracle: sequential polygon clipping, with no HPI/Delaunay calls.
std::vector<Point> reference_cell(const std::vector<Point>& sites, std::size_t owner, Point lower, Point upper) {
    std::vector<Point> result{lower,{upper.x,lower.y},upper,{lower.x,upper.y}};
    for(std::size_t other=0;other<sites.size();++other) {
        if(other==owner) continue;
        if(same_point(sites[owner],sites[other])) {
            if(other<owner) return {};
            continue;
        }
        const Point normal=sites[other]-sites[owner];
        auto value=[&](Point point) { return dot(point-sites[owner],normal)-norm2(normal)/2; };
        std::vector<Point> next;
        for(std::size_t i=0;i<result.size();++i) {
            const Point a=result[i],b=result[(i+1)%result.size()];
            const Real fa=value(a),fb=value(b);
            if(fa<=0) next.push_back(a);
            if((fa<0 && fb>0)||(fa>0 && fb<0)) next.push_back(a+(b-a)*(fa/(fa-fb)));
        }
        result.swap(next);
        if(result.empty()) break;
    }
    return result;
}

void check_nearest(const std::vector<Point>& sites,std::size_t owner,Point point) {
    assert(owner<sites.size());
    const Real chosen=norm2(point-sites[owner]);
    for(Point site:sites) {
        const Real other=norm2(point-site);
        assert(chosen<=other+1e-7L*(1+chosen+other));
    }
}

void check_diagram(const std::vector<Point>& sites,Point lower={-8,-6},Point upper={9,7}) {
    const VoronoiDiagram diagram(sites);
    const auto cells=diagram.cells(lower,upper);
    assert(cells.size()==sites.size());
    Real total=0;
    for(std::size_t owner=0;owner<sites.size();++owner) {
        std::size_t expected_owner=owner;
        for(std::size_t earlier=0;earlier<owner;++earlier)
            if(same_point(sites[owner],sites[earlier])) { expected_owner=earlier; break; }
        assert(diagram.site_owner(owner)==expected_owner);
        const auto reference=reference_cell(sites,owner,lower,upper);
        const Real area=polygon_area(cells[owner]);
        assert(close(area,polygon_area(reference)));
        total+=area;
        if(owner!=expected_owner) assert(cells[owner].empty());
        for(Point vertex:cells[owner]) {
            assert(vertex.x>=lower.x-1e-7L && vertex.x<=upper.x+1e-7L);
            assert(vertex.y>=lower.y-1e-7L && vertex.y<=upper.y+1e-7L);
            check_nearest(sites,owner,vertex);
        }
    }
    if(!sites.empty()) assert(close(total,(upper.x-lower.x)*(upper.y-lower.y)));
    for(const auto& vertex:diagram.vertices()) check_nearest(sites,vertex.site,vertex.point);
    for(const auto& edge:diagram.edges()) {
        assert(edge.first_site<sites.size() && edge.second_site<sites.size());
        assert(edge.first_site==diagram.site_owner(edge.first_site));
        assert(edge.second_site==diagram.site_owner(edge.second_site));
        assert(edge.first_site!=edge.second_site);
        assert(norm2(edge.direction)>0);
        const std::vector<Real> parameters=edge.kind==VoronoiEdgeKind::segment ?
            std::vector<Real>{0,0.2L,0.5L,1} : edge.kind==VoronoiEdgeKind::ray ?
            std::vector<Real>{0,0.1L,1,10,100} : std::vector<Real>{-100,-1,0,1,100};
        for(Real parameter:parameters) {
            const Point query=edge.origin+edge.direction*parameter;
            assert(close(norm2(query-sites[edge.first_site]),norm2(query-sites[edge.second_site])));
            check_nearest(sites,edge.first_site,query);
        }
        const auto clipped=edge.clip(lower,upper);
        if(clipped) {
            for(Point endpoint:{clipped->first,clipped->second}) {
                assert(endpoint.x>=lower.x-1e-7L && endpoint.x<=upper.x+1e-7L);
                assert(endpoint.y>=lower.y-1e-7L && endpoint.y<=upper.y+1e-7L);
                check_nearest(sites,edge.first_site,endpoint);
            }
        }
    }
    const auto largest=cp_geometry::largest_empty_circle_in_box(diagram,lower,upper);
    assert(bool(largest)==!sites.empty());
    if(largest) {
        Real expected=-1;
        for(std::size_t owner=0;owner<sites.size();++owner)
            for(Point point:reference_cell(sites,owner,lower,upper))
                expected=std::max(expected,norm2(point-sites[owner]));
        assert(close(largest->squared_radius,expected));
        check_nearest(sites,largest->nearest_site,largest->center);
    }
}

void check_edge_primitives() {
    const Point lower{-1,-1},upper{1,1};
    const VoronoiEdge line{0,1,VoronoiEdgeKind::line,{0,0},{1,0}};
    const auto clipped=line.clip(lower,upper);
    assert(clipped && close_point(clipped->first,{-1,0}) && close_point(clipped->second,{1,0}));
    auto points=line.segment_intersections({2,0},{-2,0});
    assert(points.size()==2 && close_point(points[0],{2,0}) && close_point(points[1],{-2,0}));
    points=line.segment_intersections({0,-2},{0,2});
    assert(points.size()==1 && close_point(points[0],{0,0}));
    points=line.segment_intersections({3,1},{3,1});
    assert(points.empty());
    const VoronoiEdge ray{0,1,VoronoiEdgeKind::ray,{0,0},{1,0}};
    assert(!ray.clip({-2,-1},{-1,1}));
    points=ray.segment_intersections({2,0},{-2,0});
    assert(points.size()==2 && close_point(points[0],{2,0}) && close_point(points[1],{0,0}));
    points=ray.segment_intersections({-1,-1},{-1,1});
    assert(points.empty());
    const VoronoiEdge segment{0,1,VoronoiEdgeKind::segment,{-0.5L,0},{1,0}};
    points=segment.segment_intersections({1,0},{0,0});
    assert(points.size()==2 && close_point(points[0],{0.5L,0}) && close_point(points[1],{0,0}));
    points=segment.segment_intersections({0.5L,0},{0.5L,0});
    assert(points.size()==1);
    assert(!segment.clip({2,-1},{3,1}));
    const VoronoiEdge tangent{0,1,VoronoiEdgeKind::ray,{0,0},{1,1}};
    const auto contact=tangent.clip({1,0},{2,1});
    assert(contact && close_point(contact->first,{1,1}) && close_point(contact->second,{1,1}));
}

Real reference_box_radius(const std::vector<Point>& sites,Point lower,Point upper) {
    Real best=-1;
    for(std::size_t owner=0;owner<sites.size();++owner)
        for(Point point:reference_cell(sites,owner,lower,upper)) best=std::max(best,norm2(point-sites[owner]));
    return best;
}

void check_concave_domain(const std::vector<Point>& sites) {
    const std::vector<Point> polygon{{0,0},{8,0},{8,3},{3,3},{3,9},{0,9}};
    const Real expected=std::max(reference_box_radius(sites,{0,0},{8,3}),
                                 reference_box_radius(sites,{0,3},{3,9}));
    const VoronoiDiagram diagram(sites);
    const auto answer=cp_geometry::largest_empty_circle_in_polygon(diagram,polygon);
    assert(answer && close(answer->squared_radius,expected));
    assert(cp_geometry::contains(polygon,answer->center)!=cp_geometry::Location::outside);
    check_nearest(sites,answer->nearest_site,answer->center);
    auto reversed=polygon;
    std::reverse(reversed.begin(),reversed.end());
    const auto other=cp_geometry::largest_empty_circle_in_polygon(diagram,reversed);
    assert(other && close(other->squared_radius,expected));
}

void check_tiny_cells() {
    const Real side=1e-7L;
    const Point lower{0,0},upper{side,side};
    const auto relative_close=[](Real actual,Real expected) {
        return expected>0 && std::abs(actual-expected)<=expected*1e-7L;
    };
    const VoronoiDiagram one({{0,0}});
    const auto cell=one.cell(0,lower,upper);
    assert(cell.size()>=3 && relative_close(polygon_area(cell),side*side));
    const auto answer=cp_geometry::largest_empty_circle_in_box(one,lower,upper);
    assert(answer && relative_close(answer->squared_radius,2*side*side));
    const VoronoiDiagram pair({{0,0},{side,0}});
    const auto cells=pair.cells(lower,upper);
    assert(cells.size()==2);
    for(const auto& piece:cells) assert(piece.size()>=3 && relative_close(polygon_area(piece),side*side/2));
    const auto second_answer=cp_geometry::largest_empty_circle_in_box(pair,lower,upper);
    assert(second_answer && relative_close(second_answer->squared_radius,side*side*1.25L));
    const std::vector<cp_geometry::HalfPlane> constraints{{-1,0,0},{1,0,side},{0,-1,0},{0,1,side}};
    const auto bounded=cp_geometry::half_plane_intersection(constraints,{-1,-1},{1,1});
    assert(bounded.size()>=3 && relative_close(polygon_area(bounded),side*side));
    const auto unbounded_api=cp_geometry::half_plane_intersection(constraints);
    assert(unbounded_api.kind==cp_geometry::RegionKind::POLYGON);
    assert(relative_close(polygon_area(unbounded_api.polygon),side*side));
}

void check_special() {
    check_tiny_cells();
    check_diagram({});
    check_diagram({{1,2}});
    check_diagram({{1,2},{1,2},{1,2}});
    check_diagram({{-1,0},{1,0}});
    const VoronoiDiagram pair({{-1,0},{1,0}});
    assert(pair.vertices().empty() && pair.edges().size()==1);
    assert(pair.edges()[0].kind==VoronoiEdgeKind::line);
    std::vector<Point> line;
    for(int i=-4;i<=4;++i) line.push_back({Real(i),Real(i)*2+1});
    check_diagram(line);
    const VoronoiDiagram collinear(line);
    assert(collinear.vertices().empty() && collinear.edges().size()==line.size()-1);
    for(const auto& edge:collinear.edges()) assert(edge.kind==VoronoiEdgeKind::line);
    const std::vector<Point> square{{-2,-2},{2,-2},{2,2},{-2,2}};
    check_diagram(square);
    const VoronoiDiagram cocircular(square);
    assert(cocircular.vertices().size()==1 && cocircular.edges().size()==4);
    for(const auto& edge:cocircular.edges()) assert(edge.kind==VoronoiEdgeKind::ray);
    check_diagram({{-2,-2},{2,-2},{2,2},{-2,2},{0,0}});
    check_diagram({{-25,-30},{-20,-10},{12,20},{45,-18}});
    std::vector<Point> grid;
    for(int x=-3;x<=3;++x) for(int y=-2;y<=2;++y) grid.push_back({Real(x),Real(y)});
    check_diagram(grid);
    auto scaled=grid;
    for(Point& point:scaled) point=point*0.01L;
    check_diagram(scaled,{-0.08L,-0.06L},{0.09L,0.07L});
    std::vector<Point> translated=grid;
    const Point offset{1000000000.L,-1000000000.L};
    for(Point& point:translated) point=point+offset;
    check_diagram(translated,offset+Point{-8,-6},offset+Point{9,7});
    const VoronoiDiagram triangle({{0,0},{2,0},{0,2}});
    assert(triangle.vertices().size()==1 && close_point(triangle.vertices()[0].point,{1,1}));
    const VoronoiDiagram no_sites({});
    assert(!cp_geometry::largest_empty_circle_in_polygon(no_sites,{{0,0},{1,0},{0,1}}));
    // ICPC WF2018 G statement sample 1, with transmitters at polygon vertices.
    const std::vector<Point> polygon{{0,0},{170,0},{140,30},{60,30},{0,70}};
    const auto answer=cp_geometry::largest_empty_circle_in_polygon(VoronoiDiagram(polygon),polygon);
    assert(answer && close(answer->squared_radius,2500));
    auto reversed=polygon;
    std::reverse(reversed.begin(),reversed.end());
    const auto reverse_answer=cp_geometry::largest_empty_circle_in_polygon(VoronoiDiagram(polygon),reversed);
    assert(reverse_answer && close(reverse_answer->squared_radius,answer->squared_radius));
}
}

int main() {
    review::check_edge_primitives();
    review::check_special();
    std::mt19937 random(721391);
    for(int test=0;test<260;++test) {
        const std::size_t count=1+random()%35;
        std::vector<review::Point> sites;
        for(std::size_t i=0;i<count;++i) sites.push_back({review::Real(int(random()%51)-25),review::Real(int(random()%51)-25)});
        review::check_diagram(sites);
        if(test%20==0) {
            review::check_concave_domain(sites);
            std::reverse(sites.begin(),sites.end());
            review::check_diagram(sites);
        }
    }
}
