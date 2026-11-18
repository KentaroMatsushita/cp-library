#include "geometry/voronoi.hpp"
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>
using namespace cp_geometry;

bool close(Real a,Real b) {return std::abs(a-b)<1e-7L*(1+std::abs(b));}

// Slow all-pairs clipping used only as an independent bounded-domain oracle.
Real reference_maximum(const std::vector<Point>& sites,Point lower,Point upper) {
    Real result=0;
    for(std::size_t i=0;i<sites.size();++i) {
        std::vector<Point> polygon{lower,{upper.x,lower.y},upper,{lower.x,upper.y}};
        for(std::size_t j=0;j<sites.size();++j) if(i!=j) {
            Point normal=sites[j]-sites[i];
            if(norm2(normal)==0) continue;
            const Point midpoint=sites[i]+normal/2;
            polygon=convex_cut(polygon,Line(midpoint,midpoint+perpendicular(normal)));
        }
        for(auto p:polygon) result=std::max(result,norm2(p-sites[i]));
    }
    return result;
}

int main() {
    {
        VoronoiDiagram empty({});
        assert(empty.edges().empty()&&empty.vertices().empty()&&empty.cells({0,0},{1,1}).empty());
        assert(!largest_empty_circle_in_box(empty,{0,0},{1,1}));
        VoronoiDiagram single({Point(1,1)});
        assert(close(polygon_area(single.cell(0,{0,0},{2,2})),4));
        auto circle=largest_empty_circle_in_box(single,{0,0},{2,2});
        assert(circle&&close(circle->squared_radius,2));
    }
    {
        VoronoiDiagram line({Point(-1,0),Point(1,0)});
        assert(line.edges().size()==1&&line.edges()[0].kind==VoronoiEdgeKind::line);
        auto clipped=line.edges()[0].clip({-2,-3},{2,3});
        assert(clipped&&close(clipped->first.x,0)&&close(clipped->second.x,0));
        assert(close(std::abs(clipped->first.y-clipped->second.y),6));
        auto crossings=line.edges()[0].segment_intersections({-2,1},{2,1});
        assert(crossings.size()==1&&near(crossings[0],{0,1}));
        crossings=line.edges()[0].segment_intersections({0,-2},{0,2});
        assert(crossings.size()==2);
        assert(near(crossings[0],{0,-2})&&near(crossings[1],{0,2}));
    }
    {
        std::vector<Point> square{{-1,-1},{1,-1},{1,1},{-1,1}};
        VoronoiDiagram diagram(square);
        assert(diagram.vertices().size()==1&&near(diagram.vertices()[0].point,{0,0}));
        assert(diagram.edges().size()==4);
        for(auto e:diagram.edges()) assert(e.kind==VoronoiEdgeKind::ray);
        auto best=largest_empty_circle_in_polygon(diagram,square);
        assert(best&&near(best->center,{0,0})&&close(best->squared_radius,2));
        std::reverse(square.begin(),square.end());
        best=largest_empty_circle_in_polygon(diagram,square);
        assert(best&&close(best->squared_radius,2));
    }
    {
        // Determinant exactly 1, although its two 1e31 products round equal.
        std::vector<Point> sites{{0,0},{2111485077978050.L,3416454622906707.L},
                                {3416454622906707.L,5527939700884757.L}};
        VoronoiDiagram diagram(sites);
        assert(diagram.vertices().size()==1&&diagram.edges().size()==3);
        const auto center=diagram.vertices()[0].point;
        assert(std::isfinite(center.x)&&std::isfinite(center.y)&&center.x<0&&center.y>0);
        for(const auto& edge:diagram.edges()) {
            assert(edge.kind==VoronoiEdgeKind::ray);
            const auto third=3-edge.first_site-edge.second_site;
            const Point a=sites[edge.first_site],b=sites[edge.second_site],c=sites[third];
            const __int128 ax=static_cast<long long>(b.x)-static_cast<long long>(a.x);
            const __int128 ay=static_cast<long long>(b.y)-static_cast<long long>(a.y);
            const __int128 bx=static_cast<long long>(c.x)-static_cast<long long>(a.x);
            const __int128 by=static_cast<long long>(c.y)-static_cast<long long>(a.y);
            const __int128 orientation=ax*by-ay*bx;
            const Real side=dot(edge.direction,perpendicular(b-a));
            assert((orientation>0&&side<0)||(orientation<0&&side>0));
        }
    }
    std::mt19937 engine(523782);
    const std::vector<Point> concave{{0,0},{6,0},{6,2},{2,2},{2,6},{0,6}};
    for(int trial=0;trial<600;++trial) {
        std::vector<Point> sites;
        const int n=1+engine()%25;
        for(int i=0;i<n;++i) sites.emplace_back(int(engine()%23)-8,int(engine()%23)-8);
        VoronoiDiagram diagram(sites);
        auto best=largest_empty_circle_in_box(diagram,{-3,-4},{8,9});
        assert(best&&close(best->squared_radius,reference_maximum(sites,{-3,-4},{8,9})));
        auto polygon_best=largest_empty_circle_in_polygon(diagram,concave);
        const Real expected=std::max(reference_maximum(sites,{0,0},{6,2}),
                                     reference_maximum(sites,{0,0},{2,6}));
        assert(polygon_best&&close(polygon_best->squared_radius,expected));
        assert(contains(concave,polygon_best->center)!=Location::outside);
        for(auto site:sites) assert(norm2(site-polygon_best->center)+1e-7L>=polygon_best->squared_radius);
    }
    std::cout<<"voronoi: dual segments/rays/lines, cocircular cells and 600 empty-circle comparisons passed\n";
}
