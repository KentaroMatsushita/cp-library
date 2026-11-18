#pragma once
#include "detail/macro_push.hpp"
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC push_options
#pragma GCC optimize("no-fast-math")
#endif
#include "base.hpp"
#include <algorithm>
#include <array>
#include <cfenv>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace cp_geometry {
struct DelaunayTriangulation {
    std::vector<Point> sites; // Exact distinct coordinates, lexicographically sorted.
    std::vector<std::size_t> representative; // Smallest original index of each site.
    std::vector<std::size_t> original_to_site;
    std::vector<std::array<std::size_t,3>> triangles; // CCW; smallest vertex first.
    std::vector<std::pair<std::size_t,std::size_t>> edges; // Unique, smaller index first.
    std::vector<std::size_t> hull; // CCW, retaining boundary-collinear sites.
};

namespace delaunay_detail {
// Error-free expansion arithmetic, written from the arithmetic identities.
// Components are stored from least to most significant; zero components vanish.
using Expansion = std::vector<Real>;
inline std::pair<Real,Real> two_sum(Real a,Real b) {
    const Real sum=a+b, recovered_b=sum-a;
    const Real error=(a-(sum-recovered_b))+(b-recovered_b);
    return {sum,error};
}
inline void grow(Expansion& expansion,Real value) {
    Expansion next;
    next.reserve(expansion.size()+1);
    for(Real component:expansion) {
        const auto result=two_sum(value,component);
        if(result.second!=0) next.push_back(result.second);
        value=result.first;
    }
    if(value!=0) next.push_back(value);
    expansion.swap(next);
}
inline Expansion difference(Real a,Real b) {
    const auto result=two_sum(a,-b);
    Expansion value;
    if(result.second!=0) value.push_back(result.second);
    if(result.first!=0) value.push_back(result.first);
    return value;
}
inline Expansion add(Expansion a,const Expansion& b,int multiplier=1) {
    for(Real component:b) grow(a,multiplier*component);
    return a;
}
inline Expansion multiply(const Expansion& a,const Expansion& b) {
    Expansion result;
    for(Real left:a) for(Real right:b) {
        const Real product=left*right;
        const Real error=std::fma(left,right,-product);
        if(error!=0) grow(result,error);
        if(product!=0) grow(result,product);
    }
    return result;
}
inline Expansion determinant(const Expansion& ax,const Expansion& ay,
                             const Expansion& bx,const Expansion& by) {
    return add(multiply(ax,by),multiply(ay,bx),-1);
}
inline int expansion_sign(const Expansion& expansion) {
    return expansion.empty()?0:(expansion.back()>0?1:-1);
}

class Predicates {
    std::vector<Point> points_;
public:
    explicit Predicates(const std::vector<Point>& points) {
        if(std::numeric_limits<Real>::radix!=2 || std::fegetround()!=FE_TONEAREST)
            throw std::invalid_argument("Delaunay predicates require binary round-to-nearest arithmetic");
        Real magnitude=0;
        for(Point point:points) magnitude=(std::max)({magnitude,std::abs(point.x),std::abs(point.y)});
        int exponent=0;
        if(magnitude!=0) std::frexp(magnitude,&exponent);
        // Degree-four predicates need room for all four significands and their
        // residuals. Reject exponent spreads that could underflow exact products.
        const int safe_exponent=(std::numeric_limits<Real>::min_exponent+
                                 8*std::numeric_limits<Real>::digits)/4;
        const Real threshold=std::ldexp(Real(1),safe_exponent);
        points_.reserve(points.size());
        for(Point point:points) {
            Point normalized{std::ldexp(point.x,-exponent),std::ldexp(point.y,-exponent)};
            for(auto coordinate:std::array<std::pair<Real,Real>,2>{{{point.x,normalized.x},{point.y,normalized.y}}})
                if(coordinate.first!=0 && std::abs(coordinate.second)<threshold)
                    throw std::overflow_error("Delaunay coordinate exponent spread is too large");
            points_.push_back(normalized);
        }
    }
    int orient(std::size_t i,std::size_t j,std::size_t k) const {
        const Point a=points_[i],b=points_[j],c=points_[k];
        const Real ax=a.x-c.x,ay=a.y-c.y,bx=b.x-c.x,by=b.y-c.y;
        const Real positive=ax*by,negative=ay*bx,result=positive-negative;
        const Real error=64*std::numeric_limits<Real>::epsilon()*(std::abs(positive)+std::abs(negative));
        if(std::abs(result)>error) return result>0?1:-1;
        return expansion_sign(determinant(difference(a.x,c.x),difference(a.y,c.y),
                                          difference(b.x,c.x),difference(b.y,c.y)));
    }
    // Positive when d is inside the oriented circle of CCW (a,b,c).
    int incircle(std::size_t ia,std::size_t ib,std::size_t ic,std::size_t id) const {
        const Point a=points_[ia],b=points_[ib],c=points_[ic],d=points_[id];
        const Real ax=a.x-d.x,ay=a.y-d.y,bx=b.x-d.x,by=b.y-d.y,cx=c.x-d.x,cy=c.y-d.y;
        const Real aa=ax*ax+ay*ay,bb=bx*bx+by*by,cc=cx*cx+cy*cy;
        const Real bc=bx*cy-by*cx,ca=cx*ay-cy*ax,ab=ax*by-ay*bx;
        const Real result=aa*bc+bb*ca+cc*ab;
        const Real permanent=aa*(std::abs(bx*cy)+std::abs(by*cx))+
                             bb*(std::abs(cx*ay)+std::abs(cy*ax))+
                             cc*(std::abs(ax*by)+std::abs(ay*bx));
        if(std::abs(result)>512*std::numeric_limits<Real>::epsilon()*permanent)
            return result>0?1:-1;
        const auto ex=difference(a.x,d.x),ey=difference(a.y,d.y);
        const auto fx=difference(b.x,d.x),fy=difference(b.y,d.y);
        const auto gx=difference(c.x,d.x),gy=difference(c.y,d.y);
        const auto alift=add(multiply(ex,ex),multiply(ey,ey));
        const auto blift=add(multiply(fx,fx),multiply(fy,fy));
        const auto clift=add(multiply(gx,gx),multiply(gy,gy));
        return expansion_sign(add(add(multiply(alift,determinant(fx,fy,gx,gy)),
                                      multiply(blift,determinant(gx,gy,ex,ey))),
                                  multiply(clift,determinant(ex,ey,fx,fy))));
    }
};

// Each undirected edge owns four rotated handles. Only primal handles (0,2)
// have a site origin. Onext stores a CCW ring around the origin.
class QuadEdges {
    struct Node { std::size_t next=0,origin=0; };
    std::vector<Node> nodes_;
    std::vector<unsigned char> alive_;
    std::vector<std::size_t> free_;
    const Predicates& predicates_;
    static std::size_t rot(std::size_t edge) { return (edge&~std::size_t(3))|((edge+1)&3); }
    static std::size_t invrot(std::size_t edge) { return (edge&~std::size_t(3))|((edge+3)&3); }
public:
    explicit QuadEdges(const Predicates& predicates,std::size_t count):predicates_(predicates) {
        nodes_.reserve(16*count);
        alive_.reserve(4*count);
    }
    static std::size_t reverse(std::size_t edge) { return edge^2; }
    std::size_t origin(std::size_t edge) const { return nodes_[edge].origin; }
    std::size_t destination(std::size_t edge) const { return origin(reverse(edge)); }
    std::size_t next_origin(std::size_t edge) const { return nodes_[edge].next; }
    std::size_t previous_origin(std::size_t edge) const { return rot(next_origin(rot(edge))); }
    std::size_t next_left(std::size_t edge) const { return rot(next_origin(invrot(edge))); }
    std::size_t make_edge(std::size_t a,std::size_t b) {
        std::size_t edge;
        if(free_.empty()) { edge=nodes_.size(); nodes_.resize(edge+4); alive_.push_back(1); }
        else { edge=free_.back(); free_.pop_back(); alive_[edge/4]=1; }
        nodes_[edge]={edge,a}; nodes_[edge+1]={edge+3,0};
        nodes_[edge+2]={edge+2,b}; nodes_[edge+3]={edge+1,0};
        return edge;
    }
    void splice(std::size_t a,std::size_t b) {
        const std::size_t alpha=rot(next_origin(a)),beta=rot(next_origin(b));
        std::swap(nodes_[a].next,nodes_[b].next);
        std::swap(nodes_[alpha].next,nodes_[beta].next);
    }
    std::size_t connect(std::size_t a,std::size_t b) {
        const auto edge=make_edge(destination(a),origin(b));
        splice(edge,next_left(a));
        splice(reverse(edge),b);
        return edge;
    }
    void erase(std::size_t edge) {
        splice(edge,previous_origin(edge));
        splice(reverse(edge),previous_origin(reverse(edge)));
        const std::size_t base=edge&~std::size_t(3);
        alive_[base/4]=0; free_.push_back(base);
    }
    bool left_of(std::size_t vertex,std::size_t edge) const {
        return predicates_.orient(origin(edge),destination(edge),vertex)>0;
    }
    bool right_of(std::size_t vertex,std::size_t edge) const {
        return predicates_.orient(origin(edge),destination(edge),vertex)<0;
    }
    std::pair<std::size_t,std::size_t> build(std::size_t left,std::size_t right) {
        if(right-left==2) {
            const auto edge=make_edge(left,left+1);
            return {edge,reverse(edge)};
        }
        if(right-left==3) {
            const auto a=make_edge(left,left+1),b=make_edge(left+1,left+2);
            splice(reverse(a),b);
            const int orientation=predicates_.orient(left,left+1,left+2);
            if(orientation==0) return {a,reverse(b)};
            const auto c=connect(b,a);
            return orientation>0?std::make_pair(a,reverse(b)):std::make_pair(reverse(c),c);
        }
        const auto middle=left+(right-left)/2;
        auto lhs=build(left,middle),rhs=build(middle,right);
        auto left_outer=lhs.first,left_inner=lhs.second,right_inner=rhs.first,right_outer=rhs.second;
        while(true) {
            if(left_of(origin(right_inner),left_inner)) left_inner=next_left(left_inner);
            else if(right_of(origin(left_inner),right_inner)) right_inner=next_origin(reverse(right_inner));
            else break;
        }
        auto base=connect(reverse(right_inner),left_inner);
        if(origin(left_inner)==origin(left_outer)) left_outer=reverse(base);
        if(origin(right_inner)==origin(right_outer)) right_outer=base;
        while(true) {
            auto lc=next_origin(reverse(base));
            if(right_of(destination(lc),base)) {
                while(predicates_.incircle(destination(base),origin(base),destination(lc),
                                           destination(next_origin(lc)))>0) {
                    const auto successor=next_origin(lc); erase(lc); lc=successor;
                }
            }
            auto rc=previous_origin(base);
            if(right_of(destination(rc),base)) {
                while(predicates_.incircle(destination(base),origin(base),destination(rc),
                                           destination(previous_origin(rc)))>0) {
                    const auto predecessor=previous_origin(rc); erase(rc); rc=predecessor;
                }
            }
            const bool left_valid=right_of(destination(lc),base),right_valid=right_of(destination(rc),base);
            if(!left_valid&&!right_valid) break;
            if(!left_valid || (right_valid && predicates_.incircle(destination(lc),origin(lc),origin(rc),destination(rc))>0))
                base=connect(rc,reverse(base));
            else base=connect(reverse(base),reverse(lc));
        }
        return {left_outer,right_outer};
    }
    void extract(DelaunayTriangulation& result) const {
        std::vector<unsigned char> seen(nodes_.size(),0);
        for(std::size_t group=0;group<alive_.size();++group) if(alive_[group]) {
            const auto edge=group*4;
            auto a=origin(edge),b=destination(edge);
            if(a>b) std::swap(a,b);
            result.edges.push_back({a,b});
            for(auto start:{edge,reverse(edge)}) if(!seen[start]) {
                std::vector<std::size_t> face;
                auto current=start;
                do {
                    if(seen[current]) throw std::runtime_error("invalid Delaunay face topology");
                    seen[current]=1; face.push_back(origin(current));
                    current=next_left(current);
                } while(current!=start);
                if(face.size()==3 && predicates_.orient(face[0],face[1],face[2])>0) {
                    std::rotate(face.begin(),std::min_element(face.begin(),face.end()),face.end());
                    result.triangles.push_back({face[0],face[1],face[2]});
                }
            }
        }
        std::sort(result.edges.begin(),result.edges.end());
        std::sort(result.triangles.begin(),result.triangles.end());
    }
};
} // namespace delaunay_detail

// Ordinary (unweighted) planar Delaunay triangulation. O(N log N) time and O(N)
// memory. Cocircular sites choose one deterministic triangulation, not all edges.
// Requires binary round-to-nearest arithmetic without fast-math; see docs.
inline DelaunayTriangulation delaunay_triangulation(const std::vector<Point>& points) {
#if defined(__FAST_MATH__) && (!defined(__GNUC__) || defined(__clang__))
    throw std::invalid_argument("Delaunay predicates require compilation without fast-math");
#endif
    DelaunayTriangulation result;
    result.original_to_site.resize(points.size());
    std::vector<std::size_t> order(points.size());
    std::iota(order.begin(),order.end(),std::size_t(0));
    for(Point point:points) if(!std::isfinite(point.x)||!std::isfinite(point.y))
        throw std::invalid_argument("Delaunay sites must be finite");
    std::sort(order.begin(),order.end(),[&](std::size_t i,std::size_t j) {
        if(points[i].x!=points[j].x) return points[i].x<points[j].x;
        if(points[i].y!=points[j].y) return points[i].y<points[j].y;
        return i<j;
    });
    for(std::size_t original:order) {
        Point point=points[original];
        if(result.sites.empty() || point.x!=result.sites.back().x || point.y!=result.sites.back().y) {
            result.sites.push_back(point); result.representative.push_back(original);
        }
        result.original_to_site[original]=result.sites.size()-1;
    }
    const auto count=result.sites.size();
    if(count<2) { if(count) result.hull.push_back(0); return result; }
    const delaunay_detail::Predicates predicates(result.sites);
    bool collinear=true;
    for(std::size_t i=1;i+1<count;++i) if(predicates.orient(0,count-1,i)!=0) { collinear=false; break; }
    if(collinear) {
        result.hull.resize(count); std::iota(result.hull.begin(),result.hull.end(),std::size_t(0));
        for(std::size_t i=1;i<count;++i) result.edges.push_back({i-1,i});
        return result;
    }
    // Boundary-collinear sites are kept: each boundary edge has a Voronoi ray.
    std::vector<std::size_t> lower,upper;
    for(std::size_t i=0;i<count;++i) {
        while(lower.size()>1 && predicates.orient(lower[lower.size()-2],lower.back(),i)<0) lower.pop_back();
        lower.push_back(i);
    }
    for(std::size_t i=count;i-->0;) {
        while(upper.size()>1 && predicates.orient(upper[upper.size()-2],upper.back(),i)<0) upper.pop_back();
        upper.push_back(i);
    }
    lower.pop_back(); upper.pop_back();
    result.hull=std::move(lower); result.hull.insert(result.hull.end(),upper.begin(),upper.end());
    delaunay_detail::QuadEdges graph(predicates,count);
    graph.build(0,count);
    graph.extract(result);
    return result;
}
} // namespace cp_geometry
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC pop_options
#endif
#include "detail/macro_pop.hpp"
