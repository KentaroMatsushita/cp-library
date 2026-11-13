#pragma once
#include <algorithm>
#include <cmath>
#include <deque>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

// Floating point geometry. All polygons returned here are counterclockwise.
namespace cp_geometry {
using Real = long double;
inline constexpr Real eps = 1e-12L;
inline int sign(Real x) { return (x > eps) - (x < -eps); }
struct Point {
    Real x = 0, y = 0;
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
    return sign(cross(b-a,p-a))==0 && dot(p-a,p-b)<=eps;
}
struct Line {
    Point a,b;
    Line(Point a, Point b):a(a),b(b) {
        if(near(a,b)) throw std::invalid_argument("line needs two distinct points");
    }
    Point direction() const { return b-a; }
};
// Requires nonparallel lines.
inline Point line_intersection(const Line& l, const Line& m) {
    Real d=cross(l.direction(),m.direction());
    if(sign(d)==0) throw std::invalid_argument("parallel lines");
    return l.a+l.direction()*(cross(m.a-l.a,m.direction())/d);
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
    Real s=0; for(size_t i=0;i<p.size();++i) s+=cross(p[i],p[(i+1)%p.size()]);
    return s;
}
// Exact lexicographic ordering is intentional: epsilon ordering is not transitive.
inline std::vector<Point> convex_hull(std::vector<Point> p) {
    std::sort(p.begin(),p.end(),[](Point a,Point b){return a.x!=b.x?a.x<b.x:a.y<b.y;});
    p.erase(std::unique(p.begin(),p.end(),[](Point a,Point b){return a.x==b.x&&a.y==b.y;}),p.end());
    if(p.size()<3) return p;
    std::vector<Point> h;
    for(Point q:p) {
        while(h.size()>1 && sign(cross(h.back()-h[h.size()-2],q-h.back()))<=0) h.pop_back();
        h.push_back(q);
    }
    size_t lower=h.size();
    for(int i=int(p.size())-2;i>=0;--i) {
        while(h.size()>lower && sign(cross(h.back()-h[h.size()-2],p[i]-h.back()))<=0) h.pop_back();
        h.push_back(p[i]);
    }
    h.pop_back(); return h;
}
enum Location { OUT=0, ON=1, IN=2 };
inline Location contains(const std::vector<Point>& poly,Point p) {
    bool in=false;
    for(size_t i=0;i<poly.size();++i) {
        Point a=poly[i],b=poly[(i+1)%poly.size()];
        if(on_segment(a,b,p)) return ON;
        if(a.y>b.y) std::swap(a,b);
        if(a.y<=p.y && p.y<b.y && cross(b-a,p-a)>0) in=!in;
    }
    return in?IN:OUT;
}
inline int angle_half(Point p) { return p.y<0 || (p.y==0&&p.x<0); }
inline bool angle_less(Point a,Point b) {
    if(angle_half(a)!=angle_half(b)) return angle_half(a)<angle_half(b);
    return cross(a,b)>0;
}

class ConvexPolygon {
    std::vector<Point> p_;
    std::vector<std::pair<Point,int>> edges_;
public:
    // Accepts an arbitrary point cloud; builds a strict hull in O(N log N).
    explicit ConvexPolygon(std::vector<Point> points):p_(convex_hull(std::move(points))) {
        if(p_.size()<3) return;
        for(int i=0;i<int(p_.size());++i) edges_.push_back({p_[(i+1)%p_.size()]-p_[i],i});
        std::sort(edges_.begin(),edges_.end(),[](const auto& a,const auto& b){return angle_less(a.first,b.first);});
    }
    const std::vector<Point>& vertices() const { return p_; }
    Location contains(Point q) const {
        int n=int(p_.size());
        if(n==0) return OUT;
        if(n==1) return near(p_[0],q)?ON:OUT;
        if(n==2) return on_segment(p_[0],p_[1],q)?ON:OUT;
        int a=sign(cross(p_[1]-p_[0],q-p_[0]));
        int b=sign(cross(p_.back()-p_[0],q-p_[0]));
        if(a<0||b>0) return OUT;
        if(a==0) return on_segment(p_[0],p_[1],q)?ON:OUT;
        if(b==0) return on_segment(p_[0],p_.back(),q)?ON:OUT;
        int l=1,r=n-1;
        while(r-l>1) {int m=(l+r)/2; if(cross(p_[m]-p_[0],q-p_[0])>=0) l=m; else r=m;}
        int s=sign(cross(p_[r]-p_[l],q-p_[l]));
        return s<0?OUT:s==0?ON:IN;
    }
    // Index maximizing dot(direction,vertex). O(log N), ties arbitrary.
    int support(Point direction) const {
        if(p_.empty()) throw std::invalid_argument("empty polygon");
        if(p_.size()<3 || norm2(direction)==0) {
            int k=0; for(int i=1;i<int(p_.size());++i) if(dot(direction,p_[i])>dot(direction,p_[k])) k=i;
            return k;
        }
        Point target=perpendicular(direction);
        auto it=std::lower_bound(edges_.begin(),edges_.end(),target,
            [](const auto& e,Point t){return angle_less(e.first,t);});
        return (it==edges_.end()?edges_.front():*it).second;
    }
    // Boundary intersection: zero, one, or two points, ordered along the line.
    // A coinciding polygon edge returns its two endpoints. O(log N).
    std::vector<Point> line_intersections(const Line& line) const {
        std::vector<Point> out;
        int n=int(p_.size()); if(n==0) return out;
        Point v=line.direction();
        auto f=[&](int i){return cross(v,p_[i]-line.a);};
        auto add=[&](Point q){for(Point r:out) if(near(q,r)) return; out.push_back(q);};
        if(n<3) {
            for(Point q:p_) if(sign(cross(v,q-line.a))==0) add(q);
            if(n==2 && sign(f(0))*sign(f(1))<0) add(p_[0]+(p_[1]-p_[0])*(f(0)/(f(0)-f(1))));
        } else {
            int hi=support(perpendicular(v)),lo=support(perpendicular(v)*-1);
            if(sign(f(lo))>0||sign(f(hi))<0) return out;
            for(int step:{-1,1}) {
                auto index=[&](int k){return (lo+step*k+n)%n;};
                int len=step==1?(hi-lo+n)%n:(lo-hi+n)%n;
                int l=0,r=len;
                while(r-l>1) {int m=(l+r)/2; if(f(index(m))<0) l=m; else r=m;}
                int a=index(l),b=index(r);
                if(sign(f(a))==0) add(p_[a]);
                if(sign(f(b))==0) add(p_[b]);
                if(sign(f(a))*sign(f(b))<0) add(p_[a]+(p_[b]-p_[a])*(f(a)/(f(a)-f(b))));
            }
            // Support plateaus have at most two vertices in a strict hull.
            for(int k:{lo,hi}) if(sign(f(k))==0) {
                add(p_[k]);
                for(int j:{(k+1)%n,(k+n-1)%n}) if(sign(f(j))==0) add(p_[j]);
            }
        }
        std::sort(out.begin(),out.end(),[&](Point a,Point b){return dot(v,a)<dot(v,b);});
        return out;
    }
};

// Rotating calipers on a strict CCW convex hull. O(N).
inline Real convex_diameter2(const std::vector<Point>& p) {
    int n=int(p.size()); if(n<2) return 0;
    if(n==2) return norm2(p[0]-p[1]);
    Real ans=0; int j=1;
    for(int i=0;i<n;++i) {
        Point e=p[(i+1)%n]-p[i];
        while(cross(e,p[(j+1)%n]-p[i])>cross(e,p[j]-p[i])+eps) j=(j+1)%n;
        ans=std::max({ans,norm2(p[i]-p[j]),norm2(p[(i+1)%n]-p[j])});
    }
    return ans;
}
// Clip a convex polygon by the closed left half-plane of a directed line. O(N).
inline std::vector<Point> convex_cut(const std::vector<Point>& p,const Line& l) {
    std::vector<Point> q;
    for(size_t i=0;i<p.size();++i) {
        Point a=p[i],b=p[(i+1)%p.size()];
        Real x=cross(l.direction(),a-l.a),y=cross(l.direction(),b-l.a);
        if(sign(x)>=0) q.push_back(a);
        if(sign(x)*sign(y)<0) q.push_back(a+(b-a)*(x/(x-y)));
    }
    return q;
}
// Intersect left half-planes inside an EXPLICIT finite bounding rectangle.
// Returns only positive-area intersections; empty/point/segment => {}.
// O(N log N), including sorting. The rectangle is part of the constraints.
inline std::vector<Point> half_plane_intersection(std::vector<Line> lines,
                                                Point lower,Point upper) {
    if(!(lower.x<upper.x && lower.y<upper.y)) throw std::invalid_argument("invalid bounding box");
    Point a=lower,b{upper.x,lower.y},c=upper,d{lower.x,upper.y};
    lines.insert(lines.end(),{Line(a,b),Line(b,c),Line(c,d),Line(d,a)});
    std::sort(lines.begin(),lines.end(),[](const Line& l,const Line& r){return angle_less(l.direction(),r.direction());});
    std::vector<Line> unique;
    for(const auto& l:lines) {
        if(!unique.empty() && sign(cross(unique.back().direction(),l.direction()))==0 && dot(unique.back().direction(),l.direction())>0) {
            if(cross(unique.back().direction(),l.a-unique.back().a)>0) unique.back()=l;
        } else unique.push_back(l);
    }
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
    for(size_t i=0;i<q.size();++i) result.push_back(meet(q[i],q[(i+1)%q.size()]));
    if(sign(area2(result))<=0) return {};
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
    for(auto& l:lines) l.b=l.a+l.direction()/std::sqrt(norm2(l.direction()));
    std::mt19937 rng(712367);
    std::shuffle(lines.begin(),lines.end(),rng);
    Point feasible{};
    for(size_t i=0;i<lines.size();++i) {
        const auto& l=lines[i];
        if(sign(cross(l.direction(),feasible-l.a))>=0) continue;
        Real left=-std::numeric_limits<Real>::infinity(),right=-left;
        Point v=l.direction();
        for(size_t j=0;j<i;++j) {
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
    for(size_t i=0;i<lines.size();++i) {
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
    for(size_t i=0;i<q.size();++i) result.push_back(meet(q[i],q[(i+1)%q.size()]));
    if(sign(area2(result))<=0) return {RegionKind::DEGENERATE,feasible,{}};
    return {RegionKind::POLYGON,feasible,std::move(result)};
}
} // namespace cp_geometry
