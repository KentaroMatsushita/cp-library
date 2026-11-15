#pragma once
#include "detail/macro_push.hpp"
#include "base.hpp"
#include <algorithm>
#include <utility>
#include <vector>

namespace cp_geometry {
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
    std::size_t lower=h.size();
    for(int i=int(p.size())-2;i>=0;--i) {
        while(h.size()>lower && sign(cross(h.back()-h[h.size()-2],p[i]-h.back()))<=0) h.pop_back();
        h.push_back(p[i]);
    }
    h.pop_back(); return h;
}
inline int angle_half(Point p) { return p.y<0 || (p.y==0&&p.x<0); }
inline bool angle_less(Point a,Point b) {
    if(angle_half(a)!=angle_half(b)) return angle_half(a)<angle_half(b);
    return cross(a,b)>0;
}

class ConvexPolygon {
    std::vector<Point> p_;
    std::vector<std::pair<Point,int>> edges_;
    void index_edges() {
        edges_.clear();
        if(p_.size()<3) return;
        edges_.reserve(p_.size());
        for(int i=0;i<int(p_.size());++i) edges_.push_back({p_[(i+1)%p_.size()]-p_[i],i});
        auto first=std::min_element(edges_.begin(),edges_.end(),
            [](const auto& a,const auto& b){return angle_less(a.first,b.first);});
        std::rotate(edges_.begin(),first,edges_.end());
    }
public:
    // Arbitrary point cloud: strict CCW hull and index, O(N log N).
    explicit ConvexPolygon(std::vector<Point> points):p_(convex_hull(std::move(points))) { index_edges(); }
    // O(N) construction from an already strict convex hull, CW or CCW.
    // No repeated closing vertex or intermediate collinear vertices are allowed.
    // The caller must supply a simple convex polygon, not a general point cloud.
    static ConvexPolygon from_convex_hull(std::vector<Point> hull) {
        ConvexPolygon result({});
        if(hull.size()>2) {
            if(area2(hull)<0) std::reverse(hull.begin(),hull.end());
            for(std::size_t i=0;i<hull.size();++i)
                if(sign(cross(hull[(i+1)%hull.size()]-hull[i],
                              hull[(i+2)%hull.size()]-hull[(i+1)%hull.size()]))<=0)
                    throw std::invalid_argument("strict convex hull required");
        }
        result.p_=std::move(hull);
        result.index_edges();
        return result;
    }
    const std::vector<Point>& vertices() const { return p_; }
    Location contains(Point q) const {
        int n=int(p_.size());
        if(n==0) return Location::outside;
        if(n==1) return near(p_[0],q)?Location::boundary:Location::outside;
        if(n==2) return on_segment(p_[0],p_[1],q)?Location::boundary:Location::outside;
        int a=sign(cross(p_[1]-p_[0],q-p_[0]));
        int b=sign(cross(p_.back()-p_[0],q-p_[0]));
        if(a<0||b>0) return Location::outside;
        if(a==0) return on_segment(p_[0],p_[1],q)?Location::boundary:Location::outside;
        if(b==0) return on_segment(p_[0],p_.back(),q)?Location::boundary:Location::outside;
        int l=1,r=n-1;
        while(r-l>1) {int m=(l+r)/2; if(cross(p_[m]-p_[0],q-p_[0])>=0) l=m; else r=m;}
        int s=sign(cross(p_[r]-p_[l],q-p_[l]));
        return s<0?Location::outside:s==0?Location::boundary:Location::inside;
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
    for(std::size_t i=0;i<p.size();++i) {
        Point a=p[i],b=p[(i+1)%p.size()];
        Real x=cross(l.direction(),a-l.a),y=cross(l.direction(),b-l.a);
        if(sign(x)>=0) q.push_back(a);
        if(sign(x)*sign(y)<0) q.push_back(a+(b-a)*(x/(x-y)));
    }
    return q;
}

} // namespace cp_geometry
#include "detail/macro_pop.hpp"
