namespace Geometry {
using T = ll;
inline constexpr int type(T x, T y) {
    if(x == 0 and y == 0) return 0;
    if(y < 0 or (y == 0 and x > 0)) return -1;
    return 1;
}
struct Point {
    T x, y;
    constexpr Point(T _x = 0, T _y = 0) : x(_x), y(_y) {}
    constexpr Point operator+() const noexcept { return *this; }
    constexpr Point operator-() const noexcept { return Point(-x, -y); }
    constexpr Point operator+(const Point &p) const { return Point(x + p.x, y + p.y); }
    constexpr Point operator-(const Point &p) const { return Point(x - p.x, y - p.y); }
    constexpr Point &operator+=(const Point &p) { return x += p.x, y += p.y, *this; }
    constexpr Point &operator-=(const Point &p) { return x -= p.x, y -= p.y, *this; }
    constexpr T operator*(const Point &p) const { return x * p.x + y * p.y; }
    constexpr Point &operator*=(const T &k) { return x *= k, y *= k, *this; }
    constexpr Point operator*(const T &k) const { return Point(x * k, y * k); }
    constexpr bool operator==(const Point &r) const noexcept { return r.x == x and r.y == y; }

    constexpr bool operator<(const Point &r) const {
        int L = type(x, y), R = type(r.x, r.y);
        if(L != R) return L < R;

        T X = x * r.y, Y = r.x * y;
        if(X != Y) return X > Y;
        return x < r.x;
    }
};
ostream &operator<<(ostream &os, const Point &p) { return os << p.x << " " << p.y; }
istream &operator>>(istream &is, Point &p) {
    is >> p.x >> p.y;
    return is;
}

struct Line {
    Point a, b;
    Line() = default;
    Line(Point a, Point b) : a(a), b(b) {}
    // ax + by = c
    Line(T A, T B, T C) {
        if(A == 0) {
            a = Point(0, C / B), b = Point(1, C / B);
        } else if(B == 0) {
            a = Point(C / A, 0), b = Point(C / A, 1);
        } else {
            a = Point(0, C / B), b = Point(C / A, 0);
        }
    }
};
ostream &operator<<(ostream &os, Line &p) { return os << p.a << " to " << p.b; }
istream &operator>>(istream &is, Line &p) { return is >> p.a >> p.b; }

struct Segment : Line {
    Segment() = default;
    Segment(Point a, Point b) : Line(a, b) {}
};

ostream &operator<<(ostream &os, Segment &p) { return os << p.a << " to " << p.b; }
istream &operator>>(istream &is, Segment &p) {
    is >> p.a >> p.b;
    return is;
}

struct Circle {
    Point p;
    T r;
    Circle() = default;
    Circle(Point p, T r) : p(p), r(r) {}
};

using pt = Point;
using Points = vector<pt>;
T cross(const pt &x, const pt &y) { return x.x * y.y - x.y * y.x; }
T dot(const pt &x, const pt &y) { return x.x * y.x + x.y * y.y; }

T abs2(const pt &x) { return dot(x, x); }
// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_1_C
// 点の回転方向
int ccw(const Point &a, Point b, Point c) {
    b = b - a, c = c - a;
    if(cross(b, c) > 0) return +1;   // "COUNTER_CLOCKWISE"
    if(cross(b, c) < 0) return -1;   // "CLOCKWISE"
    if(dot(b, c) < 0) return +2;     // "ONLINE_BACK"
    if(abs2(b) < abs2(c)) return -2; // "ONLINE_FRONT"
    return 0;                        // "ON_SEGMENT"
}

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_2_A
// 平行判定
bool parallel(const Line &a, const Line &b) { return (cross(a.b - a.a, b.b - b.a) == 0); }

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_2_A
// 垂直判定
bool orthogonal(const Line &a, const Line &b) { return (dot(a.a - a.b, b.a - b.b) == 0); }

bool intersect(const Line &l, const Point &p) { return abs(ccw(l.a, l.b, p)) != 1; }

bool intersect(const Line &l, const Line &m) { return !parallel(l, m); }

bool intersect(const Segment &s, const Point &p) { return ccw(s.a, s.b, p) == 0; }

bool intersect(const Line &l, const Segment &s) { return cross(l.b - l.a, s.a - l.a) * cross(l.b - l.a, s.b - l.a) <= 0; }

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_2_B
bool intersect(const Segment &s, const Segment &t) { return ccw(s.a, s.b, t.a) * ccw(s.a, s.b, t.b) <= 0 && ccw(t.a, t.b, s.a) * ccw(t.a, t.b, s.b) <= 0; }

T norm(const Point &p) { return p.x * p.x + p.y * p.y; }
Point projection(const Segment &l, const Point &p) {
    T t = dot(p - l.a, l.a - l.b) / norm(l.a - l.b);
    return l.a + (l.a - l.b) * t;
}

Point crosspoint(const Line &l, const Line &m) {
    T A = cross(l.b - l.a, m.b - m.a);
    T B = cross(l.b - l.a, l.b - m.a);
    if(A == 0 and B == 0) return m.a;
    return m.a + (m.b - m.a) * (B / A);
}

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_2_C
Point crosspoint(const Segment &l, const Segment &m) { return crosspoint(Line(l), Line(m)); }

// http://judge.u-aizu.ac.jp/onlinejudge/description.jsp?id=CGL_3_B
// 凸性判定
bool is_convex(const Points &p) {
    int n = (int)p.size();
    for(int i = 0; i < n; i++) {
        if(ccw(p[(i + n - 1) % n], p[i], p[(i + 1) % n]) == -1) return false;
    }
    return true;
}

Points convex_hull(Points &p) {
    int n = p.size(), k = 0;
    if(n <= 2) return p;
    sort(begin(p), end(p), [](pt x, pt y) { return (x.x != y.x ? x.x < y.x : x.y < y.y); });
    Points ch(2 * n);
    for(int i = 0; i < n; ch[k++] = p[i++]) {
        while(k >= 2 && cross(ch[k - 1] - ch[k - 2], p[i] - ch[k - 1]) <= 0) --k;
    }
    for(int i = n - 2, t = k + 1; i >= 0; ch[k++] = p[i--]) {
        while(k >= t && cross(ch[k - 1] - ch[k - 2], p[i] - ch[k - 1]) <= 0) --k;
    }
    ch.resize(k - 1);
    return ch;
}
T area2(const Points &p) {
    T res = 0;
    rep(i, si(p)) { res += cross(p[i], p[i == si(p) - 1 ? 0 : i + 1]); }
    return res;
}

enum { OUTTT, ONNN, INNN };

int contains(const Points &Q, const Point &p) {
    bool in = false;
    for(int i = 0; i < Q.size(); i++) {
        Point a = Q[i] - p, b = Q[(i + 1) % Q.size()] - p;
        if(a.y > b.y) swap(a, b);
        if(a.y <= 0 && 0 < b.y && cross(a, b) < 0) in = !in;
        if(cross(a, b) == 0 && dot(a, b) <= 0) return ONNN;
    }
    return in ? INNN : OUTTT;
}
} // namespace Geometry

using namespace Geometry;
