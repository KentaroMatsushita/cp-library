/**
 * Author: Simon Lindholm
 * Date: 2017-04-20
 * License: CC0
 * Source: own work
 * Description: Container where you can add lines of the form kx+m, and query maximum values at points x.
 *  Useful for dynamic programming (``convex hull trick'').
 * Time: O(\log N)
 * Status: stress-tested
 */

#define double long double
struct Line {
    mutable double k, m, p;
    bool operator<(const Line &o) const { return k < o.k; }
    bool operator<(double x) const { return p < x; }
};

template <bool ismin = true> struct CHT : multiset<Line, less<>> {
    // (for doubles, use inf = 1/.0, div(a,b) = a/b)
    const double inf = 1 / .0;
    bool isect(iterator x, iterator y) {
        if(y == end()) {
            x->p = inf;
            return false;
        }
        if(x->k == y->k)
            x->p = x->m > y->m ? inf : -inf;
        else
            x->p = (y->m - x->m) / (x->k - y->k);
        return x->p >= y->p;
    }
    void add(double k, double m) {
        if(ismin) k = -k, m = -m;
        auto z = insert({k, m, 0}), y = z++, x = y;
        while(isect(y, z)) z = erase(z);
        if(x != begin() && isect(--x, y)) isect(x, y = erase(y));
        while((y = x) != begin() && (--x)->p >= y->p) isect(x, erase(y));
    }
    double query(double x) {
        auto l = *lower_bound(x);
        ll s = 1;
        if(ismin) s = -1;
        return s * (l.k * x + l.m);
    }
};