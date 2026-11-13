namespace egcd_impl {
template <typename T, typename S> T ceil(T x, S y) {
    assert(y);
    return (y < 0 ? ceil(-x, -y) : (x > 0 ? (x + y - 1) / y : x / y));
}

template <typename T, typename S> T floor(T x, S y) {
    assert(y);
    return (y < 0 ? floor(-x, -y) : (x > 0 ? x / y : x / y - (x % y == 0 ? 0 : 1)));
}

// (x,y) は "(b / gcd(a,b), a / gcd(a,b))" で潰した時の代表元！！！！
template <typename T = long long> struct EXTgcd {
    static T INF;
    EXTgcd(T _a, T _b, T _INF = -1) : a(_a), b(_b) {
        if(_INF != -1) INF = _INF;
        if(INF == 0) INF = numeric_limits<T>::max() / 4;
        if(a == 0 && b == 0) throw std::invalid_argument("both coefficients are zero");
        g = extgcd(a, b, x, y);
        if(g < 0) g = -g, x = -x, y = -y;
        dx = b / g, dy = -a / g;
    }
    // 解は条件を満たす中で x が最小のものを取る
    pair<bool, pair<T, T>> get(const T &v, const T &min_x, T min_y = -INF, T max_x = INF, T max_y = INF) {
        if(v % g) return {false, {}};
        T k = v / g;
        T tx = x * k, ty = y * k;
        auto p = get_seg(tx, min_x, max_x, dx);
        auto q = get_seg(ty, min_y, max_y, dy);
        T L = max(p.first, q.first), R = min(p.second, q.second);
        if(L <= R) { T t = dx >= 0 ? L : R; return {true, {tx + dx * t, ty + dy * t}}; }
        return {false, {}};
    }
    T count(const T &v, const T &min_x, const T &min_y = -INF, const T &max_x = INF, const T &max_y = INF) {
        if(v % g) return 0;
        T k = v / g;
        T tx = x * k, ty = y * k;
        auto p = get_seg(tx, min_x, max_x, dx);
        auto q = get_seg(ty, min_y, max_y, dy);
        T L = max(p.first, q.first), R = min(p.second, q.second);
        if(L <= R) { return R - L + 1; }
        return 0;
    }

  private:
    T a, b, g, x, y, dx, dy;
    T extgcd(T a, T b, T &x, T &y) {
        if(!b) {
            x = 1, y = 0;
            return a;
        }
        T res = extgcd(b, a % b, y, x);
        y -= a / b * x;
        return res;
    }
    pair<T, T> get_seg(const T &now, const T &lx, const T &rx, const T &dx) {
        if(dx == 0) return lx <= now && now <= rx ? pair<T,T>{-INF,INF} : pair<T,T>{1,0};
        if(dx > 0) { return _get_seg(lx - now, rx - now, dx); }
        return _get_seg(now - rx, now - lx, -dx);
    }
    pair<T, T> _get_seg(const T &lx, const T &rx, const T &dx) { return pair<T, T>(ceil(lx, dx), floor(rx, dx)); }
};
template <typename T> T EXTgcd<T>::INF = 0;
} // namespace egcd_impl
using egcd_impl::EXTgcd;
