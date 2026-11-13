template <typename T = ll> struct Point {
    T x, y;
    constexpr Point operator+(const Point rhs) const noexcept { return Point(*this) += rhs; }
    constexpr Point operator-(const Point rhs) const noexcept { return Point(*this) -= rhs; }
    constexpr Point &operator+=(const Point rhs) noexcept {
        x += rhs.x, y += rhs.y;
        return *this;
    }
    constexpr Point &operator-=(const Point rhs) noexcept {
        x -= rhs.x, y -= rhs.y;
        return *this;
    }
    template <typename S> constexpr Point &operator*=(const S a) noexcept {
        x *= a, y *= a;
        return *this;
    }
    template <typename S> constexpr Point operator*(const S a) const noexcept { return Point(*this) *= a; }
    template <typename S> constexpr Point &operator/=(const S a) noexcept {
        x /= a, y /= a;
        return *this;
    }
    template <typename S> constexpr Point operator/(const S a) const noexcept { return Point(*this) /= a; }
    constexpr bool operator<(const Point rhs) const noexcept {
        if(!x and !y) return rhs.x != 0 or rhs.y != 0;
        if(!rhs.x and !rhs.y) return false;
        bool u = (y > 0 or (x < 0 and y == 0)), v = (rhs.y > 0 or (rhs.x < 0 and rhs.y == 0));
        if(u != v) return u;
        return x * rhs.y > y * rhs.x;
    }
    Point(T x, T y) : x(x), y(y) {}
    Point() = default;
    void print() { cout << x << " " << y << endl; }
};
template <typename T = ll> using Points = vector<Point<T>>;

// ソートは x >= 0 y > 0 から反時計回り
