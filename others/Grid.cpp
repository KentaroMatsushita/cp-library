struct Grid {
    int n, m;
    Grid(int n, int m) : n(n), m(m) {}
    struct GridSquare {
        const int &n, &m;
        int x, y;
        GridSquare(const int &n, const int &m, int x, int y) : x(x), y(y), n(n), m(m) {}
        struct itr {
            const int &n, &m, &x, &y;
            int i;
            itr(const int &n, const int &m, const int &x, const int &y) : n(n), m(m), x(x), y(y), i(0) {
                while(i < 4 and (!inc(dx4[i].fi + x, 0, n) or !inc(dx4[i].se + y, 0, m))) { i++; }
            }
            bool operator!=(const itr &) const { return i != 4; }
            void operator++() {
                while(++i < 4) {
                    if(inc(dx4[i].fi + x, 0, n) and inc(dx4[i].se + y, 0, m)) break;
                }
            }
            pair<int, int> operator*() const { return pair(x + dx4[i].fi, y + dx4[i].se); }
        };
        itr begin() const { return itr(n, m, x, y); }
        itr end() const { return itr(n, m, x, y); }
    };
    GridSquare get(int x, int y) { return GridSquare(n, m, x, y); }
    GridSquare operator()(int x, int y) { return get(x, y); }
    int id(int x, int y) { return x * m + y; }
};
