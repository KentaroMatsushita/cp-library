template <class T, class F> struct UnionFind {
    vector<T> data;
    vector<int> sz;
    int num;
    const F f;
    T unit;

    UnionFind(int n, const F &f = plus<T>(), const T &unit = T()) : num(n), f(f), unit(unit), sz(n, -1), data(n, unit) {}

    bool unite(int x, int y) {
        x = find(x), y = find(y);
        if(x == y) return false;
        num--;
        if(-sz[x] < -sz[y]) swap(x, y);
        f(data[x], data[y]);
        sz[x] += sz[y];
        sz[y] = x;
        return true;
    }

    bool same(int x, int y) { return find(x) == find(y); }

    int find(int x) {
        if(sz[x] < 0) return x;
        return (sz[x] = find(sz[x]));
    }
    int size(int x) { return -sz[find(x)]; }
    T &get(int x) { return data[find(x)]; }
    const int operator[](const int k) { return find(k); }
    vector<vector<int>> belong() {
        vector<vector<int>> res(data.size());
        for(int i = 0; i < data.size(); i++) res[find(i)].emplace_back(i);
        return res;
    }
};
