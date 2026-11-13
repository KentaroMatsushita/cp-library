struct UnionFind {
    vector<int> data;
    int num;

    UnionFind(int n) : num(n) { data.assign(n, -1); }

    bool unite(int x, int y) {
        x = find(x), y = find(y);
        if(x == y) return false;
        num--;
        if(data[x] > data[y]) swap(x, y);
        data[x] += data[y];
        data[y] = x;
        return true;
    }

    bool same(int x, int y) { return find(x) == find(y); }

    int find(int x) {
        if(data[x] < 0) return x;
        return (data[x] = find(data[x]));
    }
    int size(int x) { return -data[find(x)]; }
    const int operator[](const int k) { return find(k); }
    vector<vector<int>> belong() {
        vector<vector<int>> res(data.size());
        for(int i = 0; i < data.size(); i++) res[find(i)].emplace_back(i);
        return res;
    }
    void reset() {
        num = data.size();
        data.assign(num, -1);
    }
};

// f(x, y) := x に y の値を merge
// 参照で受け取る
template <typename T, typename F> struct UnionFind_value {
    vector<int> dat;
    int num;
    vector<T> val;
    F f;
    UnionFind_value(int n, F f, T e) : dat(n, -1), val(n, e), num(n), f(f) {}
    UnionFind_value(const vector<T> &a, F f) : dat(a.size(), -1), val(a), num(a.size()), f(f) {}

    bool unite(int x, int y) {
        x = find(x), y = find(y);
        if(x == y) return false;
        if(dat[x] > dat[y]) swap(x, y);
        dat[x] += dat[y];
        dat[y] = x;
        --num;
        f(val[x], val[y]);
        return true;
    }
    bool same(int x, int y) { return find(x) == find(y); }

    int find(int x) {
        if(dat[x] < 0) return x;
        return (dat[x] = find(dat[x]));
    }
    int size(int x) { return -dat[find(x)]; }
    const int operator[](const int k) { return find(k); }
    vector<vector<int>> belong() {
        vector<vector<int>> res(dat.size());
        for(int i = 0; i < dat.size(); i++) res[find(i)].emplace_back(i);
        return res;
    }
    T &value(int x) { return val[find(x)]; }
    void reset(T e) {
        num = dat.size();
        dat.assign(num, -1);
        val.assign(num, e);
    }
};
