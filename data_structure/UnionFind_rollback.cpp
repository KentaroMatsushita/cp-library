struct UnionFind_undo {
    vector<int> data;
    vector<int> a;
    int num;

    UnionFind_undo(int n) : num(n) {
        data.assign(n, -1);
        a.assign(n, 0);
    }
    UnionFind_undo(vi _a) : a(_a), num(si(_a)) { data.assign(si(_a), -1); }

    stack<tuple<int, int, int, int, int>> st;
    bool unite(int x, int y) {
        x = find(x), y = find(y);
        if(x == y) {
            st.emplace(-1, -1, -1, -1, -1);
            return false;
        }
        num--;
        if(data[x] > data[y]) swap(x, y);
        st.emplace(x, y, data[y], a[x], a[y]);
        chmax(a[x], a[y]);
        data[x] += data[y];
        data[y] = x;
        return true;
    }

    void rollback(int num = 1) {
        while(num--) {
            auto [x, y, dy, ax, ay] = st.top();
            st.pop();
            if(x == -1) continue;
            ++this->num;
            data[x] -= dy;
            data[y] = dy;
            a[x] = ax, a[y] = ay;
        }
    }

    bool same(int x, int y) { return find(x) == find(y); }

    int find(int x) {
        if(data[x] < 0) return x;
        return find(data[x]);
    }
    int size(int x) { return -data[find(x)]; }
};
