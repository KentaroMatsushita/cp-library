struct UnionFind {
    vector<int> data, col;
    int num;

    UnionFind(int n) : num(n), col(n) { data.assign(n, -1); }

    stack<tuple<int, int, int>> st;
    bool unite(int x, int y) {
        int cx = color(x), cy = color(y);
        x = find(x), y = find(y);
        if(x == y) {
            st.emplace(-1, -1, -1);
            st.emplace(-1, -1, -1);
            if(cx == cy) return false;
            return true;
        }
        num--;
        if(data[x] > data[y]) swap(x, y);
        data[x] += data[y];
        st.emplace(y, data[y], 0);
        data[y] = x;
        if(cx == cy) {
            col[y] ^= 1;
            st.emplace(y, 1, 1);
        } else
            st.emplace(y, 0, 1);
        return true;
    }

    void rollback(int num = 1) {
        num *= 2;
        while(num--) {
            auto [a, b, c] = st.top();
            st.pop();
            if(a == -1) continue;
            if(c) {
                col[a] ^= b;
            } else {
                ++this->num;
                int x = data[a];
                data[a] = b;
                data[x] -= data[a];
            }
        }
    }

    bool same(int x, int y) { return find(x) == find(y); }

    int color(int x) {
        int c = 0;
        while(data[x] >= 0) {
            c ^= col[x];
            x = data[x];
        }
        return c;
    }
    int find(int x) {
        if(data[x] < 0) return x;
        return find(data[x]);
    }
    int size(int x) { return -data[find(x)]; }
};
