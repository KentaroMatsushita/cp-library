constexpr int SIZE = 400;
template <typename T> struct Node {
    // 持つ情報を書く
    T v[SIZE];
    int l, r;
    Node(vector<T> &a, int l, int r) : l(l), r(r) { rep2(i, l, r - 1) v[i - l] = a[i]; }
    void eval() {
        // 遅延分を評価
        for(auto &e : v) e += lazy;

        // 遅延を 0 にセット
        lazy = 0;
    }
    void update(int L, int R, ll x) {
        if(R <= l or r <= L) return;
        if(L <= l and r <= R) {
            // 全体をまとめて変更する

        } else {
            L = clamp(L, l, r) - l, R = clamp(R, l, r) - l;
            eval();
            rep2(i, L, R - 1) v[i] += x;
            //
        }
    }
    ll solve(int L, int R) {
        if(R <= l or r <= L) return LLONG_MAX;
        if(L <= l and r <= R)
            return mi;
        else {
            ll res = LLONG_MAX;
            L = clamp(L, l, r) - l, R = clamp(R, l, r) - l;
            eval();
            rep2(i, L, R - 1) chmin(res, v[i]);
            return res;
        }
    }
};
