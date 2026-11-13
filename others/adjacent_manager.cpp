// ソートされた列 A_0, ..., に対する値
// G(F(A_0, A_1), F(A_0, A_2), ...) を管理する
// x を追加/削除 が O(log(要素数))
// 結果を出力が O(1) で出来る
// S : set にいれる型 T : 結果の型
// F : 隣接をマージ G : 総和 rG : 逆
// 重複なし
template <typename S, typename T, typename F, typename G, typename rG> struct adj_manager {
    const F f;
    const G g;
    const rG rg;
    set<S> s;
    T sum;
    // f : 隣接のマージ g : 結果のマージ rg : 逆元
    adj_manager(const F &f, const T &e, const G &g = plus<T>(), const rG &rg = minus<T>()) : f(f), sum(e), g(g), rg(rg) {}
    void add(const S &x) {
        auto [it, _] = s.emplace(x);
        if(!_) return;
        if(next(it) != s.end()) sum = g(sum, f(x, *next(it)));
        if(it != s.begin()) sum = g(sum, f(*prev(it), x));
        if(next(it) != s.end() and it != s.begin()) sum = rg(sum, f(*prev(it), *next(it)));
    }
    template <class... Args> void emplace(const Args &...args) { add(S(args...)); }
    void erase(const S &x) {
        auto it = s.find(x);
        if(it == s.end()) return;
        it = s.erase(it);
        if(it != s.end()) sum = rg(sum, f(x, *it));
        if(it != s.begin()) sum = rg(sum, f(*prev(it), x));
        if(it != s.end() and it != s.begin()) sum = g(sum, f(*prev(it), *it));
    }

    constexpr T operator()() { return sum; }
    T circle_sum() {
        if(s.empty()) return sum;
        return g(sum, f(*s.begin(), *s.rbegin()));
    }
};
template <typename S, typename T, typename F, typename G = plus<T>, typename rG = minus<T>>
auto get_adj_manager(const F &f, const T &e, const G &g = G{}, const rG &rg = rG{}) {
    return adj_manager<S, T, F, G, rG>(f, e, g, rg);
}
