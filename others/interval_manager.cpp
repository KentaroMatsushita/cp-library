// S : 持つデータの型 T : 範囲の型
template <class S, class T = int> struct IntervalManager {
    struct node {
        T l, r;
        S x;
        node(const T &l, const T &r, const S &x) : l(l), r(r), x(x) {}
        constexpr bool operator<(const node &rhs) const {
            if(l != rhs.l) return l < rhs.l;
            return r < rhs.r;
        }
    };
    const S unit;
    set<node> s;
    IntervalManager(const S &u = S()) : unit(u) {}
    IntervalManager(const vector<S> &a, const S &u = S()) : unit(u) {
        vector<node> setter;
        for(int l = 0; l < a.size(); l++) {
            int r = l;
            for(; r < a.size() and a[l] == a[r]; r++) {}
            setter.emplace_back(l, r, a[l]);
            l = r - 1;
        }
        s = set<node>(all(setter));
    }

    // x を含んだセグメントのイテレータを返す
    typename set<node>::const_iterator getIt(const T &x) const {
        auto it = s.upper_bound(node(x, numeric_limits<T>::max(), unit));
        if(it == begin(s)) return end(s);
        it = prev(it);
        if(it->l <= x and x < it->r) return it;
        return end(s);
    }

    // x を含んだセグメントの情報を持ってくる
    node getSeg(const T &x) {
        auto it = getIt(x);
        if(it != end(s)) return *it;
        return node(x, x + 1, unit);
    }

    // x 以上をはじめて含むセグメントの iterator が帰ってくる
    typename set<node>::iterator nextIt(const T &x) {
        auto it = s.upper_bound(node(x, numeric_limits<T>::max(), unit));
        if(it == begin(s)) return it;
        auto previous = prev(it);
        return previous->r > x ? previous : it;
    }

    // x に設定されてる値を取得
    S get(const T &x) const {
        auto it = getIt(x);
        if(it != end(s)) return it->x;
        return unit;
    }

    S operator[](T k) const { return get(k); }

    // [l, r) := x を set するときに消える区間加える区間ごとに関数を呼び出す
    // ただし、内包, 結合される [L, r) := x の区間も一旦消す
    template <typename ADD, typename DEL> void update(T l, T r, const S &x, const ADD &add, const DEL &del) {
        if(l >= r) return;
        auto it = s.lower_bound(node{l, 0, x});
        while(it != end(s) and it->l <= r) {
            if(it->l == r) {
                if(it->x == x) {
                    del(r, it->r, x);
                    r = it->r, it = s.erase(it);
                }
                break;
            }
            if(it->r <= r) {
                del(it->l, it->r, it->x);
                it = s.erase(it);
            } else {
                if(it->x == x) {
                    r = it->r;
                    del(it->l, it->r, it->x);
                    it = s.erase(it);
                } else {
                    del(it->l, it->r, it->x);
                    add(r, it->r, it->x);
                    node n = *it;
                    it = s.erase(it);
                    it = s.emplace_hint(it, r, n.r, n.x);
                }
            }
        }

        if(it != begin(s)) {
            it = prev(it);
            if(it->r == l) {
                if(it->x == x) {
                    del(it->l, it->r, it->x);
                    l = it->l;
                    it = s.erase(it);
                }
            } else if(l < it->r) {
                if(it->x == x) {
                    del(it->l, it->r, it->x);
                    l = min(l, it->l);
                    r = max(r, it->r);
                    it = s.erase(it);
                } else {
                    if(r < it->r) {
                        it = s.emplace_hint(next(it), r, it->r, it->x);
                        it = prev(it);
                    }
                    del(it->l, it->r, it->x);
                    if(it->l < l) add(it->l, l, it->x);
                    if(r < it->r) add(r, it->r, it->x);
                    node n = *it;
                    it = s.erase(it);
                    it = s.emplace_hint(it, n.l, l, n.x);
                }
            }
        }
        if(it != end(s)) it = next(it);
        add(l, r, x);
        s.emplace_hint(it, l, r, x);
    }

    void update(T l, T r, const S &x) {
        update(l, r, x, [](T, T, S) {}, [](T, T, S) {});
    }

    void show() {
        for(auto e : s) {
            cerr << "("
                 << "[" << e.l << ", " << e.r << "): " << e.x << ") ";
        }
        cerr << endl;
    }
};
