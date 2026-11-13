// add(s, t) : x 座標が s 以上 or y 座標が t 以下という条件を追加
// min() : 上の条件を満たす (x, y) 全体での x + y の最小値

template <typename T> struct RectangleSets {
    set<pair<T, T>> s;
    multiset<T> ans;
    RectangleSets() {
        const T inf = numeric_limits<T>::max() / 3;
        ans.emplace(0);
        s.emplace(0, inf);
        s.emplace(inf, 0);
    }
    void add(T x, T y) {
        auto rp = make_pair(x, y);
        auto it = s.lower_bound(rp);
        if(y <= it->se) return;
        bool flag = true;
        while(1) {
            it = prev(it);
            if(it->se <= y) {
                auto p = prev(it), q = next(it);
                T x = p->fi + it->se, y = it->fi + q->se;
                if(flag) ans.erase(ans.lower_bound(y));
                ans.erase(ans.lower_bound(x));
                it = s.erase(it);
                flag = false;
            } else
                break;
        }
        auto [i, _] = s.emplace(x, y);
        auto p = prev(i), q = next(i);
        ans.emplace(p->fi + y);
        ans.emplace(x + q->se);
        if(flag) ans.erase(ans.lower_bound(p->fi + q->se));
    }
    T min() { return *begin(ans); }
};
