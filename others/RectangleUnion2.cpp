// [0, a] * [0, b] の長方形の union の面積を求める
template <typename T = long long> struct RectangleUnion {
    map<T, T> mp;
    T s;
    RectangleUnion() : s(0), mp({{0, numeric_limits<T>::max()}, {numeric_limits<T>::max(), 0}}) {}
    void add(T a, T b) {

        auto it = mp.upper_bound(a);
        if(it->se >= b) return;
        T d = it->second;
        while(1) {
            it = prev(it);
            if(it->se > b) {
                s += (a - it->fi) * (b - d);
                mp.emplace_hint(it, a, b);
                break;
            } else {
                s -= (it->fi - prev(it)->fi) * (it->se - d);
                it = mp.erase(it);
            }
        }
    }
    T get() { return s; }
};