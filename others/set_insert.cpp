template <typename T> void INSERT(map<int> mp, T x, T y) {
    if(s.lower_bound(x)->second >= y) return;
    auto it = mp.upper_bound(x);
    while(it != s.begin()) {
        it = prev(it);
        if(it->second <= x.second)
            it = s.erase(it);
        else
            break;
    }
    mp[x] = y;
    return;
};