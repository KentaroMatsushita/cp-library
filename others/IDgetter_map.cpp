template <class T> struct ID {
    int n;
    vector<T> a;
    map<T, int> mp;
    ID() : mp(), a(), n(0) {}
    ID(const vector<T> &v) {
        n = v.size();
        a.resize(n);
        rep(i, n) a[i] = v[i];
        rep(i, n) mp[a[i]] = i;
    }
    int id(const T &t) {
        if(mp.count(t)) return mp[t];
        a.emplace_back(t);
        mp[t] = n++;
        return n - 1;
    }
    vector<T> get() { return a; }
};