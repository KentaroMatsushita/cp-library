template <typename T> int lis(const vector<T> &v, bool strict = true) {
    vector<T> res;
    for(auto &&e : v) {
        auto it = (strict ? lower_bound(res.begin(), res.end(), e) : upper_bound(res.begin(), res.end(), e));
        if(end(res) == it)
            res.emplace_back(e);
        else
            *it = e;
    }
    return res.size();
}

template <typename T> vector<int> lis_restore(const vector<T> &v, bool strict = true) {
    vector<int> res, pre(v.size(), -1);
    res.reserve(v.size());
    for(int i = 0; i < v.size(); i++) {
        auto it = (strict ? lower_bound(res.begin(), res.end(), i, [&](int x, int y) { return v[x] < v[y]; })
                          : upper_bound(res.begin(), res.end(), i, [&](int x, int y) { return v[x] < v[y]; }));
        int pos = int(it - res.begin());
        if(pos > 0) pre[i] = res[pos - 1];
        if(end(res) == it)
            res.emplace_back(i);
        else
            *it = i;
    }
    vector<int> ret;
    ret.reserve(res.size());
    int tmp = res.empty() ? -1 : res.back();
    while(tmp != -1) {
        ret.emplace_back(tmp);
        tmp = pre[tmp];
    }
    reverse(begin(ret), end(ret));
    return ret;
}
