vector<int> ciclecut(vector<pair<int, int>> &E) {
    if(E.empty()) return {};
    int M = 0;
    for(auto &&[x, y] : E) M = max({M, x, y});
    static vector<int> l, r;
    if(l.size() <= M) l.resize(M + 1), r.resize(M + 1);
    for(auto &&[x, y] : E) { l[x] = r[x] = l[y] = r[y] = -1; }
    for(auto &&[x, y] : E) {
        (l[x] == -1 ? l[x] : r[x]) = y;
        (l[y] == -1 ? l[y] : r[y]) = x;
    }
    vi res = {E[0].fi, E[0].se};
    res.reserve(E.size());
    while(res.size() < E.size()) {
        int to = l[res.back()] ^ r[res.back()] ^ res[res.size() - 2];
        res.emplace_back(to);
    }
    return res;
}
