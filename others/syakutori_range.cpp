// a[i] に対して a[i] + b[j] \in [L, R) を満たす j の範囲を返す
// a と b はソート済みであること
template <class T> vector<pair<int, int>> syakutoriRange(const vector<T> &a, const vector<T> &b, T L, T R) {
    int n = a.size();
    vector<pair<int, int>> res(n);
    int l = b.size(), r = b.size();
    for(int i = 0; i < n; i++) {
        while(l and b[l - 1] + a[i] >= L) l--;
        while(r and b[r - 1] + a[i] >= R) r--;
        res[i] = pair(l, r);
    }
    return res;
}