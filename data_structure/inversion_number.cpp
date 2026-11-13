template <typename T> struct BinaryIndexedTree {
    vector<T> data;

    BinaryIndexedTree(int sz) { data.assign(++sz, 0); }

    T sum(int k) {
        T ret = 0;
        for(++k; k > 0; k -= k & -k) ret += data[k];
        return (ret);
    }

    void add(int k, T x) {
        for(++k; k < data.size(); k += k & -k) data[k] += x;
    }
};
#define BIT BinaryIndexedTree

template <typename T> long long inversion(vector<T> a, bool is_zip = true) {
    if(!is_zip) {
        vector<T> id(a.size());
        for(int i = 0; i < a.size(); ++i) id[i] = a[i];
        sort(id.begin(), id.end());
        for(int i = 0; i < a.size(); ++i) a[i] = lower_bound(id.begin(), id.end(), a[i]) - id.begin();
    }
    BIT<int> b(a.size() + 1);
    long long res = (long long)a.size() * (a.size() - 1) / 2;
    for(int i = 0; i < a.size(); ++i) {
        res -= b.sum(a[i]);
        b.add(a[i], 1);
    }
    return res;
}