#pragma once
#include "monotone-minima.hpp"

// a : 上凸
template <typename T> vector<T> convex_max_conv(const vector<T> &a, const vector<T> &b) {
    if(empty(a) or empty(b)) return {};
    auto arg = monotone_minima(a.size() + b.size() - 1, b.size(), [&](int i, int j, int k) {
        if(i < k) return false;
        if(i - j >= a.size()) return true;
        return a[i - j] + b[j] <= a[i - k] + b[k];
    });
    vector<T> res(a.size() + b.size() - 1);
    rep(i, res.size()) res[i] = a[i - arg[i]] + b[arg[i]];
    return res;
}

// a : 下凸
template <typename T> vector<T> convex_min_conv(const vector<T> &a, const vector<T> &b) {
    if(empty(a) or empty(b)) return {};
    auto arg = monotone_minima(a.size() + b.size() - 1, b.size(), [&](int i, int j, int k) {
        if(i < k) return false;
        if(i - j >= a.size()) return true;
        return a[i - j] + b[j] >= a[i - k] + b[k];
    });
    vector<T> res(a.size() + b.size() - 1);
    rep(i, res.size()) res[i] = a[i - arg[i]] + b[arg[i]];
    return res;
}
