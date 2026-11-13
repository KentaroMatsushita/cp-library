// f(v) : v から一回の操作でできる状態の vector
template <typename T, typename F> vector<T> restore_action(T init, T target, F f) {
    map<T, int> mp;
    map<T, T> pre;
    queue<T> q;
    q.emplace(init);
    while(!empty(q)) {
        auto now = q.front();
        q.pop();
        if(now == target) {
            vector<T> res;
            while(true) {
                res.emplace_back(now);
                if(now == init) {
                    REV(res);
                    return res;
                }
                now = pre[now];
            }
        }
        int c = mp[now];
        fore(nxt, f(now)) {
            if(!mp.count(nxt)) {
                mp[nxt] = c + 1;
                pre[nxt] = now;
                q.emplace(nxt);
            }
        }
    }
    return {};
}