template <typename T> struct Parallel_binary_search {
    // クエリの個数
    int Q;
    vector<T> ok, ng, mid;
    vector<int> idx, nidx;
    Parallel_binary_search(int Q, T OK, T NG) : Q(Q), ok(Q, OK), ng(Q, NG), mid(Q), idx(Q), nidx(Q) { iota(begin(idx), end(idx), 0); }

    // v は イベントの vector
    // T との大小関係が定義されている必要がある
    // v < T(x) は v のイベントは x より先に処理していいか
    // p(event) はイベントを処理する関数
    // f(i, x) は今 i 番目のクエリが x で行けるかを判定する
    // "x 未満" のイベントは全て処理されている
    // init は二分探索の開始毎に呼び出される
    template <typename EVENT, typename P, typename F, typename INIT>
    vector<T> build(const vector<EVENT> &v, const P &p, const F &f, const INIT &init, bool greater = false) {
        while(true) {
            bool done = true;
            for(int i = 0; i < Q; i++)
                if(abs(ok[i] - ng[i]) > 1) done = false;
            if(done) return ok;
            for(int i = 0; i < Q; i++) mid[i] = (ok[i] + ng[i]) >> 1;

            init();
            int pos = 0;
            for(int i = 0; i < Q; i++) {
                while(pos < v.size() and v[pos] < mid[idx[i]]) { p(v[pos++]); }
                int j = i;
                while(j < Q and ok[idx[i]] == ok[idx[j]] and ng[idx[i]] == ng[idx[j]]) j++;
                int l = i, r = j;
                for(int k = i; k < j; k++) {
                    bool can = f(idx[k], mid[idx[k]]);
                    if(can ^ (ng[idx[k]] < ok[idx[k]]) ^ greater)
                        nidx[--r] = idx[k];
                    else
                        nidx[l++] = idx[k];
                    (can ? ok[idx[k]] : ng[idx[k]]) = mid[idx[k]];
                }
                bool flag = true;
                i = j - 1;
            }
            swap(idx, nidx);
        }
    }
};

// https://atcoder.jp/contests/joi2012ho/submissions/27295749
