template <typename T, typename Container = vector<T>> struct imos_cyclic {
    int n;
    Container a;
    imos_cyclic(int n) : n(n) {
        if constexpr(is_same_v<Container, vector<T>>) { a.resize(n + 1); }
    }
    template <typename S> void add(S l, S r, T val) {
        a[0] -= (l / n + 1) * val;
        a[0] += (r / n + 1) * val;
        a[l % n] += val;
        a[r % n] -= val;
    }

    void build() {
        assert(!exchange(is_built, true));
        if constexpr(is_same_v<Container, vector<T>>) {
            for(int i = 0; i < n; i++) a[i + 1] += a[i];
        } else {
            T cum = 0;
            for(auto &[x, val] : a) val += exchange(cum, cum + val);
        }
    }
    Container get_vec() {
        assert(is_built);
        return a;
    }

  private:
    bool is_built = false;
};