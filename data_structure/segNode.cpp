namespace segnode_impl {
struct SegNodes {
  private:
    int ceil_pow2(int n) {
        int x = 0;
        while((1U << x) < (unsigned int)(n)) x++;
        return x;
    }
    int log;
    struct single {
        const int x;
        struct single_element_itr {
            int x;
            single_element_itr(const int &x) : x(x) {}
            int &operator*() { return x; }
            bool operator!=(const single_element_itr &) const { return x; }
            void operator++() { x >>= 1; }
        };
        single(int x) : x(x) {}
        single_element_itr begin() const { return single_element_itr(x); }
        single_element_itr end() const { return single_element_itr(x); }
    };

    struct Range {
        const int l, r;
        Range(int l, int r) : l(l), r(r) {}
        struct range_itr {
            int l, r;
            bool left;
            void run() {
                if(!(l < r)) return;
                while(true) {
                    if(left) {
                        if(l & 1)
                            return;
                        else
                            left = false;
                    }
                    if(!left) {
                        if(r & 1)
                            return;
                        else
                            left = true;
                    }
                    l >>= 1, r >>= 1;
                }
            }
            range_itr(int l, int r) : l(l), r(r), left(true) { run(); }
            int operator*() const { return left ? l : r - 1; }
            void operator++() {
                if(left)
                    l++;
                else
                    r--;
                run();
            }
            bool operator!=(const range_itr &) { return l < r; }
        };
        range_itr begin() const { return range_itr(l, r); }
        range_itr end() const { return range_itr(l, r); }
    };

  public:
    int n, sz;
    SegNodes(int n) : n(n) {
        sz = 1;
        log = ceil_pow2(n);
        sz = 1 << log;
    }
    // x を含むノードを返すイテレータ
    single get_nodes(int x) const { return single(x + sz); }
    Range get_nodes(int l, int r) const { return Range(l + sz, r + sz); }
    single operator()(int x) const { return get_nodes(x); }
    Range operator()(int l, int r) const { return get_nodes(l, r); }
    // ノードの個数
    int size() { return sz * 2; }
};
} // namespace segnode_impl
using segnode_impl::SegNodes;