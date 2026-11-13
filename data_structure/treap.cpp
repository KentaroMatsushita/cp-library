template <class Monoid = int, class OperatorMonoid = Monoid> struct Treap {
    using F = function<Monoid(Monoid, OperatorMonoid)>;                 // 作用素を適用
    using G = function<OperatorMonoid(OperatorMonoid, OperatorMonoid)>; // 作用素同士のマージ
    const Monoid M1;
    const OperatorMonoid OM0;
    const F f;
    const G g;
    struct node {
        int pri;
        Monoid key;
        OperatorMonoid lazy; // 遅延分
        node *l, *r;
        int size;
        node() = default;
        node(Monoid key, OperatorMonoid OM) : key(key), pri(rnd(INT_MAX)), size(1), l(nullptr), r(nullptr), lazy(OM) {}
    };

    vector<node> nodes;
    int ptr;

    inline node *allocate(Monoid key) {
        if(ptr >= int(nodes.size())) throw std::out_of_range("Treap node pool exhausted");
        return &(nodes[ptr++] = node(key, OM0));
    }

    Treap(int n, const F &f, const G &g, const Monoid M1, const OperatorMonoid OM0) : f(f), g(g), M1(M1), OM0(OM0), ptr(0), nodes(n) {}

    node *root = nullptr;
    int size(node *t) { return t != nullptr ? t->size : 0; }

    void propagate(node *t) {
        if(t == nullptr) return;
        t->key = f(t->key, t->lazy);
        if(t->l) { t->l->lazy = g(t->l->lazy, t->lazy); }
        if(t->r) { t->r->lazy = g(t->r->lazy, t->lazy); }
        t->lazy = OM0;
    }

    node *refresh(node *t) {
        propagate(t);
        if(t == nullptr) return t;
        t->size = 1 + size(t->l) + size(t->r);
        return t;
    }

    pair<node *, node *> split(node *&t, Monoid k) {
        propagate(t);
        if(!t)
            return {t, t};
        else if(t->key < k) {
            auto p = split(t->r, k);
            t->r = p.fi;
            return {refresh(t), p.se};
        } else {
            auto p = split(t->l, k);
            t->l = p.se;
            return {p.fi, refresh(t)};
        }
    }

    node *merge(node *a, node *b) {
        propagate(a);
        propagate(b);
        if(!a or !b) return a ? a : b;
        if(a->pri < b->pri) {
            a->r = merge(a->r, b);
            return refresh(a);
        } else {
            b->l = merge(a, b->l);
            return refresh(b);
        }
    }

    node *count(node *t, Monoid k) {
        propagate(t);
        if(!t) return nullptr;
        if(k < t->key) return count(t->l, k);
        if(k == t->key) return t;
        return count(t->r, k);
    }

    node *nth(node *t, int n) {
        propagate(t);
        if(!t) return nullptr;
        if(n <= size(t->l)) return nth(t->l, n);
        if(n == size(t->l) + 1) return t;
        return nth(t->r, n - size(t->l) - 1);
    }

    void del(node *&t) {
        propagate(t);
        if(!t) return;
        if(t->l) del(t->l);
        if(t->r) del(t->r);
        t = nullptr;
    }

    // 一個ずつmergeしていく
    void send(node *&from, node *&to, Monoid k) {
        propagate(from);
        if(!from) return;
        if(from->l) send(from->l, to, k);
        if(from->r) send(from->r, to, k);
        Monoid M = from->key;
        get<0>(M) -= get<0>(k);
        auto [a, b] = split(to, M);
        get<1>(M)++;
        to = merge(merge(a, allocate(M)), b);
        from = nullptr;
    }

    bool count(Monoid k) { return count(root, k) != nullptr; }

    node *insert(Monoid k) {
        auto [l, r] = split(root, k);
        auto c = allocate(k);
        return root = merge(merge(l, c), r);
    }

    // [k,k2) を削除
    bool erase(Monoid k, Monoid k2) {
        auto [a,b] = split(root,k);
        auto [c,d] = split(b,k2);
        bool removed = c != nullptr;
        root = merge(a, d);
        return removed;
    }
    Monoid nth(int n) {
        n++;
        node *ans = nth(root, n);
        return ans ? ans->key : M1;
    }
    void solve(int k) {
        Monoid K(k, 0, 0);
        auto [a, b] = split(root, K);
        auto [c, d] = split(b, Monoid(k * 2, 0, 0));
        send(c, a, K);
        if(d) { d->lazy = OperatorMonoid(k, 1); }
        root = merge(a, d);
    }
};
