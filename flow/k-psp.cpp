namespace suisen {
template <typename Cost> struct ProjectSelection {
    using Item = std::size_t;

  private:
    template <typename T> using CostFunc = std::array<T, 2>;

  public:
    using CostFunc0 = Cost;
    using CostFunc1 = CostFunc<CostFunc0>;
    using CostFunc2 = CostFunc<CostFunc1>;
    using CostFunc3 = CostFunc<CostFunc2>;

    ProjectSelection() = default;

    explicit ProjectSelection(std::size_t n) : _num(n), _cost1(n) {}

    void add_cost(CostFunc0 cost) { _cost0 += cost; }
    void add_profit(CostFunc0 profit) { add_cost(-profit); }

    void add_cost_0(Item i, Cost cost) { add_cost(i, CostFunc1{cost, 0}); }
    void add_cost_1(Item i, Cost cost) { add_cost(i, CostFunc1{0, cost}); }
    void add_profit_0(Item i, Cost profit) { add_cost(i, CostFunc1{-profit, 0}); }
    void add_profit_1(Item i, Cost profit) { add_cost(i, CostFunc1{0, -profit}); }

    void add_cost(Item i, CostFunc1 cost) {
        assert(i < _num);
        _cost1[i][0] += cost[0];
        _cost1[i][1] += cost[1];
    }
    void add_profit(Item i, CostFunc1 profit) { add_cost(i, neg(profit)); }

    void add_cost_01(Item i, Item j, Cost cost) {
        assert(i < _num);
        assert(j < _num);
        assert(i != j);
        add_edge(i, j, cost);
    }

    void add_cost_10(Item i, Item j, Cost cost) { add_cost_01(j, i, cost); }
    void add_cost_not_same(Item i, Item j, Cost cost) { add_cost(i, j, CostFunc2{0, cost, cost, 0}); }
    void add_profit_00(Item i, Item j, Cost profit) { add_cost(i, j, CostFunc2{-profit, 0, 0, 0}); }
    void add_profit_11(Item i, Item j, Cost profit) { add_cost(i, j, CostFunc2{0, 0, 0, -profit}); }
    void add_profit_same(Item i, Item j, Cost profit) { add_cost(i, j, CostFunc2{-profit, 0, 0, -profit}); }

    // cost: submodular (i.e., cost[0][1] + cost[1][0] >= cost[0][0] + cost[1][1])
    void add_cost(Item i, Item j, CostFunc2 cost) {
        assert(i < _num);
        assert(j < _num);
        assert(i != j);
        add_cost(cost[0][0]);
        add_cost(i, CostFunc1{0, cost[1][0] - cost[0][0]});
        add_cost(j, CostFunc1{0, cost[1][1] - cost[1][0]});
        add_cost_01(i, j, (cost[0][1] + cost[1][0]) - (cost[0][0] + cost[1][1]));
    }
    void add_profit(Item i, Item j, const CostFunc2 &profit) { add_cost(i, j, neg(profit)); }

    // cost: submodular (i.e., cost(X) + cost(Y) >= cost(X | Y) + cost(X & Y))
    void add_cost(Item i, Item j, Item k, CostFunc3 cost) {
        assert(i < _num);
        assert(j < _num);
        assert(k < _num);
        assert(i != j);
        assert(j != k);
        assert(k != i);
        const Cost A = cost[0][0][0], B = cost[0][0][1];
        const Cost C = cost[0][1][0], D = cost[0][1][1];
        const Cost E = cost[1][0][0], F = cost[1][0][1];
        const Cost G = cost[1][1][0], H = cost[1][1][1];
        const Cost P = (A + D + F + G) - (B + C + E + H);
        if(P >= 0) {
            const Cost P1 = F - B, P2 = G - E, P3 = D - C;
            const Cost P12 = (C + E) - (A + G), P23 = (B + C) - (A + D), P31 = (B + E) - (A + F);
            add_cost(A);
            add_cost(i, CostFunc1{0, P1});
            add_cost(j, CostFunc1{0, P2});
            add_cost(k, CostFunc1{0, P3});
            add_cost_01(i, j, P12);
            add_cost_01(j, k, P23);
            add_cost_01(k, i, P31);
            add_profit_all_1(std::array<Item, 3>{i, j, k}, P);
        } else {
            const Cost P1 = C - G, P2 = B - D, P3 = E - F;
            const Cost P21 = (D + F) - (B + H), P32 = (F + G) - (E + H), P13 = (D + G) - (C + H);
            add_cost(H);
            add_cost(i, CostFunc1{P1, 0});
            add_cost(j, CostFunc1{P2, 0});
            add_cost(k, CostFunc1{P3, 0});
            add_cost_01(j, i, P21);
            add_cost_01(k, j, P32);
            add_cost_01(i, k, P13);
            add_profit_all_0(std::array<Item, 3>{i, j, k}, -P);
        }
    }

    template <typename Container> auto add_profit_all_0(const Container &is, Cost profit) -> decltype(static_cast<Item>(*is.begin()), is.end(), void()) {
        assert(profit >= 0);
        if(is.size() == 0) {
            return add_profit(profit);
        } else if(is.size() == 1) {
            return add_profit(is[0], CostFunc1{profit, 0});
        } else if(is.size() == 2) {
            return add_profit_00(is[0], is[1], profit);
        }
        add_profit(profit);
        const Item aux = _num + _num_aux++;
        add_edge(S, aux, profit);
        for(Item i : is) { add_edge(aux, i, profit); }
    }
    template <typename Container> auto add_profit_all_1(const Container &is, Cost profit) -> decltype(static_cast<Item>(*is.begin()), is.end(), void()) {
        assert(profit >= 0);
        if(is.size() == 0) {
            return add_profit(profit);
        } else if(is.size() == 1) {
            return add_profit(is[0], CostFunc1{0, profit});
        } else if(is.size() == 2) {
            return add_profit_11(is[0], is[1], profit);
        }
        add_profit(profit);
        const Item aux = _num + _num_aux++;
        for(Item i : is) { add_edge(i, aux, profit); }
        add_edge(aux, T, profit);
    }

    auto min_cost() {
        atcoder::mf_graph<Cost> g(_num + _num_aux + 2);
        const std::size_t s = _num + _num_aux, t = s + 1;
        make_edges_for_cost1();
        for(const auto &[i, j, cost] : _edges) {
            std::size_t u = (i == S) ? s : (i == T) ? t : i;
            std::size_t v = (j == S) ? s : (j == T) ? t : j;
            g.add_edge(u, v, cost);
        }
        Cost ans = _cost0 + g.flow(s, t);
        auto x = g.min_cut(s);
        x.resize(_num);
        for(std::size_t i = 0; i < _num; ++i) { x[i] = not x[i]; }
        return std::make_pair(ans, x);
    }
    auto max_profit() {
        auto res = min_cost();
        res.first = -res.first;
        return res;
    }

  private:
    static constexpr std::size_t S = ~0;
    static constexpr std::size_t T = ~1;

    std::size_t _num = 0;
    std::size_t _num_aux = 0;
    Cost _cost0 = 0;
    std::vector<CostFunc1> _cost1;
    std::vector<std::tuple<Item, Item, Cost>> _edges{};

    void make_edges_for_cost1() {
        for(std::size_t i = 0; i < _num; ++i) {
            CostFunc1 &cost = _cost1[i];
            if(cost[0] <= cost[1]) {
                add_cost(cost[0]);
                add_edge(S, i, cost[1] - cost[0]);
            } else {
                add_cost(cost[1]);
                add_edge(i, T, cost[0] - cost[1]);
            }
            cost = {0, 0};
        }
    }
    void add_edge(std::size_t i, std::size_t j, Cost cost) {
        assert(cost >= 0);
        assert(i == S or i == T or i < _num + _num_aux);
        assert(j == S or j == T or j < _num + _num_aux);
        assert(i != j);
        if(cost == 0) return;
        _edges.emplace_back(i, j, cost);
    }

    static constexpr Cost neg(Cost cost) { return -cost; }
    template <typename T> static constexpr CostFunc<T> neg(const CostFunc<T> &cost) { return {neg(cost[0]), neg(cost[1])}; }
};
} // namespace suisen

namespace suisen {
template <typename Cost> struct KProjectSelection {
    using Item = std::size_t;
    using CostFunc0 = Cost;
    using CostFunc1 = std::vector<CostFunc0>;
    using CostFunc2 = std::vector<CostFunc1>;

    KProjectSelection() = default;

    template <typename Container, std::enable_if_t<std::is_integral_v<typename Container::value_type>, std::nullptr_t> = nullptr>
    explicit KProjectSelection(const Container &ks) : _num(ks.size()), _ks(ks.begin(), ks.end()), _x(_num) {
        std::size_t id = 0;
        for(std::size_t i = 0; i < _num; ++i) {
            assert(_ks[i] >= 1);
            _x[i].resize(_ks[i]);
            _x[i][0] = ~0;
            for(std::size_t d = 1; d < _ks[i]; ++d) { _x[i][d] = id++; }
        }
        // [x_i < d] <===> _x[i][d]
        // [x_i = d] <===> ~_x[i][1] & ~_x[i][2] & ... & ~_x[i][d] & _x[i][d+1] & _x[i][d+2] & ... & _x[i][k_i]
        _psp = ProjectSelection<Cost>(id);
        for(std::size_t i = 0; i < _num; ++i) {
            for(std::size_t d = 1; d < _ks[i] - 1; ++d) { _psp.add_cost_10(_x[i][d], _x[i][d + 1], INF); }
        }
    }
    KProjectSelection(std::size_t n, std::size_t k) : KProjectSelection(std::vector<std::size_t>(n, k)) {}

    void add_cost(CostFunc0 cost) { _psp.add_cost(cost); }
    void add_profit(CostFunc0 profit) { _psp.add_profit(profit); }

    void add_cost(Item i, const CostFunc1 &cost) {
        assert(i < _num);
        assert(cost.size() == _ks[i]);
        _psp.add_cost(cost[_ks[i] - 1]);
        for(std::size_t d = 1; d < _ks[i]; ++d) { _psp.add_cost(_x[i][d], {0, cost[d - 1] - cost[d]}); }
    }
    void add_profit(Item i, CostFunc1 profit) {
        for(auto &p : profit) p = -p;
        add_cost(i, profit);
    }

    // cost: Monge (i.e., cost[i][j]+cost[i+1][j+1] <= cost[i+1][j]+cost[i][j+1])
    void add_cost(Item i, Item j, CostFunc2 cost) {
        assert(i < _num);
        assert(j < _num);
        assert(i != j);

        assert(cost.size() == _ks[i]);

        CostFunc1 cost_i(_ks[i]), cost_j(_ks[j]);
        for(std::size_t di = 0; di < _ks[i]; ++di) {
            assert(cost[di].size() == _ks[j]);
            cost_i[di] = cost[di][0];
            for(std::size_t dj = 0; dj < _ks[j]; ++dj) { cost[di][dj] -= cost_i[di]; }
        }
        for(std::size_t dj = 0; dj < _ks[j]; ++dj) {
            cost_j[dj] = cost[0][dj];
            for(std::size_t di = 0; di < _ks[i]; ++di) { cost[di][dj] -= cost_j[dj]; }
        }
        add_cost(i, cost_i);
        add_cost(j, cost_j);

        for(std::size_t di = 0; di < _ks[i]; ++di) assert(cost[di][0] == 0);
        for(std::size_t dj = 0; dj < _ks[j]; ++dj) assert(cost[0][dj] == 0);

        for(std::size_t di = 1; di < _ks[i]; ++di) {
            for(std::size_t dj = 1; dj < _ks[j]; ++dj) {
                Cost cost_00 = cost[di][dj] - cost[di][dj - 1] - cost[di - 1][dj] + cost[di - 1][dj - 1];
                // Monge
                assert(cost_00 <= 0);
                _psp.add_profit_00(_x[i][di], _x[j][dj], -cost_00);
            }
        }
    }

    std::pair<Cost, std::vector<int>> min_cost() {
        auto [ans, x_bin] = _psp.min_cost();
        std::vector<int> x(_num);
        for(std::size_t i = 0; i < _num; ++i)
            for(std::size_t di = 1; di < _ks[i]; ++di) { x[i] += not x_bin[_x[i][di]]; }
        return {ans, x};
    }

    std::pair<Cost, std::vector<int>> max_profit() {
        auto res = min_cost();
        res.first = -res.first;
        return res;
    }

  private:
    static constexpr Cost INF = std::numeric_limits<Cost>::max() / 2;

    std::size_t _num;
    std::vector<std::size_t> _ks;
    std::vector<std::vector<std::size_t>> _x;
    ProjectSelection<Cost> _psp;
};
} // namespace suisen
