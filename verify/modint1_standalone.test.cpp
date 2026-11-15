#include "modular/modint1.cpp"

// This include is intentionally first: the library must supply its own headers.
#include <cassert>
#include <random>
#include <sstream>

#ifdef mint
#error "modint1 must not define a mint macro"
#endif
#ifdef vmint
#error "modint1 must not define a vmint macro"
#endif

constexpr modular::mint constant_value(1000000009LL);
static_assert(constant_value.a == 2);

int main() {
    modular::mint x = -7, y = 3;
    assert((x + y).a == modular::MOD - 4);
    assert((x * y / y) == x);
    assert(modular::modint<17>(5) / modular::modint<17>(3) == modular::modint<17>(13));
    assert(modular::C(20, 7).a == 77520);
    assert(modular::P(7, 3).a == 210);

    std::mt19937 rng(238);
    for(int n : {0, 1, 3, 99, 100, 137}) {
        modular::Poly a(n), b(n + 3);
        for(auto &value : a) value = rng();
        for(auto &value : b) value = rng();
        const auto product = a * b;
        modular::Poly expected(n ? a.size() + b.size() - 1 : 0);
        for(std::size_t i = 0; i < a.size(); ++i)
            for(std::size_t j = 0; j < b.size(); ++j) expected[i + j] += a[i] * b[j];
        assert(product == expected);
    }

    modular::Poly f{2, 3, 4};
    assert(modular::TaylorShift(f, 2) == modular::Poly({24, 19, 4}));
    assert(modular::eval(f, 2) == modular::mint(24));
    std::istringstream input("-1000000008");
    modular::mint parsed;
    input >> parsed;
    assert(parsed.a == modular::MOD - 1);
    std::ostringstream output;
    output << f;
    assert(output.str() == "2 3 4 ");
}
