#include "modular/modint1.cpp"
#include "my_template.hpp"

// Verify compatibility when the everyday template is included after the library.
int main() {
    mint x = 7;
    vmint values{1, 2, 3};
    Poly f = values;
    assert(eval(f, x) == mint(162));
    assert(C(5, 2) == mint(10));
    vi indices(3);
    rep(i, 3) indices[i] = i;
    assert(indices == vi({0, 1, 2}));
}
