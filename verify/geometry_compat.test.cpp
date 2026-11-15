#include "my_template.hpp"
#include <sstream>
#define min(a,b) forbidden_min_macro(a,b)
#define max(a,b) forbidden_max_macro(a,b)
#define ON 7123
#include "geometry/geometry.hpp"
using namespace cp_geometry;
#ifndef min
#error caller min macro was removed
#endif
#ifndef max
#error caller max macro was removed
#endif
static_assert(ON==7123, "caller ON macro changed");
#undef min
#undef max
#undef ON

int main() {
    std::istringstream input("7");
    auto* previous=std::cin.rdbuf(input.rdbuf());
    // This expanded to IN(value), ambiguous with the old enum constant IN.
    INT(value);
    std::cin.rdbuf(previous);
    assert(value==7);
    ConvexPolygon square({{0,0},{2,0},{2,2},{0,2}});
    assert(square.contains({1,1})==Location::inside);
    assert(square.contains({0,1})==Location::boundary);
    assert(square.contains({3,1})==Location::outside);
    std::vector<int> values{3,1,2}; std::sort(all(values)); assert(values.front()==1);
}
