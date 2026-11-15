#include "geometry/geometry.hpp"
#include "my_template.hpp"
#include <sstream>
using namespace cp_geometry;
int main() {
    std::istringstream input("11");
    auto* previous=std::cin.rdbuf(input.rdbuf());
    INT(value);
    std::cin.rdbuf(previous);
    assert(value==11);
    std::vector<HalfPlane> inequalities{{-1,0,0},{0,-1,0},{1,0,2},{0,1,2}};
    auto region=half_plane_intersection(inequalities);
    assert(region.kind==RegionKind::POLYGON);
    assert(std::abs(polygon_area(region.polygon)-4)<1e-10L);
}
