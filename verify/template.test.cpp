#include "my_template.hpp"
int main(){
    assert(min(3000000000LL,2)==2 && min(3,4000000000LL)==3);
    assert(counter({}).empty());
    std::ostringstream edge_output; edge<int> e(2,3); edge_output<<e; assert(edge_output.str()=="2");
    i128 value=-(i128(1)<<126);value+=value;
    ostringstream output;output<<value;
    assert(output.str()=="-170141183460469231731687303715884105728");
    istringstream input(output.str());i128 restored;input>>restored;assert(restored==value);
    vector<int> a{-1,inf<int>,7},copy=a;
    ostringstream capture;auto* previous=cout.rdbuf(capture.rdbuf());OUT2(a);cout.rdbuf(previous);assert(a==copy);
    cout<<"common template regressions passed\n";
}
