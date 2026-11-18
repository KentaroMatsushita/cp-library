#!/usr/bin/env python3
"""Compile reusable components independently and link multiple clients."""
import os
from pathlib import Path
import subprocess
import tempfile

ROOT=Path(__file__).resolve().parents[1]
headers=['geometry/base.hpp','geometry/convex.hpp','geometry/halfplane.hpp',
         'geometry/circle.hpp','geometry/algorithms.hpp','geometry/delaunay.hpp',
         'geometry/voronoi.hpp','geometry/power_diagram.hpp','geometry/geometry.hpp',
         '.vscode/geometry/geometry.hpp','.vscode/geometry/Points.cpp',
         'tree/HLD.cpp','tree/tree_compressing.cpp','tree/tree_util.hpp','data_structure/Wavelet_Matrix.cpp',
         'modular/modint1.cpp']
compiler=os.environ.get('CXX','g++')
flags=['-std=c++17','-iquote',str(ROOT),'-Wall','-Wextra','-Wno-sign-compare']
with tempfile.TemporaryDirectory(prefix='cp-headers-') as folder:
    temporary=Path(folder)
    for header in headers:
        source=temporary/'header.cpp'
        source.write_text(f'#include "{header}"\nint main() {{}}\n')
        subprocess.run([compiler,*flags,'-fsyntax-only',str(source)],check=True)
        print('Standalone:',header,flush=True)
    # Detect non-inline functions / header globals that break ordinary linking.
    source_a=temporary/'a.cpp'
    source_b=temporary/'b.cpp'
    source_a.write_text('#include "geometry/geometry.hpp"\nlong double area() {return cp_geometry::circle_area({{0,0},1});}\n')
    source_b.write_text('#include "geometry/geometry.hpp"\nlong double area(); int main() {return area()>0?0:1;}\n')
    output=temporary/'linked'
    subprocess.run([compiler,*flags,str(source_a),str(source_b),'-o',str(output)],check=True)
    subprocess.run([str(output)],check=True)
    # Competition submissions also paste the source instead of including it.
    pasted=temporary/'pasted.cpp'
    pasted.write_text((ROOT/'modular/modint1.cpp').read_text()+
        '\nint main() { modular::mint a=2; return (a*a).val()==4 ? 0 : 1; }\n')
    subprocess.run([compiler,*flags,str(pasted),'-o',str(output)],check=True)
    subprocess.run([str(output)],check=True)
    source_a.write_text('#include "modular/modint1.cpp"\nmodular::mint choose_two() {return modular::C(5,2);}\n')
    source_b.write_text('#include "modular/modint1.cpp"\nmodular::mint choose_two(); int main() {return choose_two()==10?0:1;}\n')
    subprocess.run([compiler,*flags,str(source_a),str(source_b),'-o',str(output)],check=True)
    subprocess.run([str(output)],check=True)
print('Standalone headers and multiple translation units passed.')
