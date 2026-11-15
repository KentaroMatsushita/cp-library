"""Explicit publication boundary shared by site generation and verification."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
DIRECTORIES = ('Strings', 'data_structure', 'dp', 'flow', 'fraction', 'graph',
               'linear_algebra', 'math', 'modular', 'monoid', 'numerical',
               'others', 'tree', 'fastIO', 'geometry', '.vscode/geometry')
ROOT_FILES = ('my_template.hpp', 'cht.cpp', 'monge_rowmin.cpp', 'monotone-minima.cpp',
              'Wavelet_Matrix_xor.cpp')
SUPPORT_FILES = ('ac-library/atcoder/maxflow', 'ac-library/atcoder/maxflow.hpp',
                 'ac-library/atcoder/internal_queue', 'ac-library/atcoder/internal_queue.hpp',
                 'ac-library/atcoder/LICENSE',
                 'geometry/detail/macro_push.hpp', 'geometry/detail/macro_pop.hpp')

def sources():
    result = []
    for directory in DIRECTORIES:
        for p in (ROOT / directory).glob('*'):
            if p.is_file() and not p.is_symlink() and p.suffix in ('.cpp', '.hpp', '.h'):
                result.append(p.relative_to(ROOT).as_posix())
    result += [p for p in ROOT_FILES if (ROOT / p).is_file()]
    return sorted(set(result))

def catalog():
    result=json.loads((ROOT / 'docs/catalog.json').read_text())
    for notes in sorted((ROOT / 'docs').glob('notes*.json')):
        for path, values in json.loads(notes.read_text()).items():
            result[path]=dict(zip(('title','summary','complexity','contract','usage'),values))
            result[path]['tests']=[]
    for path, changes in json.loads((ROOT/'docs/overrides.json').read_text()).items():
        result[path].update(changes)
    for test, paths in json.loads((ROOT/'docs/coverage.json').read_text()).items():
        for path in paths:
            assert path in result, path
            if test not in result[path]['tests']: result[path]['tests'].append(test)
    return result

def fingerprint():
    paths = sources() + list(SUPPORT_FILES) + [p.relative_to(ROOT).as_posix() for p in (ROOT / 'verify').glob('*.test.cpp')]
    paths += ['tools/test.py', 'tools/library.py', 'tools/check_headers.py']
    h = hashlib.sha256()
    for path in sorted(paths):
        h.update(path.encode() + b'\0' + (ROOT / path).read_bytes() + b'\0')
    return h.hexdigest()
