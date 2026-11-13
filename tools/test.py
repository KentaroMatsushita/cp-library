#!/usr/bin/env python3
"""Run deterministic offline verification; failures propagate to CI."""
import argparse
import os
from pathlib import Path
import subprocess
import tempfile
import json
from datetime import datetime, timezone
from library import fingerprint

ROOT = Path(__file__).resolve().parents[1]
parser = argparse.ArgumentParser()
parser.add_argument('--sanitize', action='store_true')
parser.add_argument('--filter', default='*.test.cpp', help='Optional test filename glob')
args = parser.parse_args()
report = ROOT / '.verify-results/latest.json'
report.parent.mkdir(exist_ok=True)
report.unlink(missing_ok=True)
before = fingerprint()
selected = sorted((ROOT / 'verify').glob(args.filter))
if not selected:
    raise RuntimeError('No tests matched the requested filter.')
flags = ['-std=c++17', '-O1', '-g', '-Wall', '-Wextra', '-Wno-sign-compare', '-iquote', str(ROOT)]
flags += ['-iquote', str(ROOT/'ac-library')]
if args.sanitize:
    flags += ['-fsanitize=address,undefined', '-fno-omit-frame-pointer', '-fno-pie', '-no-pie']
with tempfile.TemporaryDirectory(prefix='cp-library-') as temp:
    for source in selected:
        binary = str(Path(temp) / source.stem)
        print(f'Compile {source.name}', flush=True)
        subprocess.run([os.environ.get('CXX', 'g++'), *flags, str(source), '-o', binary], check=True)
        subprocess.run([binary], check=True, timeout=120)
if before != fingerprint():
    raise RuntimeError('Sources changed while testing; rerun verification.')
report.write_text(json.dumps({'fingerprint': before, 'passed': True,
    'sanitizers': args.sanitize, 'compiler': os.environ.get('CXX', 'g++'),
    'completed_at': datetime.now(timezone.utc).isoformat(),
    'tests': [p.relative_to(ROOT).as_posix() for p in selected]
}, indent=2))
print('All offline verification passed.', flush=True)
