#!/usr/bin/env python3
"""Export the current curated sources without copying the private Git history."""
import argparse
import json
from pathlib import Path
import shutil
from library import ROOT, sources, SUPPORT_FILES

p=argparse.ArgumentParser()
p.add_argument('--output',required=True)
args=p.parse_args()
out=Path(args.output).resolve()
if out==ROOT or ROOT.is_relative_to(out):
    raise ValueError('Choose a separate export directory, not the original repository or a parent.')
marker=out/'.public-export.json'
if out.exists() and any(out.iterdir()) and not marker.exists():
    raise ValueError('Nonempty output is not a previous curated export.')
paths=sources()+list(SUPPORT_FILES)+['README.md','.gitignore','.github/workflows/library.yml']
for directory in ('docs','tools','verify'):
    paths += [f.relative_to(ROOT).as_posix() for f in (ROOT/directory).rglob('*')
              if f.is_file() and not f.is_symlink() and '__pycache__' not in f.parts]
paths=sorted(set(paths))
out.mkdir(parents=True,exist_ok=True)
old=json.loads(marker.read_text())['files'] if marker.exists() else []
for stale in set(old)-set(paths):
    target=(out/stale).resolve()
    if not target.is_relative_to(out) or '.git' in target.relative_to(out).parts:
        raise ValueError('Invalid stale path in export manifest.')
    target.unlink(missing_ok=True)
for relative in paths:
    origin=ROOT/relative
    target=(out/relative).resolve()
    if not target.is_relative_to(out): raise ValueError(relative)
    target.parent.mkdir(parents=True,exist_ok=True)
    shutil.copyfile(origin,target)
marker.write_text(json.dumps({'files':paths},ensure_ascii=False,indent=2),encoding='utf-8')
print(f'Exported {len(paths)} files, without private history, to {out}')
