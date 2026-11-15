#!/usr/bin/env python3
"""Check local documentation links and the publication boundary."""
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import urlsplit, unquote
import json
import argparse
from library import sources, SUPPORT_FILES

parser=argparse.ArgumentParser()
parser.add_argument('--site',default='_site')
args=parser.parse_args()
root=Path(args.site).resolve()
class Links(HTMLParser):
    def __init__(self): super().__init__(); self.links=[]
    def handle_starttag(self,tag,attrs):
        for k,v in attrs:
            if k in ('href','src') and v: self.links.append(v)

failures=[]
for page in root.rglob('*.html'):
    p=Links();p.feed(page.read_text())
    for link in p.links:
        url=urlsplit(link)
        if url.scheme or url.netloc or not url.path: continue
        target=(page.parent/unquote(url.path)).resolve()
        if not target.is_relative_to(root) or not target.exists(): failures.append((str(page.relative_to(root)),link))
published={p.relative_to(root/'source').as_posix() for p in (root/'source').rglob('*') if p.is_file()}
assert published==set(sources())|set(SUPPORT_FILES), 'Published sources differ from the explicit allowlist'
assert not failures, failures
manifest=json.loads((root/'manifest.json').read_text())
assert set(manifest['sources'])|set(manifest['support_files'])==published
print(f'Documentation links and allowlist passed: {len(manifest["sources"])} sources and {len(SUPPORT_FILES)} support files.')
