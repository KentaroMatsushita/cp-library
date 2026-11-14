#!/usr/bin/env python3
"""Dependency-free static documentation builder. Only allowlisted sources are copied."""
import argparse
import html
import json
from pathlib import Path
import re
import shutil
import hashlib
from urllib.parse import quote
from library import ROOT, catalog, sources, fingerprint

parser=argparse.ArgumentParser()
parser.add_argument('--output',default='_site')
args=parser.parse_args()
out=Path(args.output).resolve()
out.mkdir(parents=True,exist_ok=True)
files=sources()
metadata=catalog()
assert set(metadata)<=set(files), 'Catalog refers to a source outside the publication allowlist'
for src, item in metadata.items():
    for key in ('title','summary','complexity','contract','usage','tests'):
        assert key in item, (src,key)
    for test in item['tests']:
        assert (ROOT/test).is_file(), test
report_path=ROOT/'.verify-results/latest.json'
report=json.loads(report_path.read_text()) if report_path.exists() else {}
current=report.get('passed') and report.get('fingerprint')==fingerprint()
e=html.escape
def slug(path):
    return re.sub(r'[^A-Za-z0-9_-]', '_', path)+'-'+hashlib.sha256(path.encode()).hexdigest()[:8]+'.html'
def status(path):
    if path not in metadata: return '未レビュー', 'pending'
    tests=metadata[path]['tests']
    if tests and current and all(t in report.get('tests',[]) for t in tests): return 'テスト対象・成功', 'passed'
    if tests: return 'テスト結果なし／更新後未実行', 'pending'
    return '説明あり・未テスト', 'pending'

def layout(title,body,prefix=''):
    return f'''<!doctype html><html lang="ja"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1"><title>{e(title)} · cp_library</title>
<link rel="stylesheet" href="{prefix}assets/style.css"></head><body>
<header><a class="brand" href="{prefix}index.html">cp<span>_</span>library</a><nav><a href="{prefix}index.html">ライブラリ一覧</a><a href="{prefix}guide.html">使い方・検証</a><a href="https://github.com/KentaroMatsushita/cp-library">GitHub ↗</a></nav></header>
<main>{body}</main><footer>cp_library</footer>
<script src="{prefix}assets/main.js" defer></script></body></html>'''

groups={}
for path in files: groups.setdefault(path.split('/')[0] if '/' in path else 'template',[]).append(path)
verified=sum(status(p)[1]=='passed' for p in files)
body=f'''<h1>cp_library</h1>
<p>競技プログラミング用の C++ ライブラリ。</p>
<section class="catalog"><aside><h2>カテゴリ</h2>{''.join(f'<a href="#cat-{e(k)}">{e(k)} <span>{len(v)}</span></a>' for k,v in groups.items())}</aside><div>
<label class="search">ライブラリを探す<input id="search" type="search" placeholder="名前・用途・ファイル名で検索" autocomplete="off"></label>
<label class="filter"><input id="documented" type="checkbox"> 現在のテストで確認したものだけ表示</label>
<p class="verification-note">テスト状況・前提条件・使い方は各ファイルのページに記載しています。</p><p id="empty" hidden>該当するライブラリがありません。</p>'''
for group, paths in groups.items():
    body+=f'<section class="group" id="cat-{e(group)}"><h2>{e(group)}</h2><ul>'
    for path in paths:
        item=metadata.get(path,{})
        title=item.get('title',Path(path).stem)
        text=item.get('summary','既存実装の目録。説明・コードレビュー・テストは未整備です。')
        label,cls=status(path)
        body+=f'''<li class="entry" data-documented="{int(cls=='passed')}" data-search="{e((path+' '+title+' '+text).lower())}"><a href="library/{slug(path)}"><span class="entry-name">{e(title)}</span> <code>({e(path)})</code> <span class="badge {cls}">{e(label)}</span></a></li>'''
    body+='</ul></section>'
body+='</div></section>'
(out/'index.html').write_text(layout('ライブラリ一覧',body),encoding='utf-8')
(out/'library').mkdir(exist_ok=True)
for path in files:
    source=(ROOT/path).read_text(errors='replace')
    item=metadata.get(path,{})
    title=item.get('title',Path(path).stem)
    label,cls=status(path)
    body=f'<p class="breadcrumb"><a href="../index.html">ライブラリ</a> / {e(path)}</p><article><h1>{e(title)}</h1><span class="badge {cls}">{e(label)}</span>'
    if item:
        body+=f'<p class="lead">{e(item["summary"])}</p><h2>計算量</h2><p>{e(item["complexity"])}</p><h2>前提条件・戻り値</h2><p>{e(item["contract"])}</p><h2>使い方</h2><pre><code>{e(item["usage"])}</code></pre>'
        if 'detail' in item:
            detail=(ROOT/'docs'/item['detail']).resolve()
            assert detail.is_relative_to((ROOT/'docs').resolve())
            body+=detail.read_text()
        body+='<h2>検証対象</h2>'
        if item['tests']:
            body+='<ul>'+''.join(f'<li><a href="../tests/{quote(Path(t).name)}.html">{e(t)}</a></li>' for t in item['tests'])+'</ul>'
        else: body+='<p>動作テストはまだ登録されていません。説明の整備と動作検証は区別しています。</p>'
    else:
        body+='<p class="notice">既存実装の目録です。コードを読んだ説明の整備・レビュー・テストは未完了です。ファイル名だけから API の意味や計算量を推定していません。</p>'
        comments=[s.strip().lstrip('/ *') for s in source.splitlines()[:50] if s.strip().startswith(('//',' *'))]
        if comments: body+='<h2>ソースに記載されたコメント</h2><pre>'+e('\n'.join(comments))+'</pre>'
    includes=re.findall(r'^\s*#include\s*"([^"]+)"',source,re.M)
    if includes:
        body+='<h2>明示された依存</h2><ul>'
        for dep in includes:
            relative=(Path(path).parent/dep).as_posix()
            resolved=relative if relative in files else dep
            body+=f'<li><a href="{slug(resolved)}">{e(dep)}</a></li>' if resolved in files else f'<li><code>{e(dep)}</code>（外部または公開対象外）</li>'
        body+='</ul>'
    body+=f'<h2>ソースコード</h2><p><a href="../source/{quote(path)}" download>ファイルをダウンロード</a> · <a href="https://github.com/KentaroMatsushita/cp-library/blob/main/{quote(path)}">GitHub で見る ↗</a></p><button class="copy" type="button" data-copy="source">コードをコピー</button><pre><code id="source">{e(source)}</code></pre></article>'
    (out/'library'/slug(path)).write_text(layout(title,body,'../'),encoding='utf-8')
    target=out/'source'/path
    target.parent.mkdir(parents=True,exist_ok=True)
    target.write_text(source,encoding='utf-8')
(out/'tests').mkdir(exist_ok=True)
for test in sorted((ROOT/'verify').glob('*.test.cpp')):
    body=f'<article><h1>{e(test.name)}</h1><p>公開テストの実装。固定シードで再現可能な検証を行います。</p><pre><code>{e(test.read_text())}</code></pre></article>'
    (out/'tests'/(test.name+'.html')).write_text(layout(test.name,body,'../'),encoding='utf-8')
guide=(ROOT/'docs/guide.html').read_text()
report_text=f'現在のソースに一致するテスト成功結果があります。実行時刻: {report.get("completed_at","")}' if current else '現在のソースに一致するテスト結果がありません。テストを実行すると表示が更新されます。'
(out/'guide.html').write_text(layout('使い方・検証',f'<article><h1>使い方・検証</h1><p class="notice">{e(report_text)}</p>{guide}</article>'),encoding='utf-8')
shutil.copytree(ROOT/'docs/assets',out/'assets',dirs_exist_ok=True)
(out/'.nojekyll').touch()
(out/'manifest.json').write_text(json.dumps({'sources':files,'documented':sorted(metadata),'verified':verified,'fingerprint':fingerprint()},ensure_ascii=False,indent=2),encoding='utf-8')
print(f'Built {len(files)} source pages; {len(metadata)} documented; {verified} covered by current passing tests -> {out}')
