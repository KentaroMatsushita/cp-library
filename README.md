# cp_library

競技プログラミング用 C++ ライブラリ。日本語の使い方・前提条件・計算量と、再現可能な自動検証を整備しています。

- [説明サイト](https://kentaromatsushita.github.io/cp-library/)
- [公開リポジトリ](https://github.com/KentaroMatsushita/cp-library)
- [自動テスト](https://github.com/KentaroMatsushita/cp-library/actions/workflows/library.yml)

公開版には明示したライブラリと説明・テストだけを含めます。元の作業用リポジトリの Git 履歴、コンテスト回答、入力・出力、実行ファイルは移しません。

## 幾何

`geometry/geometry.hpp` は C++17 で単独 include できます。

```cpp
#include "geometry/geometry.hpp"
using namespace cp_geometry;

int main() {
    ConvexPolygon p({{0,0},{4,0},{4,3},{0,3}});
    auto where = p.contains({1,1}); // IN
    auto hit = p.line_intersections(Line({-1,1},{5,1})); // (0,1), (4,1)
    auto diameter_squared = convex_diameter2(p.vertices()); // 25
}
```

凸包と探索用インデックスの前処理は `O(N log N)`。内外判定、指定方向の最遠点、直線との境界交点は各 `O(log N)`。半平面交差には矩形内の `O(N log N)` 版と、空・非有界・多角形・退化を区別する一般版（期待 `O(N log N)`、最悪 `O(N²)`）があります。

座標は `long double`、絶対誤差 `1e-12` です。一般版の退化した領域では実行可能点だけを返し、線分の両端は返しません。詳細は説明サイトまたは `docs/geometry.html` を参照してください。

## 検証・サイト生成

必要環境は Python 3.10+ と GCC/Clang の C++17。追加の Python パッケージや外部ジャッジ接続は不要です。一部の既存 SIMD ライブラリのテストには AVX2/BMI2/POPCNT が必要です。

```sh
python3 tools/test.py --sanitize
python3 tools/build_docs.py
python3 tools/check_docs.py
python3 -m http.server 8000 --directory _site
```

テストは素朴解とのランダム比較、境界ケースの回帰、AddressSanitizer/UndefinedBehaviorSanitizer を使います。全 API の網羅検証ではなく、各ページに対象のテストを表示します。説明があることとテスト成功は別の状態です。未実装・未完成の既存ファイルも、その状態を明記して残しています。

`CXX=clang++ python3 tools/test.py --sanitize` でコンパイラを変更、`--filter geometry.test.cpp` でテストを絞れます。テスト失敗時には古い成功結果を削除し、ソース・テストのハッシュが一致する場合だけ成功表示を出します。

既存の `.cpp` の多くは共通テンプレートやマクロに依存します。同名の別実装をまとめて include せず、各ページの前提条件を確認してください。引用ヘッダのパス指定には `-iquote /path/to/cp-library` を使えます。

## GitHub Actions と Pages

push / pull request でテスト、ページ生成、内部リンクと公開対象の検査を実行します。テスト失敗時には公開しません。

Pages 公開は初回だけ Settings → Pages で Source を **GitHub Actions** にし、Settings → Secrets and variables → Actions → Variables に **PUBLISH_LIBRARY_PAGES=true** を設定します。以後は main の成功ビルドから自動公開します。非公開の元リポジトリではこの変数を設定せず、検証だけ走らせられます。

公開版を直接更新する場合は通常通り commit / push してください。元の作業フォルダから公開用スナップショットを作る場合は、検証後に次を実行します。

```sh
python3 tools/export_public.py --output /path/to/separate-public-export
```

出力は公開対象のコード・説明・検証に限定し、元の `.git` はコピーしません。公開対象ディレクトリは `tools/library.py`、説明は `docs/catalog.json` と `docs/notes*.json`、検証範囲は `docs/coverage.json` で管理します。

## 出典

既存コードに記載された元作者・出典・ライセンスのコメントは維持しています。AtCoder Library の最大流ヘッダと必要な内部ヘッダを検証用依存として同梱し、そのライセンスは `ac-library/atcoder/LICENSE` に保持します。リポジトリ全体へ新しい一律のライセンスを付ける変更はしていません。
