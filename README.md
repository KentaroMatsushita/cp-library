# cp_library

競技プログラミング用 C++ ライブラリ。日本語の使い方・前提条件・計算量と、再現可能な自動検証を整備しています。

- [説明サイト](https://kentaromatsushita.github.io/cp-library/)
- [公開リポジトリ](https://github.com/KentaroMatsushita/cp-library)
- [自動テスト](https://github.com/KentaroMatsushita/cp-library/actions/workflows/library.yml)

公開版には明示したライブラリと説明・テストだけを含めます。元の作業用リポジトリの Git 履歴、コンテスト回答、入力・出力、実行ファイルは移しません。

## 幾何

C++17、`cp_geometry` 名前空間。用途に応じて `geometry/base.hpp`（基本）、`convex.hpp`（凸多角形）、`halfplane.hpp`（半平面）、`circle.hpp`（円）、`algorithms.hpp`（最近点対等）を単独 include できます。全体は `geometry/geometry.hpp`。共通テンプレートには依存しません。

```cpp
#include "geometry/geometry.hpp"
namespace cg=cp_geometry;

int main() {
    cg::ConvexPolygon p({{0,0},{4,0},{4,3},{0,3}});
    auto where=p.contains({1,1}); // cg::Location::inside
    auto hit=p.line_intersections(cg::Line({-1,1},{5,1}));
    std::vector<cg::HalfPlane> h{{-1,0,0},{0,-1,0},{1,0,4},{0,1,3}};
    auto region=cg::half_plane_intersection(h); // a*x+b*y<=c
    cg::Circle a({0,0},1), b({1,0},1);
    auto lens=cg::circle_overlap_area(a,b);
    auto clipped=cg::circle_polygon_intersection_area(a,p.vertices());
}
```

内外判定・直線交点は各 O(log N)。既存の厳密凸包は `from_convex_hull` で O(N) 構築、任意点集合からは O(N log N)。円と多角形の共通面積 O(N)、円と円の共通面積 O(1)、最近点対 O(N log N)、厳密凸包の Minkowski 和 O(N+M) です。

半平面交差は `std::vector<cg::HalfPlane>{{a,b,c}, ...}` を受け取ります。矩形内版 O(N log N)、一般版は期待 O(N log N)・最悪 O(N²) で、空・非有界・有界多角形・有界退化を区別します。退化領域では実行可能点だけを返します。

判定定数の旧 `IN/ON/OUT` は入力関数と競合するため廃止し、`Location::inside/boundary/outside` に変更しました。普段の `my_template.hpp` と両方の include 順でテストしています。座標は long double、絶対誤差 1e-12。精度・退化・前提条件は各説明ページを参照してください。

## Voronoi / Delaunay

```cpp
#include "geometry/voronoi.hpp"
#include "geometry/power_diagram.hpp"
namespace cg=cp_geometry;

std::vector<cg::Point> sites{{-1,0},{1,0},{0,2}};
cg::VoronoiDiagram diagram(sites);
auto cells=diagram.cells({-5,-5},{5,5}); // 通常の全セル O(N log N)
auto largest=cg::largest_empty_circle_in_box(diagram,{-5,-5},{5,5});
for(const auto& edge:diagram.edges()) {
    auto clipped=edge.clip({-5,-5},{5,5}); // 線分・半直線・直線に対応
}
std::vector<cg::PowerSite> weighted{{{-1,0},4},{{1,0},0}};
auto cell=cg::power_cell(weighted,0,{-5,-5},{5,5}); // 重み付き1セル O(N log N)
```

通常の Delaunay/Voronoi は O(N log N) 時間・O(N) メモリ。重複点の代表添字を保持し、共線・共円にも対応します。Power Diagram の距離は `|p-site|²-weight` で、円周までの距離による図ではありません。Power の全セルは O(N² log N) です。

最大空円の領域制約は **中心の位置だけ** に適用します。円板全体が領域内に収まる最大内接円とは異なります。矩形版 O(N log N)、M頂点の単純多角形版 O(NM + N log N)。関連問題は Galaxy War（5th UCup Wuhan H）、Panda Preserve（WF2018 G）、Euclidean MST（Library Checker）などです。通常の頂点・交点座標は long double の近似値です。詳細は各ページを参照してください。

## 検証・サイト生成

必要環境は Python 3.10+ と GCC/Clang の C++17。追加の Python パッケージや外部ジャッジ接続は不要です。一部の既存 SIMD ライブラリのテストには AVX2/BMI2/POPCNT が必要です。

```sh
python3 tools/check_headers.py
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
