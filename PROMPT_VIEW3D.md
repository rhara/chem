# chem.view3d.render_protein の再現プロンプト

`chem`リポジトリの`chem`パッケージ内に、py3Dmolによる構造表示用のサブパッケージ`view3d`(`chem.view3d`)を追加し、単一のPDBファイルをインタラクティブに表示する関数を実装するための指示。`chem.protein`の`SOLVENT_AND_IONS`を再利用する。

## 背景

`notebooks/example_proteins.ipynb`の「View one of the downloaded structures」セルには、カートゥーン表示+リガンドのstick表示を行うpy3Dmol用のローカル関数`render`がベタ書きされていた。これを`chem.view3d.render_protein`として切り出し、notebook側は

```python
from chem import view3d
view3d.render_protein(pdb_filename)
```

の形で呼べるようにする。

## 要件

`src/chem/view3d/render.py` に関数`render_protein`を実装し、`src/chem/view3d/__init__.py`で`from .render import render_protein`として再エクスポートする。以下の形で呼び出せること:

```python
from chem import view3d
view3d.render_protein(path, exclude=SOLVENT_AND_IONS, width=600, height=500)
```

- `path`: PDBファイルへのパス
- `exclude`: リガンドのstick表示から除外するHETコードの集合。デフォルトは`chem.protein.SOLVENT_AND_IONS`。呼び出し側で`SOLVENT_AND_IONS | {"NAG", "TYS"}`のように追加のコードを合わせて渡せば、構造固有の非リガンドHETATM(糖鎖修飾、修飾残基など)も除外できる
- `width`/`height`: ビューアのピクセルサイズ
- 戻り値: なし(`None`)。ビューを`view.show()`で明示的に表示し、続けてビューの下にキャプションを表示する副作用のみを持つ。呼び出し側は`view3d.render_protein(path)`と呼ぶだけでよく、`.show()`を連鎖させたりnotebookセルの最終式として使う必要はない(そうすると二重表示になる)

## 実装方法

1. `path`のPDBテキストを読み込む
2. `HETATM`で始まる行から残基名(18-20列目、0-indexedで`line[17:20]`)を集め、`exclude`との差集合をソートして`ligand_resnames`とする(この部分は`_ligand_resnames(pdb_text, exclude)`としてテスト可能な形で切り出す)
3. キャプション用の付随情報を集める(いずれもテスト可能なヘルパー関数に切り出す):
   - `_chain_ids(pdb_text)`: `ATOM`で始まる行のchain列(22列目、0-indexedで`line[21]`)から重複を除きソートしたリスト。`ATOM`行が無ければ空リスト
   - `_resolution(pdb_text)`: legacy PDB形式のヘッダにある`REMARK   2 RESOLUTION.    N.NN ANGSTROMS.`行を正規表現`^REMARK\s+2\s+RESOLUTION\.\s+([\d.]+)\s+ANGSTROMS\.`(`re.MULTILINE`)でパースし、`"N.NN Å"`の形式で返す。マッチしなければ`"N/A"`(NMR構造、AlphaFold予測構造、`chem.protein.align`の出力(Bio.PDBの`PDBIO`はヘッダ/REMARKを保持しないため)はいずれもこのケースになる)
   - `_caption(path, pdb_text, ligand_resnames)`: `path`のファイル名(拡張子除く)をPDB IDとして、`f"PDB ID: {pdb_id} &nbsp;|&nbsp; Chain: {chains} &nbsp;|&nbsp; Ligand: {ligands} &nbsp;|&nbsp; Resolution: {resolution}"`の形式の文字列を組み立てる。`chains`は`_chain_ids`の結果をカンマ区切りにしたもの(空なら`"N/A"`)、`ligands`は`ligand_resnames`をカンマ区切りにしたもの(空なら`"none"`)
4. ビュー構築部分は`_build_view(pdb_text, ligand_resnames, width, height)`としてテスト可能な形で切り出す:
   - `py3Dmol.view(width=width, height=height)`を作り、`addModel(pdb_text, "pdb")`
   - `setStyle({"cartoon": {"color": "spectrum", "colorscheme": "roygb"}})` -- "spectrum"単体だとN末端が紫がかった3Dmol.jsデフォルトのsinebowになるため、`colorscheme="roygb"`で青(N)→水色→緑→黄→橙→赤(C)の通常の配色にする
   - `ligand_resnames`が空でなければ`addStyle({"resn": ligand_resnames}, {"stick": {"color": "magenta"}})`(カートゥーンだけではリガンドが描画されないため)。単色`"color": "magenta"`を使う -- カートゥーンの`roygb`スペクトルに対してインパクトのある濃い色でコントラストを出すため。`"colorscheme": "yellowCarbon"`は不採用(スペクトルの黄色部分と衝突する)、また`"pinkCarbon"`は3Dmol.jsに存在しない
   - `zoomTo()`して`view`を返す
5. `render_protein`本体: `_build_view(...)`で得た`view`に対して`view.show()`を呼びビューを表示し、その直後に`IPython.display.display(HTML(f"<b>{caption}</b>"))`でキャプションをビューの下に表示する(`display`/`HTML`は`IPython.display`からimport)。`view`は返さない

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- 新規の外部依存追加は不要(既存の`py3dmol`のみ使用)
- `src/chem/view3d/`として`chem`パッケージのサブパッケージに配置(`pyproject.toml`の`[tool.setuptools.packages.find]`が`chem.view3d`を自動検出)

## サンプルノートブックの更新

`notebooks/example_proteins.ipynb`の「View one of the downloaded structures」セルを、ベタ書きのpy3Dmol呼び出しから`view3d.render_protein`呼び出しに置き換える。`ipywidgets.Output()`のコンテキスト内で使うが、`render_protein`は戻り値を持たず自身で表示まで完結するため、`.show()`は連鎖させず`view3d.render_protein(os.path.join("data", pdb_filename), exclude=_display_exclude)`とだけ呼ぶ。`_display_exclude = SOLVENT_AND_IONS | {"NAG", "TYS", "MRD"}`(グリコシル化糖鎖・スルホチロシン・結晶化添加剤MRD)はnotebook側にそのまま残す(トロンビン構造セット固有の除外リストのため、関数のデフォルトには含めない)。

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- テストは`tests/test_view3d.py`にネットワーク・ブラウザ不要な範囲(`_ligand_resnames`の除外ロジック、`_chain_ids`のchain収集ロジック、`_resolution`のREMARK 2パースロジック(値あり/なし双方)、`_caption`の文字列組み立て、`_build_view`が返す`py3Dmol.view`の内部JS文字列に期待するスタイルコマンドが含まれる/含まれないことの確認、`py3Dmol.view.show`と`chem.view3d.render.display`の両方を`monkeypatch`で差し替えて「ビュー表示→キャプション表示」の順序と`render_protein`が`None`を返すことを確認)のみ追加する。実際のブラウザ上での表示確認は手動で行う
