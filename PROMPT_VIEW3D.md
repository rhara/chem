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
- 戻り値: `py3Dmol.view`インスタンス(内部で`.show()`は呼ばない。notebookのセルの最終式として使うか、呼び出し側で`.show()`を呼んで表示する)

## 実装方法

1. `path`のPDBテキストを読み込む
2. `HETATM`で始まる行から残基名(18-20列目、0-indexedで`line[17:20]`)を集め、`exclude`との差集合をソートして`ligand_resnames`とする(この部分は`_ligand_resnames(pdb_text, exclude)`としてテスト可能な形で切り出す)
3. `py3Dmol.view(width=width, height=height)`を作り、`addModel(pdb_text, "pdb")`
4. `setStyle({"cartoon": {"color": "spectrum", "colorscheme": "roygb"}})` -- "spectrum"単体だとN末端が紫がかった3Dmol.jsデフォルトのsinebowになるため、`colorscheme="roygb"`で青(N)→水色→緑→黄→橙→赤(C)の通常の配色にする
5. `ligand_resnames`が空でなければ`addStyle({"resn": ligand_resnames}, {"stick": {"colorscheme": "yellowCarbon"}})`(カートゥーンだけではリガンドが描画されないため)
6. `zoomTo()`して`view`を返す(`.show()`は呼ばない -- 戻り値をnotebookセルの最終式にすれば`_repr_html_`経由で自動表示される)

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- 新規の外部依存追加は不要(既存の`py3dmol`のみ使用)
- `src/chem/view3d/`として`chem`パッケージのサブパッケージに配置(`pyproject.toml`の`[tool.setuptools.packages.find]`が`chem.view3d`を自動検出)

## サンプルノートブックの更新

`notebooks/example_proteins.ipynb`の「View one of the downloaded structures」セルを、ベタ書きのpy3Dmol呼び出しから`view3d.render_protein`呼び出しに置き換える。`ipywidgets.Output()`のコンテキスト内で使うため、この呼び出しでは戻り値に対して明示的に`.show()`を呼ぶ(`view3d.render_protein(os.path.join("data", pdb_filename), exclude=_display_exclude).show()`)。`_display_exclude = SOLVENT_AND_IONS | {"NAG", "TYS", "MRD"}`(グリコシル化糖鎖・スルホチロシン・結晶化添加剤MRD)はnotebook側にそのまま残す(トロンビン構造セット固有の除外リストのため、関数のデフォルトには含めない)。

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- テストは`tests/test_view3d.py`にネットワーク・ブラウザ不要な範囲(`_ligand_resnames`の除外ロジック、`render_protein`が返す`py3Dmol.view`の内部JS文字列に期待するスタイルコマンドが含まれる/含まれないことの確認)のみ追加する。実際のブラウザ上での表示確認は手動で行う
