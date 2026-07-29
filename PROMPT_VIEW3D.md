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
view3d.render_protein(
    path, exclude=SOLVENT_AND_IONS, width=600, height=500,
    coloring="spectrum", bfactor_range=(50, 90),
)
```

- `path`: PDBファイルへのパス
- `exclude`: リガンドのstick表示から除外するHETコードの集合。デフォルトは`chem.protein.SOLVENT_AND_IONS`。呼び出し側で`SOLVENT_AND_IONS | {"NAG", "TYS"}`のように追加のコードを合わせて渡せば、構造固有の非リガンドHETATM(糖鎖修飾、修飾残基など)も除外できる
- `width`/`height`: ビューアのピクセルサイズ
- `coloring`: カートゥーンの配色方式。`"spectrum"`(デフォルト) -- 残基位置でN末端→C末端をレインボー表示。`"bfactor"` -- ファイルのB-factor列(例: AlphaFoldが格納するper-residue pLDDT信頼度)でレインボー表示
- `bfactor_range`: `coloring="bfactor"`のときのグラデーションの`(min, max)`(`"spectrum"`では無視)。デフォルトはAlphaFoldのpLDDT信頼度の慣習に合わせた`(50, 90)`。結晶構造の温度因子(B-factor)を使う場合は構造自体のB-factor範囲を渡す
- `COLORINGS = ("spectrum", "bfactor")`をモジュールレベル定数として定義し、`coloring`がこれ以外の値なら`ValueError`
- 戻り値: なし(`None`)。ビュー+キャプションを表示する副作用のみを持つ。呼び出し側は`view3d.render_protein(path)`と呼ぶだけでよく、`.show()`を連鎖させたりnotebookセルの最終式として使う必要はない(そうすると二重表示になる)
- レイアウト: ビューは薄いグレー(`#ccc`)の枠で囲む -- 3Dmol.jsのマウス操作(回転・ズーム・パン)が効く範囲そのものを枠として可視化し、ユーザーが操作領域を把握できるようにする。キャプションはビューの下ではなく右側に、1プロパティ1行で表示する

## 実装方法

1. `path`のPDBテキストを読み込む
2. `HETATM`で始まる行から残基名(18-20列目、0-indexedで`line[17:20]`)を集め、`exclude`との差集合をソートして`ligand_resnames`とする(この部分は`_ligand_resnames(pdb_text, exclude)`としてテスト可能な形で切り出す)
3. キャプション用の付随情報を集める(いずれもテスト可能なヘルパー関数に切り出す):
   - `_chain_ids(pdb_text)`: `ATOM`で始まる行のchain列(22列目、0-indexedで`line[21]`)から重複を除きソートしたリスト。`ATOM`行が無ければ空リスト
   - `_resolution(pdb_text)`: legacy PDB形式のヘッダにある`REMARK   2 RESOLUTION.    N.NN ANGSTROMS.`行を正規表現`^REMARK\s+2\s+RESOLUTION\.\s+([\d.]+)\s+ANGSTROMS\.`(`re.MULTILINE`)でパースし、`"N.NN Å"`の形式で返す。マッチしなければ`"N/A"`(NMR構造、AlphaFold予測構造、`chem.protein.align`の出力(Bio.PDBの`PDBIO`はヘッダ/REMARKを保持しないため)はいずれもこのケースになる)
   - `_caption_lines(path, pdb_text, ligand_resnames)`: `path`のファイル名(拡張子除く)をPDB IDとして、`["PDB ID: {pdb_id}", "Chain: {chains}", "Ligand: {ligands}", "Resolution: {resolution}"]`の4行からなるリストを返す(1プロパティ1行、単一の文字列に結合しない)。`chains`は`_chain_ids`の結果をカンマ区切りにしたもの(空なら`"N/A"`)、`ligands`は`ligand_resnames`をカンマ区切りにしたもの(空なら`"none"`)
4. ビュー構築部分は`_build_view(pdb_text, ligand_resnames, width, height, coloring, bfactor_range)`としてテスト可能な形で切り出す:
   - `py3Dmol.view(width=width, height=height)`を作り、`addModel(pdb_text, "pdb")`
   - `coloring == "bfactor"`なら`setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": bfactor_range[0], "max": bfactor_range[1]}}})`
   - それ以外(`"spectrum"`)なら`setStyle({"cartoon": {"color": "spectrum", "colorscheme": "roygb"}})` -- "spectrum"単体だとN末端が紫がかった3Dmol.jsデフォルトのsinebowになるため、`colorscheme="roygb"`で青(N)→水色→緑→黄→橙→赤(C)の通常の配色にする
   - `ligand_resnames`が空でなければ`addStyle({"resn": ligand_resnames}, {"stick": {"color": "magenta"}})`(カートゥーンだけではリガンドが描画されないため)。単色`"color": "magenta"`を使う -- カートゥーンの`roygb`スペクトルに対してインパクトのある濃い色でコントラストを出すため。`"colorscheme": "yellowCarbon"`は不採用(スペクトルの黄色部分と衝突する)、また`"pinkCarbon"`は3Dmol.jsに存在しない
   - `zoomTo()`して`view`を返す
5. `render_protein`本体: 冒頭で`coloring`を`COLORINGS`と照合しなければ`ValueError`。`_build_view(...)`で`view`を作った後、表示は次の2段階で行う:
   - `uuid.uuid4().hex`から一意な`frame_id`(例: `f"chem-view3d-{uuid.uuid4().hex}"`)を作り、`display(HTML(...))`で「幅・高さを`width`/`height`に合わせ、`border:1px solid #ccc;`を付けた空の`<div id="{frame_id}">`」と「`_caption_lines`の各行を`<div><b>{line}</b></div>`にした`<div>`」を`display:flex; align-items:flex-start; gap:16px;`の親`<div>`で横並びにしたHTMLを表示する。枠を先に幅・高さ・ボーダー込みで表示しておくことで、3Dmol.jsが非同期でCDNから読み込まれる間もレイアウトが安定する
   - 続けて`view.insert(frame_id)`を呼ぶ -- py3Dmolの公開API`insert()`は内部で自分のビュー用divを生成した上で、`document.getElementById(frame_id).append(...)`するスクリプトも含めて`publish_display_data`する。これにより直前の`display()`で作った空枠の中に実際の3Dmolビューが挿入される(`view.show()`や`view._make_html()`のような非公開実装には依存しない)
   - `view`は返さない

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- 新規の外部依存追加は不要(既存の`py3dmol`のみ使用)
- `src/chem/view3d/`として`chem`パッケージのサブパッケージに配置(`pyproject.toml`の`[tool.setuptools.packages.find]`が`chem.view3d`を自動検出)

## サンプルノートブックの更新

`notebooks/example_proteins.ipynb`の2つのセルを`view3d.render_protein`呼び出しに置き換える:

- 「View one of the downloaded structures」: ベタ書きのpy3Dmol呼び出しを`view3d.render_protein`呼び出しに置き換える。`ipywidgets.Output()`のコンテキスト内で使うが、`render_protein`は戻り値を持たず自身で表示まで完結するため、`.show()`は連鎖させず`view3d.render_protein(os.path.join("data", pdb_filename), exclude=_display_exclude)`とだけ呼ぶ(`coloring`は指定せず既定の`"spectrum"`のまま)。`_display_exclude = SOLVENT_AND_IONS | {"NAG", "TYS", "MRD"}`(グリコシル化糖鎖・スルホチロシン・結晶化添加剤MRD)はnotebook側にそのまま残す(トロンビン構造セット固有の除外リストのため、関数のデフォルトには含めない)
- 「View the predicted structure, colored by pLDDT confidence」: ベタ書きのpy3Dmol呼び出し(`setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": 50, "max": 90}}})`)を`view3d.render_protein(os.path.join("af_data", pdb_files[0]), coloring="bfactor", width=600, height=600)`に置き換える(`bfactor_range`は既定の`(50, 90)`のままでpLDDTの慣習と一致するため省略可)

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- テストは`tests/test_view3d.py`にネットワーク・ブラウザ不要な範囲(`_ligand_resnames`の除外ロジック、`_chain_ids`のchain収集ロジック、`_resolution`のREMARK 2パースロジック(値あり/なし双方)、`_caption_lines`が返すリストの内容、`_build_view`が返す`py3Dmol.view`の内部JS文字列に`coloring="spectrum"`/`"bfactor"`(指定した`bfactor_range`込み)それぞれで期待するスタイルコマンドが含まれる/含まれないことの確認、`render_protein`が不正な`coloring`で`ValueError`になることの確認、`chem.view3d.render.display`と`py3Dmol.view.insert`の両方を`monkeypatch`で差し替えて「枠+キャプション表示→ビューのinsert」の順序・`display`に渡されたHTMLに含まれる`frame_id`が`insert`に渡された`containerid`と一致すること・`render_protein`が`None`を返すことを確認)のみ追加する。実際のブラウザ上での表示確認は手動で行う
