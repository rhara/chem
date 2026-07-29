# chem.protein.align / chem.protein.find_pocket の再現プロンプト

`chem`リポジトリの`chem`パッケージ内に、構造解析用のサブパッケージ`protein`(`chem.protein`)を追加し、複数構造のシーケンス・3Dアラインメント(`align`)と、リガンド近傍のfpocketポケット・構成残基特定(`find_pocket`)を実装するための指示。`chem.rcsb`・`chem.alphafold`でダウンロードした構造をそのまま入力にできることを想定する。

## パッケージ構成

- `src/chem/protein/__init__.py`: `from .structural_align import align` / `from .pocket import find_pocket` として再エクスポート
- `src/chem/protein/structural_align.py`: `align`の実装
- `src/chem/protein/pocket.py`: `find_pocket`の実装

**注意**: 実装ファイル名を公開関数名と同じにしない(例: `align.py`に`align`関数を置かない)。パッケージの`__init__.py`で`from .align import align`のように再エクスポートすると、`chem.protein`パッケージの`align`属性がサブモジュールからその関数へ上書きされてしまい、以降`import chem.protein.align`のようなドット区切りでのサブモジュールアクセスが(内部テストなどで)壊れる。既存の`chembl`/`rcsb`/`alphafold`サブパッケージが実装ファイルを常に`fetch.py`(公開関数名と別名)にしているのも同じ理由。

## chem.protein.align

```python
from chem import protein
protein.align(structures, reference=None, chain=None, outdir="aligned")
```

- `structures`: 同一ターゲットの構造ファイルパスのリスト(PDB/CIF混在可。`chem.rcsb`/`chem.alphafold`でダウンロードしたファイルをそのまま渡せる)
- `reference`: アラインメントの基準構造。`structures`へのインデックス(int)またはパス文字列。省略時は`structures[0]`。**`structures`配列に含まれている必要はない** — 含まれていてもいなくても、`outdir`には常に1回だけ書き出される(含まれる場合はアラインメントループ側でスキップする)
- `chain`: 各構造で使うチェーンIDを明示指定(省略時は下記の自動選択)
- `outdir`: 座標変換後の構造を書き出す先のディレクトリ(存在しなければ作成)
- 戻り値: `{path: rmsd}`(シーケンスマッチしたCA原子上のRMSD。reference自身は`0.0`)

### アルゴリズム

1. `reference`を解決し、Bio.PDB(拡張子`.cif`/`.mmcif`なら`MMCIFParser`、それ以外は`PDBParser`)で読み込む
2. 各構造について「主鎖(primary polymer chain)」を選択する:
   - `chain`引数が指定されていればそのチェーンIDを使う(見つからなければ`ValueError`)
   - 指定がなければ、標準アミノ酸残基を最も多く含むチェーンを自動選択する。**単純に「ファイル内最初のチェーン」を採用してはいけない** — 例えばトロンビンはL鎖(軽鎖、~36残基)がH鎖(重鎖、~259残基)より先に現れることがあり、「最初に見つかった一定長以上のチェーン」という閾値ベースの選択だと軽鎖を誤って選んでしまう。残基数最大のチェーンを選ぶことで正しくH鎖(触媒ドメイン)が選ばれることを確認済み
   - 「標準アミノ酸残基」の判定は`Bio.PDB.Polypeptide.is_aa(residue, standard=True)`**だけでは不十分**: `is_aa`はresnameのみを見ており、hetero flag(`residue.id[0]`)をチェックしない。共有結合したペプチド模倣リガンド(D-アミノ酸を含むペプチド性阻害剤など)が蛋白質と同じチェーンIDを使い、かつ標準アミノ酸名(例: `PRO`)をHETATMとして持つ場合、`is_aa`だけの判定だとリガンド側の残基を誤って蛋白質配列に取り込んでしまう。実データ(トロンビン+ペプチド性阻害剤複合体 `6YHG`、リガンド側の`PRO H 307`がHETATMとして紐づいていた)でこの誤りによりRMSDが約0.3Å→約1.9Åに悪化する不具合を発見・修正済み。`_is_polymer_residue(r) = r.id[0] == " " and is_aa(r, standard=True)`という追加のhetero flagチェックを`_select_chain`・`_chain_seq_and_ca`の両方に適用する
3. 選択したチェーンから、CA原子を持つ標準アミノ酸残基(上記の`_is_polymer_residue`基準)のみを使って1文字シーケンスと対応するCA原子リストを作る
4. `Bio.Align.PairwiseAligner`(`mode="global"`、`open_gap_score=-10`、`extend_gap_score=-0.5`、`match_score=2`、`mismatch_score=-1`)でreferenceシーケンスと各構造のシーケンスをグローバルアラインメントし、`alignment.aligned`からギャップのないマッチ位置ペアを取り出し、対応するCA原子ペアのリストを作る
5. マッチしたCA原子ペアが3組未満なら(`_MIN_MATCHED_RESIDUES = 3`)、その構造はアラインメント不可としてスキップ(quiet時以外は`stderr`に警告を出し、処理は継続)
6. `Bio.PDB.Superimposer`でマッチしたCA原子ペアに対してKabsch法によるフィッティングを行い、得られた回転・並進を**構造の全原子**(タンパク質だけでなくリガンド・水も含む)に適用する。これにより、align出力後もリガンドの相対位置が保たれ、後続の`find_pocket`にそのまま使える
7. `Bio.PDB.PDBIO`で、referenceも含めて各構造を`outdir`に`{入力ファイルのstem}.pdb`として書き出す(**入力がCIFでも常にPDB形式で出力する** — fpocketがPDB形式を要求するための一貫性)。referenceの書き出しはメインループの前に1回だけ行い、`structures`側のループでは`path == ref_path`のときスキップする(referenceが`structures`に含まれない場合でも正しく書き出される)

## chem.protein.find_pocket

```python
from chem import protein
protein.find_pocket(structure, ligand=None, outdir=None)
```

- `structure`: PDBファイルのパス(fpocketはレガシーPDB形式を要求する。`chem.protein.align`の出力がそのまま使える)
- `ligand`: リガンドの指定方法。以下3種類すべてをサポートする:
  - `None`(デフォルト): `structure`内のHETATM残基から、水・イオン・一般的な結晶化添加剤(下記リスト)を除いた中で最も原子数の多い残基を自動検出する
  - 1〜3文字のPDB HETコード文字列(例: `"STI"`): `structure`内でそのresnameに一致するHETATM残基を使う(複数該当時は最初のものを使う)
  - 既存ファイルへのパス(`.pdb`/`.sdf`/`.mol`/`.mol2`): RDKit(`Chem.MolFromPDBFile`/`MolFromMolFile`/`MolFromMol2File`、拡張子で分岐)でリガンドの3D座標を読み込む。蛋白質ファイル自体にリガンドが含まれない場合(ドッキングポーズなど)向け
  - 判定順序: `os.path.isfile(ligand)`が真なら外部ファイル、そうでなく`^[A-Za-z0-9]{1,3}$`にマッチすればHETコード、どちらでもなければ`ValueError`
- `outdir`: fpocketの生ログ出力(`pockets/`、`*_info.txt`等)を保持するディレクトリ。`None`なら一時ディレクトリを使い、処理後に破棄する
- 戻り値: 以下を含む辞書
  - `pocket_id`: fpocketのポケット番号
  - `score` / `druggability_score` / `volume`: `_info.txt`から取り出した簡易フィールド
  - `residues`: `[{"chain":..., "resnum":..., "resname":...}, ...]`(該当ポケットを構成する残基、出現順、重複なし)
  - `info`: 該当ポケットのfpocketスコア情報全体(`_info.txt`の全フィールドを`{項目名: 値}`の辞書にしたもの)

### アルゴリズム

1. `structure`を作業ディレクトリ(`outdir`があればそこ、なければ一時ディレクトリ)にコピーし、`fpocket -f {basename}`をサブプロセスとして実行する(`cwd`を作業ディレクトリにすることで、元ファイルの場所を汚さない)。`fpocket`実行ファイルが見つからない場合はconda-forgeでのインストール方法を含む分かりやすい`RuntimeError`にする
2. `{stem}_out/pockets/pocket{N}_atm.pdb`(各ポケットを構成する蛋白質原子)と`{stem}_out/{stem}_info.txt`(ポケットごとのスコア情報)がfpocketの出力として得られる
3. `ligand`引数からリガンド原子の3D座標(numpy配列)を解決する(上記参照)
4. 各`pocket{N}_atm.pdb`をBio.PDBでパースし、その原子群とリガンド原子群との最小距離(全原子対distanceのmin)を計算する。最小距離が最も小さいポケット番号を採用する(3PTB(トリプシン)+BENリガンドで検証済み: 最もfpocketスコアが高いPocket 1が同時に最もBENに近く、既知の活性部位Asp189・Ser195を含むことを確認)
5. 採用したポケットの`pocket{N}_atm.pdb`から一意な`(chain, resnum, resname)`のリストを抽出する
6. `{stem}_info.txt`をパースする。フォーマットは`Pocket N :`という見出し行の後に、タブ区切りの`項目名 : \t値`行が続き、空行でポケットが区切られる。数値に変換できるものは`float`にする

### `SOLVENT_AND_IONS`(自動検出で除外するHETコード、非網羅的)

```
HOH, WAT, DOD,
NA, CL, K, MG, CA, ZN, MN, FE, FE2, CO, NI, CU, CD, LI, RB, CS, BR, IOD, NH4,
SO4, PO4, GOL, EDO, PEG, PG4, 1PE, P6G, MPD, FMT, ACT, DMS, TRS, BME, EPE, HEPES, IPA, UNX, UNL
```

`chem/protein/pocket.py`のモジュールレベル定数で、`chem.protein.SOLVENT_AND_IONS`として`chem.protein`パッケージレベルでも再エクスポートされる(`find_pocket`内部だけでなく、notebookでの表示用にリガンドらしきHETATMを判定する用途にも再利用できるよう公開)。

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- 新規依存: `biopython`(PyPI版あり、`pyproject.toml`の`dependencies`に追加)、`numpy`(同様)、`fpocket`(PyPIなし、`environment.yml`のconda-forge依存に追加。`mamba install -n chem -c conda-forge fpocket`または`mamba env update -f environment.yml`でインストール)
- `gemmi`は既存環境に(他パッケージ経由で)入っているが、今回は使わずBio.PDBのみで実装する

## サンプルノートブック

`notebooks/example_proteins.ipynb`の末尾に、`chem.rcsb`でダウンロードした複数のトロンビン構造をアラインし、py3Dmolで重ね書き表示するサンプルセルを追加する。fpocketのサンプルは別途、リガンド入り構造(例: トリプシン+ベンザミジン `3PTB`)で`find_pocket`を実行し、選ばれたポケットの残基をハイライト表示するセルを追加する。

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- テストは`tests/test_protein_align.py`・`tests/test_protein_pocket.py`にネットワーク・外部バイナリ(fpocket)不要な範囲(シーケンスマッチングロジック、リガンド自動検出・HETコード判定・ファイル判定の分岐、ポケット選択の距離計算、`_info.txt`パーサ、fpocket未インストール時のエラーメッセージ)のみ追加する。実際のBio.PDB構造アラインメントとfpocket実行は、簡易的な合成PDBテキスト(固定カラム位置で手書きしたATOM/HETATMレコード)を使ったオフラインテストと、実データでの手動実行確認(3PTB+BENでPocket 1が選ばれ、Asp189・Ser195が残基リストに含まれることを確認済み)の組み合わせで検証する
