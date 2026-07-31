# chem.protein.align / chem.protein.find_pocket / chem.protein.list_pockets の再現プロンプト

`chem`リポジトリの`chem`パッケージ内に、構造解析用のサブパッケージ`protein`(`chem.protein`)を追加し、複数構造のシーケンス・3Dアラインメント(`align`)と、リガンド近傍のfpocketポケット・構成残基特定(`find_pocket`)、リガンドを持たない構造向けの全候補ポケット列挙(`list_pockets`)を実装するための指示。`chem.rcsb`・`chem.alphafold`でダウンロードした構造をそのまま入力にできることを想定する。

## パッケージ構成

- `src/chem/protein/__init__.py`: `from .structural_align import align` / `from .pocket import find_pocket, list_pockets` として再エクスポート
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
- `chain`: 各非reference構造で使うチェーンIDを明示指定(省略時は下記の自動選択)。referenceの自身のチェーンにも適用される
- `outdir`: 座標変換後の構造を書き出す先のディレクトリ(存在しなければ作成)
- 戻り値: `{path: {"rmsd": ..., "identity": ...}}`(シーケンスマッチしたCA原子上のRMSDと配列一致度。reference自身は`{"rmsd": 0.0, "identity": 1.0}`)

### アルゴリズム

1. `reference`を解決し、Bio.PDB(拡張子`.cif`/`.mmcif`なら`MMCIFParser`、それ以外は`PDBParser`)で読み込む
2. referenceの「主鎖(primary polymer chain)」を選択する: `chain`引数が指定されていればそのチェーンID(見つからなければ`ValueError`)、なければ標準アミノ酸残基を最も多く含むチェーンを自動選択する(`_select_chain`。比較対象がreference自身しかないため、後述のmobile構造向けとは違いサイズだけで選ぶ)。**単純に「ファイル内最初のチェーン」を採用してはいけない** — 例えばトロンビンはL鎖(軽鎖、~36残基)がH鎖(重鎖、~259残基)より先に現れることがあり、「最初に見つかった一定長以上のチェーン」という閾値ベースの選択だと軽鎖を誤って選んでしまう
   - 「標準アミノ酸残基」の判定は`Bio.PDB.Polypeptide.is_aa(residue, standard=True)`**だけでは不十分**: `is_aa`はresnameのみを見ており、hetero flag(`residue.id[0]`)をチェックしない。共有結合したペプチド模倣リガンド(D-アミノ酸を含むペプチド性阻害剤など)が蛋白質と同じチェーンIDを使い、かつ標準アミノ酸名(例: `PRO`)をHETATMとして持つ場合、`is_aa`だけの判定だとリガンド側の残基を誤って蛋白質配列に取り込んでしまう。実データ(トロンビン+ペプチド性阻害剤複合体 `6YHG`、リガンド側の`PRO H 307`がHETATMとして紐づいていた)でこの誤りによりRMSDが約0.3Å→約1.9Åに悪化する不具合を発見・修正済み。`_is_polymer_residue(r) = r.id[0] == " " and is_aa(r, standard=True)`という追加のhetero flagチェックを`_select_chain`・`_chain_seq_and_ca`の両方に適用する
3. `reference`のチェーンから、CA原子を持つ標準アミノ酸残基(上記の`_is_polymer_residue`基準)のみを使って1文字シーケンスと対応するCA原子リストを作る(`ref_seq`/`ref_ca`)
4. 各非reference構造について、使うチェーンを選ぶ:
   - `chain`引数が指定されていればそのチェーンID(`_select_chain`)をそのまま使い、`_matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca)`を1回呼ぶ
   - 指定がなければ`_best_matching_chain_alignment(model, ref_seq, ref_ca)`を使う: モデル内の(標準アミノ酸残基を持つ)全チェーンそれぞれについて`_chain_seq_and_ca`→`_matched_ca_pairs`を実行し、`identity * マッチ位置数`(=一致した残基の実数)が最大のチェーンの`(ref_pts, mob_pts, identity)`を採用する。**なぜ「残基数最大のチェーン」ではダメか**: 実データ(トロンビン+Protein C Inhibitor複合体`3B9F`)で、セルピン阻害剤のチェーン(356残基)がトロンビン自身の重鎖(253残基)より残基数が多く、単純なサイズ比較だと誤って阻害剤側のチェーンを選んでしまう不具合を発見・修正済み。しかも「マッチ位置数(ギャップなし長)」で比較しても解決しない — `open_gap_score=-10`が`mismatch_score=-1`よりずっと重いため、無関係な配列同士でもアライナーはギャップを開かずミスマッチのまま繋げようとし、実際このケースでは阻害剤チェーンの342残基がギャップなしで"マッチ"扱いになった(トロンビン重鎖の253残基より多い)。一致度(identity)で見ると阻害剤チェーンは0.219、トロンビン重鎖は0.996と歴然の差があり、これが正しい判定基準になる
5. (`_matched_ca_pairs`内)`Bio.Align.PairwiseAligner`(`mode="global"`、`open_gap_score=-10`、`extend_gap_score=-0.5`、`match_score=2`、`mismatch_score=-1`)でreferenceシーケンスと対象シーケンスをグローバルアラインメントし、`alignment.aligned`からギャップのないマッチ位置ペアを取り出し、対応するCA原子ペアのリストを作る。同時に、各マッチ位置についてreference側と対象側のアミノ酸1文字が一致するかどうかも数え、「一致数 / マッチ位置数」を`identity`として返す(戻り値は`(ref_pts, mob_pts, identity)`の3-tuple)。**注意**: ここでの「マッチ」はギャップを挟まない位置という意味であり、アミノ酸が同一という意味ではない。つまりギャップ(挿入・欠損)になった位置だけが`identity`の分母・分子どちらからも除外され、変異(置換)はギャップにならない限り分母には数えられるが分子には数えられない
6. マッチしたCA原子ペアが3組未満なら(`_MIN_MATCHED_RESIDUES = 3`)、その構造はアラインメント不可としてスキップ(quiet時以外は`stderr`に警告を出し、処理は継続)
7. `Bio.PDB.Superimposer`でマッチしたCA原子ペアに対してKabsch法によるフィッティングを行い、得られた回転・並進を**構造の全原子**(タンパク質だけでなくリガンド・水も含む)に適用する。これにより、align出力後もリガンドの相対位置が保たれ、後続の`find_pocket`にそのまま使える。`results[path] = {"rmsd": round(float(sup.rms), 3), "identity": round(identity, 3)}`
8. `Bio.PDB.PDBIO`で、referenceも含めて各構造を`outdir`に`{入力ファイルのstem}.pdb`として書き出す(**入力がCIFでも常にPDB形式で出力する** — fpocketがPDB形式を要求するための一貫性)。referenceの書き出しはメインループの前に1回だけ行い、`structures`側のループでは`path == ref_path`のときスキップする(referenceが`structures`に含まれない場合でも正しく書き出される)

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
  - `spheres`: `[{"x":..., "y":..., "z":..., "radius":...}, ...]`(fpocketが空洞の形状を近似するために生成したアルファ球の中心座標・半径。`residues`が「ポケットを構成する蛋白質原子」であるのに対しこちらは「空洞そのものの形状」。py3Dmolの`addSphere`を球ごとに1回呼べば、リガンド近傍を構成残基のstickではなく充填された体積として表示できる)
  - `info`: 該当ポケットのfpocketスコア情報全体(`_info.txt`の全フィールドを`{項目名: 値}`の辞書にしたもの)

### アルゴリズム

1. `structure`を作業ディレクトリ(`outdir`があればそこ、なければ一時ディレクトリ)にコピーし、`fpocket -f {basename}`をサブプロセスとして実行する(`cwd`を作業ディレクトリにすることで、元ファイルの場所を汚さない)。`fpocket`実行ファイルが見つからない場合はconda-forgeでのインストール方法を含む分かりやすい`RuntimeError`にする
2. `{stem}_out/pockets/pocket{N}_atm.pdb`(各ポケットを構成する蛋白質原子)、同ディレクトリの`pocket{N}_vert.pqr`(fpocketが空洞形状を近似するために生成したアルファ球群。PQR形式で1行1球、列は`serial, atom名, resname("STP"), resseq, x, y, z, charge, radius`)、`{stem}_out/{stem}_info.txt`(ポケットごとのスコア情報)がfpocketの出力として得られる
3. `ligand`引数からリガンド原子の3D座標(numpy配列)を解決する(上記参照)
4. 各`pocket{N}_atm.pdb`をBio.PDBでパースし、その原子群とリガンド原子群との最小距離(全原子対distanceのmin)を計算する。最小距離が最も小さいポケット番号を採用する(3PTB(トリプシン)+BENリガンドで検証済み: 最もfpocketスコアが高いPocket 1が同時に最もBENに近く、既知の活性部位Asp189・Ser195を含むことを確認)
5. 採用したポケットについて、`_pocket_result(pocket_id, atm_path, info)`ヘルパーで結果辞書を組み立てる: `pocket{N}_atm.pdb`から一意な`(chain, resnum, resname)`のリスト(`residues`)を抽出し、`atm_path`と同じディレクトリの`pocket{N}_vert.pqr`(`_atm.pdb`→`_vert.pqr`と文字列置換したパス)を`_parse_pocket_spheres`でパースして`spheres`とする(ファイルが存在しなければ空リスト)
6. `{stem}_info.txt`をパースする。フォーマットは`Pocket N :`という見出し行の後に、タブ区切りの`項目名 : \t値`行が続き、空行でポケットが区切られる。数値に変換できるものは`float`にする

## chem.protein.list_pockets

```python
from chem import protein
protein.list_pockets(structure, outdir=None, druggability_thres=0.1)
```

AlphaFold予測構造のようにリガンドが一切結合していない構造に対して、`find_pocket`(リガンド近傍の1つだけを選ぶ)の代わりに、fpocketが検出した候補ポケットを一覧で返す。

- `structure`: PDBファイルのパス(`find_pocket`と同様)
- `outdir`: fpocketの生ログ出力を保持するディレクトリ。`None`なら一時ディレクトリを使い処理後に破棄する(`find_pocket`と同様)
- `druggability_thres`: 保持する`druggability_score`の下限(inclusive)。fpocketは典型的な構造で数十個の低品質な空洞をほぼ0点で報告するため、デフォルト(`0.1`)でそれらと`druggability_score`が全く付与されなかったポケットを除外する。`None`を渡せばフィルタなしで全ポケットを返す
- 戻り値: `find_pocket`の戻り値と全く同じ形の辞書(`pocket_id`/`score`/`druggability_score`/`volume`/`residues`/`spheres`/`info`)の、保持されたポケットのみのリスト。**`druggability_score`降順**でソートする

### アルゴリズム

1. `find_pocket`と同じ`_run_fpocket`でfpocketを実行し、`{stem}_info.txt`を`_parse_info_txt`でパースする(リガンド解決は行わない)
2. `pockets/pocket{N}_atm.pdb`をすべて(`_pocket_atm_paths`でid→パスの辞書として)走査し、各ポケットについて`find_pocket`と共通の`_pocket_result(pocket_id, atm_path, info)`ヘルパーで結果辞書を組み立てる(`find_pocket`もこのヘルパーを使うようリファクタリングする)
3. `druggability_thres`が`None`でなければ、`druggability_score`が`None`のポケット、および`druggability_thres`未満のポケットを除外する
4. `druggability_score`降順(`druggability_thres=None`で`None`が残った場合は最後)でソートして返す

### `SOLVENT_AND_IONS`(自動検出で除外するHETコード、非網羅的)

```
HOH, WAT, DOD,
NA, CL, K, MG, CA, ZN, MN, FE, FE2, CO, NI, CU, CD, LI, RB, CS, BR, IOD, NH4,
SO4, PO4, GOL, EDO, PEG, PG4, 1PE, P6G, MPD, FMT, ACT, DMS, TRS, BME, EPE, HEPES, IPA, UNX, UNL
```

`chem/protein/pocket.py`のモジュールレベル定数で、`chem.protein.SOLVENT_AND_IONS`として`chem.protein`パッケージレベルでも再エクスポートされる。`find_pocket`のリガンド自動検出、`chem.ligand.list_ligand_codes`/`list_ligand_instances`のデフォルト`exclude`で使う — いずれも「これは薬理学的に意味のあるリガンドか」を判定する用途なので、結合イオンや結晶化添加剤も除外する

### `WATER`(表示用、水のみ除外)

```
HOH, WAT, DOD
```

同じく`chem/protein/pocket.py`のモジュールレベル定数(`SOLVENT_AND_IONS`の真部分集合)で、`chem.protein.WATER`として再エクスポートされる。`chem.view3d.render_protein`のデフォルト`exclude`はこちら — 「表示する価値があるか」という基準では、結合イオン(例: `ZN`、`NA`)や結晶化添加剤も薬らしくはないが構造の一部として見えていた方が有用なことが多く、`SOLVENT_AND_IONS`ほど広く除外すべきではないため、水だけを除く別の定数として分けている

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- 新規依存: `biopython`(PyPI版あり、`pyproject.toml`の`dependencies`に追加)、`numpy`(同様)、`fpocket`(PyPIなし、`environment.yml`のconda-forge依存に追加。`mamba install -n chem -c conda-forge fpocket`または`mamba env update -f environment.yml`でインストール)
- `gemmi`は既存環境に(他パッケージ経由で)入っているが、今回は使わずBio.PDBのみで実装する

## サンプルノートブック

`notebooks/alphafold_pocket_thrb_human.ipynb`のRCSBダウンロードセルの直後に、ダウンロードした複数のトロンビン構造をAlphaFold予測構造を参照にアラインする独立セクション(タイトルmd + コード。コードは`align()`実行と`align_df`(rmsd/identity列)の表示のみ)を置く。重ね書きpy3Dmolビューア(トグルボタンで表示構造を選ぶウィジェット)はさらに別の独立セクション(独自のH2タイトルmd + コード)として、アラインメントセクションの直後に続ける -- 1セルに両方を詰め込まない。fpocketのサンプルは別途、リガンド入り構造(例: トリプシン+ベンザミジン `3PTB`)で`find_pocket`を実行し、選ばれたポケットの残基をハイライト表示するセルを追加する。`list_pockets`のサンプルは、AlphaFold予測構造セクション(リガンドが存在しない)に以下3セルを追加する: (1) `pandas.DataFrame`で「pocket_id / score / druggability_score / volume / n_residues」の表として全候補ポケットを表示、(2) `druggability_score >= 0.2`の候補ポケットを、半透明cartoonの上にポケットごとに異なる色の構成残基stickでハイライトして可視化(色とpocket_id/druggability_scoreの凡例付き)、(3) 同じ候補ポケットを、構成残基のstickの代わりに`spheres`フィールド(fpocketのアルファ球)を`py3Dmol.addSphere`で球ごとに描画し、空洞を充填された体積として可視化(こちらもポケットごとに色分け・凡例付き)。既存セルは書き換えず、新規セルとして追記する。

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- テストは`tests/test_protein_align.py`・`tests/test_protein_pocket.py`にネットワーク・外部バイナリ(fpocket)不要な範囲(シーケンスマッチングロジック、`_matched_ca_pairs`が返す`identity`について同一配列で`1.0`・ギャップを含む場合はギャップ位置を分母/分子どちらからも除外・ギャップなしミスマッチを含む場合は分母に数えて分子には数えないこと、`align()`の戻り値が`{"rmsd":..., "identity":...}`の形でreferenceは`{"rmsd": 0.0, "identity": 1.0}`になること、`_best_matching_chain_alignment`が「マッチ位置数は多いが一致度が低い大きなチェーン」より「短くても一致度が高いチェーン」を選ぶこと(`3B9F`相当の合成データで再現)、`align()`をエンドツーエンドで実行してもサイズ最大の無関係なチェーンではなく正しいチェーンが選ばれ`identity`が高くなること、リガンド自動検出・HETコード判定・ファイル判定の分岐、ポケット選択の距離計算、`_info.txt`パーサ、fpocket未インストール時のエラーメッセージ、`_pocket_atm_paths`/`_pocket_result`の単体動作(`pocket{N}_vert.pqr`が存在する/しない両方のケース)、`_parse_pocket_spheres`が合成PQRテキストから`x`/`y`/`z`/`radius`を正しく取り出すこと、`list_pockets`を`_run_fpocket`を`monkeypatch`で偽の出力ディレクトリに差し替えて実行し`druggability_thres=None`なら全ポケットが`druggability_score`降順(値なしは最後)で返ること・デフォルト(`0.1`)では値なし/閾値未満のポケットが除外されること・`druggability_thres`を明示指定すればその閾値で絞り込まれること、`WATER`が`SOLVENT_AND_IONS`の真部分集合であること)のみ追加する。実際のBio.PDB構造アラインメントとfpocket実行は、簡易的な合成PDBテキスト(固定カラム位置で手書きしたATOM/HETATMレコード)を使ったオフラインテストと、実データでの手動実行確認(3PTB+BENでPocket 1が選ばれAsp189・Ser195が残基リストに含まれること、AlphaFold予測構造(リガンド無し)で`list_pockets`が複数候補を返し、そのうち`druggability_score >= 0.2`の3件を`spheres`経由でHTMLに書き出しブラウザで表示確認したところ、cartoon上にポケットごとに色分けされた充填体積として描画されることを確認済み)の組み合わせで検証する
