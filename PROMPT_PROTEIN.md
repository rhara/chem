# chem.protein.summary / chem.protein.get_fasta / chem.protein.sequence_align / chem.protein.align / chem.protein.compute_transform / chem.protein.apply_transform / chem.protein.identity_matrix / chem.protein.find_pocket / chem.protein.list_pockets / chem.protein.residues_near / chem.protein.split の再現プロンプト

`chem`リポジトリの`chem`パッケージ内に、構造解析用のサブパッケージ`protein`(`chem.protein`)を追加し、UniProtエントリのアノテーション取得(`summary`)とFASTA配列取得(`get_fasta`)、複数のPDB/CIF構造の観測配列をUniProt canonical配列の座標系に1次元アラインメントする(`sequence_align`)、複数構造のシーケンス・3Dアラインメント(`align`)、そのアラインメントの回転・並進だけを取り出して他の構造ファイルに使い回す(`compute_transform`/`apply_transform`)、任意の構造集合の総当たり配列一致度マトリクス(`identity_matrix`)、リガンド近傍のfpocketポケット・構成残基特定(`find_pocket`)、リガンドを持たない構造向けの全候補ポケット列挙(`list_pockets`)、リガンドファイルからの単純な距離ベース活性部位残基抽出(`residues_near`)、構造ファイルをリガンドフリーの蛋白質PDBと各HETATM残基のSDFに分割する(`split`)を実装するための指示。`chem.rcsb`・`chem.alphafold`でダウンロードした構造をそのまま入力にできることを想定する。

## パッケージ構成

- `src/chem/protein/__init__.py`: `from .annotation import get_fasta, summary` / `from .sequence_align import sequence_align` / `from .structural_align import align, apply_transform, compute_transform, identity_matrix` / `from .pocket import SOLVENT_AND_IONS, WATER, find_pocket, list_pockets, residues_near` / `from .splitting import split` として再エクスポート
- `src/chem/protein/annotation.py`: `summary`・`get_fasta`の実装
- `src/chem/protein/sequence_align.py`: `sequence_align`の実装(UniProt REST APIへの`requests.get`、`Bio.Align.PairwiseAligner`によるcanonicalへの1次元アラインメント、ProDyによる構造からのCA配列抽出を1ファイルにまとめる。`structural_align.py`とは別ファイルにする — こちらはBio.PDBではなくProDyでCA配列を読み、3D重ね合わせは一切行わない全く別の処理系統のため)
- `src/chem/protein/structural_align.py`: `align`・`compute_transform`・`apply_transform`・`identity_matrix`の実装(いずれも`_load_structure`/`_select_chain`/`_chain_seq_and_ca`/`_matched_ca_pairs`を共有するため同じファイルに置く。`compute_transform`は`align`のアラインメント〜Kabsch fitまでを再利用し、構造の書き出しだけをしない版。`apply_transform`はその結果の`rotation`/`translation`を任意の構造ファイルに適用するだけの薄い関数)
- `src/chem/protein/pocket.py`: `find_pocket`・`list_pockets`・`residues_near`の実装(`residues_near`は`find_pocket`が使う`_load_external_ligand_coords`をリガンド座標読み込みに再利用する)
- `src/chem/protein/splitting.py`: `split`の実装(ファイル名を`split.py`にしない理由は下記の注意点と同じ)

**注意**: 実装ファイル名を公開関数名と同じにしない(例: `align.py`に`align`関数を置かない、`summary.py`に`summary`関数を置かない)。パッケージの`__init__.py`で`from .align import align`のように再エクスポートすると、`chem.protein`パッケージの`align`属性がサブモジュールからその関数へ上書きされてしまい、以降`import chem.protein.align`のようなドット区切りでのサブモジュールアクセスが(内部テストなどで)壊れる。既存の`chembl`/`rcsb`/`alphafold`サブパッケージが実装ファイルを常に`fetch.py`(公開関数名と別名)にしているのも同じ理由。

## chem.protein.summary

```python
from chem import protein
protein.summary(id)
```

`chem.rcsb`/`chem.alphafold`で構造をダウンロードしたり`chem.blast.blastp`でホモログを探したりする前段として、対象タンパク質そのものの素性を素早く確認するための、UniProtエントリのアノテーション取得関数。

- `id`: UniProtアクセッション(例: `"Q8IZL9"`)、エントリ名/mnemonic(例: `"CDK20_HUMAN"`)、またはChEMBL target id(例: `"CHEMBL3559690"`)。`chem.ids.resolve_uniprot_accession_any`でUniProtアクセッションに解決してから使う
- 戻り値: 以下のキーを持つ`dict`(UniProtに該当データがない項目は`None`、`has_*`系は`False`になる。例外は投げない)
  - `entry_name`, `accession`: UniProtのエントリ名(mnemonic)とアクセッション
  - `protein_name`, `gene_name`(synonymsがあれば`"CDK20 (synonyms: CCRK, CDCH)"`のように付記), `organism`, `sequence_length`, `ec_number`
  - `family`: `SIMILARITY`コメント(蛋白質ファミリー分類)
  - `function`: `FUNCTION`コメント
  - `subcellular_location`: `SUBCELLULAR LOCATION`コメントの各locationをカンマ結合
  - `kinase_domain_range`: `description`に"kinase"(大文字小文字無視)を含む最初の`Domain`フィーチャーの`"{start}-{end}"`残基範囲(非キナーゼやドメイン未注釈のエントリでは`None`)
  - `active_site_residue`: 最初の`Active site`フィーチャーの残基番号(未注釈なら`None`)
  - `n_pdb_xrefs`: UniProtエントリ内の`PDB`クロスリファレンス数(RCSB Search APIによるライブな件数ではない点に注意 — 実際のPDB構造数が必要なら`chem.rcsb`側を使う)
  - `has_alphafold_model`, `has_bindingdb_entry`: `AlphaFoldDB`/`BindingDB`クロスリファレンスの有無
  - `chembl_target_id`: `ChEMBL`クロスリファレンスのid(なければ`None`)
  - `pharos_development_level`: Pharosのtarget development level(`Tclin`/`Tchem`/`Tbio`/`Tdark`)に簡単な説明を付記した文字列(例: `"Tbio (biology characterized, no known drug/chemical probe)"`)。UniProtに`Pharos`クロスリファレンスがなければ`None`
  - `protein_existence`, `annotation_score`: UniProt自身のエビデンスレベル・アノテーション充実度スコア

### アルゴリズム

1. `id`を`chem.ids.resolve_uniprot_accession_any`でUniProtアクセッションに解決する(ChEMBL target id/UniProtアクセッション/エントリ名のいずれも受け付ける)
2. `https://rest.uniprot.org/uniprotkb/{accession}.json`をfetchする
3. 上記の各プロパティをJSONの該当フィールド(`proteinDescription`/`genes`/`comments`/`features`/`uniProtKBCrossReferences`等)から抽出し、辞書にまとめて返す(抽出ロジックは`_extract_properties(entry)`という、ネットワークを叩かない純粋関数に分離する — テストではこちらに直接、手組みのJSON相当の`dict`を渡す)

## chem.protein.get_fasta

```python
from chem import protein
protein.get_fasta(id, email="user@example.com")
```

UniProtエントリの配列をFASTA文字列として取得する。`summary`と同じくID解決を共有する軽量な関数。

- `id`: `summary`と同じ3形式(UniProtアクセッション/エントリ名/ChEMBL target id)を受け付ける
- `email`: 省略可能。UniProt REST APIは必須にしていないが、[UniProt自身のAPI利用ガイドライン](https://www.uniprot.org/help/api)がリクエストの`User-Agent`ヘッダーに連絡先メールアドレスを含めることを推奨しているため、`User-Agent: chem/{chem.__version__} ({email})`として送る。デフォルトは`"user@example.com"`のプレースホルダー
- 戻り値: `https://rest.uniprot.org/uniprotkb/{accession}.fasta`のレスポンステキストそのまま(FASTA形式の文字列)

## chem.protein.sequence_align

```python
from chem import protein

result = protein.sequence_align(
    accession,
    structures,
    canonical_feature_type=None,
    canonical_feature_description=None,
    marker_feature_types=None,
    chain=None,
)
```

`align`/`identity_matrix`は同じターゲットの構造同士を互いに比較する(reference構造、または構造集合内の他の構造)が、`sequence_align`はUniProtの**canonical配列**(WTの完全長配列、または特定のfeature範囲)を唯一の基準にして、複数構造の観測配列(構造中に実際に見えている残基そのもの)を全て同じ座標系に載せる。WT構造・変異体構造・電子密度が部分的にしか見えていない構造を横並びで比較したい用途に向く。3D重ね合わせ(Kabsch fit)や構造ファイルの書き出しは一切行わない — アラインメント済みのデータを返すだけで、そこからテキスト表示・色付きHTML・変異点一覧などを組み立てるのは呼び出し側の仕事とする(「データはライブラリ、表示はカスタムスクリプト」という設計方針。リポジトリ所有者からの明示的な指示: 「出力はカスタムなスクリプトで行うという設計にしたい」)。

- `accession`: UniProtアクセッション(例: `"P07550"`)
- `structures`: PDB/CIF構造ファイルパスのリスト(`chem.rcsb.download_structures`の出力をそのまま渡せる。`align`と同じく、IDではなくファイルパスを受け取る設計 — ダウンロードは`chem.rcsb`側の責務という既存方針を踏襲する)
- `canonical_feature_type` / `canonical_feature_description`: canonical配列として切り出すUniProt featureを指定する。例: シグナルペプチドが切断される分泌蛋白なら`("Chain", "Beta-lactamase TEM")`のように成熟鎖のfeatureを指定してその範囲だけを使う、多ドメイン蛋白の1ドメインだけを見たいなら`("Domain", "Bromo 1")`のように指定する。`canonical_feature_type=None`(デフォルト)ならUniProtのfull-length配列をそのまま使う(`canonical_feature_description`だけを指定してもtypeがNoneなら無視される)
- `marker_feature_types`: マーカーとして集めたいUniProt feature typeのタプル。酵素の触媒残基なら`("Active site",)`、受容体のリガンド結合ポケットなら`("Binding site",)`、非触媒だが機能的に重要な残基なら`("Site",)`など、対象蛋白の性質に応じて呼び出し側が選ぶ。`None`(デフォルト)なら空集合(マーカーなし)を返す — 全ての蛋白が触媒残基を持つわけではない(受容体やリーダードメインなど)ため、決め打ちの単一feature typeにしない
- `chain`: 全構造で使うチェーンIDを明示指定(省略時は下記の自動選択)
- 戻り値: 以下のキーを持つ`dict`
  - `protein_name`, `organism`: UniProtエントリから
  - `canonical_seq`: (featureで切り出した場合はその範囲の)canonical配列
  - `feature_start`, `feature_end`: `canonical_seq`のUniProt全長配列内での1-based範囲(`canonical_feature_type=None`なら`(1, len(full_seq))`)
  - `marker_positions`: `marker_feature_types`から集めた、`canonical_seq`内(1-based)の位置集合。`marker_feature_types=None`なら空集合
  - `raw_sequences`: `{path: 配列}`。各構造で選ばれたチェーンのCA由来配列を、アラインメント前のそのまま
  - `sequences`: `{path: 配列}`。`raw_sequences`の各配列をcanonicalの座標に収めたもの — 全ての`path`について`len(sequences[path]) == len(canonical_seq)`が成り立ち、その構造でcanonicalのその位置が観測されていなければ`-`になる

### アルゴリズム

1. `requests.get(f"{UNIPROT_API}/{accession}.json")`でUniProtエントリを取得する(`chem.ids`のようなID解決は行わない — `accession`は既にUniProtアクセッションであることを前提とする)
2. `canonical_feature_type`が`None`なら`full_seq`全体、そうでなければ`entry["features"]`から`type == canonical_feature_type`かつ(`canonical_feature_description`が指定されていれば)`description == canonical_feature_description`に一致する最初のfeatureの`location.start/end`で`full_seq`をスライスし、`canonical_seq`・`feature_start`・`feature_end`とする
3. `marker_feature_types`が`None`または空なら`marker_positions = set()`。そうでなければ`entry["features"]`から`type in marker_feature_types`かつ`feature_start <= location.start.value <= feature_end`のfeatureを集め、各`location.start.value - feature_start + 1`(canonical内1-based相対位置)の集合を返す
4. `Bio.Align.PairwiseAligner`を1つ構築する(`mode="global"`、`open_gap_score=-10`、`extend_gap_score=-0.5`、`match_score=2`、`mismatch_score=-3`)。**`mismatch_score`の値は慎重に選ぶ必要がある**:
   - `-1`(Biopythonのデフォルトに近い、ごく僅かなペナルティ)だと、T4リゾチームのような長い(~160残基)融合パートナーがcanonical上のある範囲(例: 細胞内ループ3、~30残基)を置き換えている場合、正しく1つの連続したギャップにならず、接合部で数残基が「たまたま似ている」候補にミスマッチとして吸収され、本来ならギャップになるべき位置に偽の変異が大量に出現する(実データ: β2アドレナリン受容体の複数のT4L融合構造で確認・修正した不具合)
   - 対策として`mismatch_score`を`open_gap_score`と同値(`-10`)まで引き上げると、この融合部分のスミアリングは解消するが、**別の、より深刻な不具合が生じる**: 既に開いている(隣接する)ギャップのすぐ隣に本物の点変異がある場合、そのギャップを1残基分延長するコストは`extend_gap_score`(`-0.5`)しかかからないため、`mismatch_score`をそこまで引き上げると本物の変異(例: `S262D`、`C265F`)がミスマッチとして報告される代わりにギャップへ静かに吸収されて消えてしまう(実データのβ2AR比較で発見。`open_gap_score`と同値にした結果、複数構造で既知の変異点が消失することを確認して`-3`に差し戻した)
   - `-3`は、実際に扱う規模(数十残基以上の融合ブロック)でのスミアリング防止と、ギャップ隣接の点変異を保持することの両立点として選んだ値。原理的な限界として、canonicalで置き換えられる範囲が非常に短い(4残基程度以下)場合はなお僅かにスミアリングし得るが、実際に検証した構造(TEM-1 β-ラクタマーゼ、BRD4、β2ARのT4リゾチーム/ナノボディ/Gタンパク質融合)では十分機能する
5. 各`structures`について、ProDy(`prody.parsePDB`、拡張子が`.cif`/`.mmcif`なら`prody.parseMMCIF`)で読み込み、`prody.HierView`で全ポリマー鎖の`chain.select("protein and name CA")`から1文字配列を集める(`_load_ca_sequences`)。これはATOMレコードに実際に座標を持つ残基のみで、SEQRES(構築物全体の設計配列)は一切見ない — 「生のPDBファイルで実際に見えている残基」をそのまま反映する
6. チェーン選択: `chain`引数が指定されていればそのチェーンID(無ければ`ValueError`)。指定が無ければ、その構造の全チェーンそれぞれについて`aligner.align(canonical_seq, seq).score`を計算し、最もスコアの高いチェーンを採用する(`_best_matching_chain_sequence`)。**単純な「chain A」固定ではない** — GPCR構造は同じchain IDでも構造ごとにT4リゾチーム融合・ナノボディ・ヘテロ三量体Gタンパク質サブユニットなど、標的蛋白以外のポリマー鎖が混ざっており、これらを取り違えないための自動選択が必須になる
7. 選ばれた生配列を`canonical_seq`にアラインメントし、canonical座標に収める(`_align_to_canonical`)。`aligner.align(canonical_seq, query_seq)`は同点最適解を複数返し得るため(下記参照)、`itertools.islice`で最大`_MAX_TIED_ALIGNMENTS=500`個まで列挙し、その中から`len(alignment.aligned[0])`(マッチしたブロック数)が最小のものを選ぶ。選んだアラインメントの`gapped_canonical`/`gapped_query`から、`gapped_canonical`が`-`でない列だけを残して(`gapped_query`側の対応する文字、`-`ならそのまま)連結したものを返す — これにより長さが必ず`len(canonical_seq)`と一致し、query側がcanonicalに無い挿入(融合パートナーなど)を持っていた場合はその分の文字がそのまま失われる(canonical上に対応する位置が無いため)
8. **同点タイブレークの理由**: canonicalの構築物に含まれない長い未observed区間(例: β2ARのC末端側~70残基の細胞内ドメイン、常にどの結晶構造でも解けない)の直前で、query配列の最後の残基がcanonicalのその区間より**後**にある残基と偶然同じアミノ酸だった場合(例: 両方とも`L`)、「区間の直前で正しく終わる」アラインメントと「区間の直後まで飛んで(間に大きなギャップを開けて)そこにマッチさせる」アラインメントが完全に同スコアになり得る。Biopythonはこの場合どちらを返すか保証しない。実データ(2RH1)で、本来canonical位置342で終わるはずの観測配列が、たまたま同じ`L`である位置413(配列の本当の末尾)に誤ってマッチしていたことを、PDBファイル自身の残基番号(`resnum`が29から342まで連続していることを直接確認)で検証し発見した。マッチブロック数最小のタイブレークはこれを正しい側(canonicalの342で終わる、1つの連続したギャップ)に解決する。ただし**それでも解決しない残るケース**がある: canonical自身が同一アミノ酸の繰り返しで終わる(または始まる)場合(例: `...NDSLL`の`L,L`)、query側の対応する残基数がその繰り返しの長さより短ければ、「末尾1文字だけ一致」と「末尾2文字とも一致」が完全に同スコア・同ブロック数になり得、これは配列だけからは原理的に決定不能(6MXTで確認)。実害は末端の`del`範囲が±1残基ずれる程度に限られる

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

## chem.protein.compute_transform / chem.protein.apply_transform

```python
from chem import protein

t = protein.compute_transform(mobile, reference, mobile_chain=None, reference_chain=None)
protein.apply_transform(structure_path, t["rotation"], t["translation"], outpath)
```

`align`は複合体の**全チェーン**を1つの構造ファイルとしてまとめて重ね合わせるが、`chem.protein.split`が書き出したような「チェーンごとに1ファイル」のセットでは、キナーゼ本体はシーケンスアラインメントできても、結合パートナー(サイクリンなど)やリガンドSDFには比較対象となる配列そのものが無い。`compute_transform`はキナーゼ側1本のアラインメント〜Kabsch fitの結果(回転・並進)だけを取り出して返し、`apply_transform`(蛋白質PDB用)/`chem.ligand.apply_transform`(リガンドSDF用)でそのまま他ピースに使い回すことで、複合体全体を「あたかも一括でアラインしたかのように」共通座標系に再構成できるようにする。

- `mobile`, `reference`: PDB/CIFファイルパス
- `mobile_chain`, `reference_chain`: 使うチェーンIDを明示指定(省略時は`align`と同じ自動選択 — referenceはサイズ最大の主鎖、mobileはreferenceとの一致度が最も高いチェーン)
- 戻り値: `{"rotation": 3x3のネストリスト, "translation": 3要素リスト, "rmsd": float, "identity": float}`。`Bio.PDB.Superimposer`/`Atom.transform`の規約(`new_coord = old_coord @ rotation + translation`、行ベクトル座標)そのまま。`rmsd`はシーケンスマッチしたCA原子上のÅ単位RMSD、`identity`は`align`と同じ定義。マッチ残基が`_MIN_MATCHED_RESIDUES`未満なら`ValueError`

### アルゴリズム(`compute_transform`)

1. `align`のステップ1〜7と全く同じロジック(`_load_structure`→チェーン選択→`_chain_seq_and_ca`→(`mobile_chain`指定時は`_matched_ca_pairs`を直接、省略時は`_best_matching_chain_alignment`)→`Bio.PDB.Superimposer.set_atoms`)を実行するが、**構造ファイルへの適用・書き出しは一切行わない**
2. `sup.rotran`から得た`(rot, tran)`を`{"rotation": rot.tolist(), "translation": tran.tolist(), "rmsd": round(float(sup.rms), 3), "identity": round(identity, 3)}`として返す

`apply_transform(structure_path, rotation, translation, outpath)`:

- `structure_path`: 変換対象のPDB/CIFファイル(`compute_transform`の計算に使ったファイルそのものである必要はない — 例えば`chem.protein.split`が書き出した、キナーゼとは別チェーンのパートナー蛋白質PDB)
- `rotation`, `translation`: `compute_transform`の戻り値の`"rotation"`/`"translation"`をそのまま渡す
- `outpath`: 出力PDBファイルパス(親ディレクトリが無ければ作成)
- 戻り値: `outpath`

### アルゴリズム(`apply_transform`)

1. `_load_structure(structure_path)`で読み込み、`rotation`/`translation`を`numpy`配列(`dtype="f"`)に変換
2. 構造内の**全原子**(`structure.get_atoms()`)に対して`atom.transform(rot, tran)`(Biopythonの規約通り`new_coord = old_coord @ rotation + translation`)を適用する
3. `Bio.PDB.PDBIO`で`outpath`に書き出す(**入力がCIFでも常にPDB形式で出力する** — `align`と同じ一貫性)

## chem.protein.identity_matrix

```python
from chem import protein
protein.identity_matrix(structures, chain=None)
```

リポジトリ所有者から「splitされた蛋白のidentityのマトリクスを作りたい」という要望を受けて実装。`align`が計算する配列一致度(`identity`)を、reference/構造的重ね合わせ(Kabsch fit・PDB書き出し)なしに、渡された構造集合の**全ペア**について計算する。典型的な用途は`chem.protein.split`の出力(チェーンごとに1ファイルになった蛋白質PDB群)を渡し、どのチェーンが同一蛋白質か(識別度~1.0)、無関係か(~0付近)を一覧できるようにすること。

- `structures`: PDB/CIF構造ファイルパスのリスト
- `chain`: 全構造に共通で使うチェーンIDを明示指定(省略時は各構造ごとに標準アミノ酸残基を最も多く含むチェーンを自動選択 -- `align`のreference選択と同じ基準)。`chem.protein.split`の出力は既に1ファイル1チェーンになっていることが多いので、通常は指定不要
- 戻り値: `{path_i: {path_j: identity, ...}, ...}`の辞書の辞書。使用可能なチェーンを持つ構造同士の全ペア(自分自身との組み合わせ`path_i == path_i`は`1.0`)について1エントリずつ持ち、対称(`matrix[a][b] == matrix[b][a]`)。`identity`は`align`と同じ定義(ギャップなしマッチ位置に対する一致率)の`float`、小数点3桁に丸め済み。配列比較に使えるチェーンがない構造(標準アミノ酸残基が皆無、または指定した`chain`が存在しない)は、quiet時以外は`stderr`に警告を出した上でスキップし、行・列どちらにもマトリクスへ含めない(`align`の「アラインメント不可な構造はスキップ」と同じ方針)

### アルゴリズム

1. 各`structures`について`_load_structure`で読み込み、`_select_chain(model, chain)`でチェーンを選択(`chain`引数が`None`なら自動選択、指定されていてそのチェーンIDが存在しなければ`ValueError`)し、`_chain_seq_and_ca`でシーケンス・CA原子リストを取得する。シーケンスが空(標準アミノ酸残基が1つもない)場合も`ValueError`とする。いずれかで`ValueError`が起きたパスは、quiet時以外`stderr`に警告を出してその構造をスキップし(`seqs`辞書に追加しない)、処理は継続する
2. 生き残った構造(`seqs`に登録された構造)全ての組み合わせ(`i < j`のペア)について`_matched_ca_pairs(seq_i, ca_i, seq_j, ca_j)`(`align`が使うのと全く同じヘルパー)を呼び、返り値の3番目(`identity`)だけを使う(CA座標ペアは構造的重ね合わせをしないので不要)。結果を`matrix[path_i][path_j]`と`matrix[path_j][path_i]`の両方に同じ丸め値で格納する(対称性を保証)
3. 各構造自身についても`matrix[path][path] = 1.0`をあらかじめ設定しておく(自明な自己一致なので、わざわざアラインメントを実行しない)

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

## chem.protein.residues_near

```python
from chem import protein

protein.residues_near(structure_path, ligands, radius=8.0)
```

`find_pocket`/`list_pockets`はfpocketの空洞検出に基づくが、`residues_near`はそれより軽量な、単純な距離ベースの活性部位残基抽出。リガンド自体が既に`structure_path`と同じ座標系の別ファイルとして存在するケース(典型的には`chem.protein.split`が書き出した蛋白質PDB+各リガンドSDF、必要なら`chem.protein.apply_transform`/`chem.ligand.apply_transform`で共通座標系に移した後)に向く。

- `structure_path`: PDB/CIFファイル。水・リガンド・イオンなどのHETATM残基は(`radius`以内であっても)対象外 — 常にポリマー(ATOM)残基のみを返す
- `ligands`: リガンドファイルパス1つ、またはそのリスト(`.pdb`/`.sdf`/`.mol`/`.mol2`、3D座標は`structure_path`と同じ座標系にあること)
- `radius`: 距離閾値(Å、inclusive)。デフォルト`8.0`
- 戻り値: `[{"chain", "resnum", "icode", "resname"}, ...]`(`icode`は挿入コードが無ければ`""`)、`(chain, resnum, icode)`順にソート済み — `find_pocket`の`"residues"`と同じ形

### アルゴリズム

1. `ligands`が単一パス(`str`/`os.PathLike`)なら1要素のリストにし、各パスを`find_pocket`と共通の`_load_external_ligand_coords`(RDKit、拡張子で`Chem.MolFromPDBFile`/`MolFromMolFile`/`MolFromMol2File`を使い分け)で読み込み、全リガンドの原子座標を`numpy.concatenate`で1つの配列にまとめる
2. `structure_path`を`_load_structure`で読み込み、`residue.id[0] == " "`(標準残基のhetero flag)の原子だけを対象に`Bio.PDB.NeighborSearch`を構築する(対象原子が1つも無ければ空リストを返す)
3. リガンド原子座標1つずつについて`ns.search(coord, radius)`で近傍原子を検索し、ヒットした原子の親残基を`{(chain_id, resnum, icode): residue}`の辞書に集める(重複は自動的に排除される)
4. キーでソートし、`{"chain", "resnum", "icode", "resname"}`のリストに変換して返す

## chem.protein.split

```python
from chem import protein
protein.split(structure_path, all_chains=False, remove_water=False, outdir="split")
```

構造ファイルを、ドッキング等の下流処理用に (1) リガンドフリーの蛋白質PDBと (2) 水以外の各HETATM残基インスタンス**全て**をSDF化したファイル群、に分割する。「全て」が重要な点で、結合次数テンプレートが見つからないインスタンス(共有結合した糖鎖・ペプチド様リガンドなど)であってもSDF自体は必ず書き出す — 当初の実装は`chem.ligand.load_ligand`が`ValueError`を投げたインスタンスをスキップしていたが、リポジトリ所有者から「NAGもSDFに書き出してほしい。水以外の全てのHETATMを書き出す」という明示的な訂正が入り、結合次数を復元できない場合は代わりに生の(単結合のみの)結合情報で書き出すよう変更した

- `structure_path`: PDB/CIF構造ファイルのパス
- `all_chains`: `True`にすると、蛋白質PDBをチェーンごとに分割せず、全チェーンをまとめた1ファイルにする。デフォルト`False`(=デフォルトはチェーンごとの分割 — 当初`split_chains=True`でオプトインする設計にしていたが、リポジトリ所有者から「デフォルトをチェーン単位のsplitにし、まとめたい場合だけ`all_chains=True`を指定する形にしたい」という明示的な訂正が入り、デフォルトの意味を反転した)
- `remove_water`: `True`にすると、蛋白質PDBから水(HETATM `HOH`/`WAT`/`DOD`)も取り除く。デフォルト`False`(結晶水は残す)。リポジトリ所有者から「splitの際の水分子は削除するオプション(`remove_water=False`)を追加してほしい」という明示的な要望が入り追加した — デフォルトは既存動作(水を残す)のまま、オプトインで取り除けるようにする形
- `outdir`: 出力先ディレクトリ(存在しなければ作成)。デフォルト`"split"`
- 戻り値: 以下のキーを持つ`dict`
  - `"protein"`: `all_chains=False`(デフォルト)なら`{チェーンID: パス}`の辞書、`True`なら蛋白質PDBのパス(文字列)。デフォルト(`remove_water=False`)では結晶水は残し、それ以外のHETATM残基(実際のリガンド・イオン・糖鎖修飾・結晶化添加剤など)は全て取り除く — それらは代わりに`"ligands"`側でSDF化されるため。`remove_water=True`なら水も同様に取り除く
  - `"ligands"`: `[{"path", "code", "chain", "resnum", "icode", "bond_orders_restored"}, ...]`。水以外のHETATM残基*インスタンス*ごとに**必ず**1エントリ(`chem.ligand.list_ligand_instances`と同じ粒度 — 同じコードのリガンドが複数チェーンに結合していれば別エントリになる。スキップは一切ない)。各分子の3D座標は`structure_path`そのまま。`chem.ligand.load_ligand`がPDB Chemical Component Dictionaryのテンプレートと照合して結合次数・芳香族性を復元できた場合は`bond_orders_restored=True`でその分子を、できなかった場合(共有結合した糖鎖修飾やペプチド様リガンドが、遊離型テンプレートに対して結合部位の原子が欠けている場合や、密度が不完全な残基など)は`bond_orders_restored=False`で、原子座標間の距離から単純に推定した単結合のみ・芳香族性なしの生の結合情報(quietでない限り警告を出す)を、それぞれSDFとして書き出す。`remove_water`の値に関わらず水そのものがSDF化されることは一切ない(そもそも`list_ligand_instances(..., exclude=WATER)`で列挙対象から除外している)

### アルゴリズム

1. `structure_path`をBio.PDBで読み込み、`_hetero_residues`(`chem.protein.pocket`内、標準残基・水を除く全HETATM残基を返す既存のヘルパー)で「除外すべき残基」の集合を作る。`remove_water=True`なら、読み込み済みの`model`を直接走査して`residue.id[0] == "W"`(Bio.PDBの水のhetero flag)の残基も同じ除外集合に追加する(`_hetero_residues`は内部で`structure_path`を再パースするため、既に読み込み済みの`structure`/`model`から集めた残基オブジェクトを同じ集合に混ぜても、`Bio.PDB.Entity`の`__eq__`/`__hash__`が`full_id`(モデル・チェーン・残基番号)ベースであるため問題なく機能する)
2. `Bio.PDB.PDBIO`と、上記の除外残基集合をチェックする`Select`サブクラス(`accept_residue`で`residue not in exclude_residues`、`accept_chain`で`all_chains=False`時のみ対象チェーンに絞り込み)を使い、蛋白質PDBを書き出す。`all_chains=True`なら`{outdir}/{stem}_protein.pdb`に1ファイル、デフォルト(`False`)ならチェーンごとに`{outdir}/{stem}_protein_{chain_id}.pdb`
3. `chem.ligand.list_ligand_instances(structure_path, exclude=chem.protein.WATER)`で水以外の全HETATM残基インスタンスを列挙する(`WATER`を渡すことで、`list_ligand_instances`のデフォルト`exclude=SOLVENT_AND_IONS`とは異なり、イオンや糖鎖修飾も除外せず全て対象にする — ユーザーからの明示的な要望: 例`1R1H`の`BIR`(低分子阻害剤)・`NAG`(糖鎖)・`ZN`(イオン)は全て`split`の対象とする)
4. 各インスタンスについて`chem.ligand.load_ligand(structure_path, code, chain=..., resnum=..., icode=...)`を呼ぶ
   - 成功すれば、その分子(結合次数復元済み)を使い`bond_orders_restored=True`
   - `ValueError`(テンプレート不一致)が起きた場合は、`stderr`に警告を出し(quiet時以外)、`chem.ligand.extract`の内部ヘルパー`_pick_ligand_residue`+`_residue_to_raw_mol`(`load_ligand`が内部で結合次数復元の**前**に使っているのと同じ2ステップ)を直接呼んで、結合次数復元前の生のRDKit分子(PDB由来の距離ベース単結合のみ、芳香族性なし)を取得し、`bond_orders_restored=False`とする。**例外を全体に伝播させたりインスタンスを丸ごとスキップしたりはしない** — 必ず何らかのSDFを書き出す
   - どちらの場合も`{outdir}/{stem}_ligand_{code}_{chain}{resnum}{icode}.sdf`に`Chem.SDWriter`で1分子を書き出し、`{"path", "bond_orders_restored", **instance}`を結果リストに追加する
5. `chem.ligand.extract`(このモジュールが依存する)は`chem.protein.pocket`をトップレベルでimportしているため、`splitting.py`側で`chem.ligand.extract`をトップレベルでimportすると、どちらのサブパッケージが先にimportされるかによって循環importで壊れる(片方が初期化途中のもう片方から未定義のシンボルをimportしようとする)。これを避けるため、`list_ligand_instances`/`load_ligand`/`_pick_ligand_residue`/`_residue_to_raw_mol`は`split`関数の**内部**で(呼び出し時に初めて)importする(遅延import)

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
- 新規依存: `biopython`(PyPI版あり、`pyproject.toml`の`dependencies`に追加)、`numpy`(同様)、`prody`(同様。`sequence_align`がCA配列抽出に使う。以前から`chem`環境には入っていたが`pyproject.toml`/`environment.yml`のどちらにも未宣言だったため、これを機に`dependencies`へ追加した)、`fpocket`(PyPIなし、`environment.yml`のconda-forge依存に追加。`mamba install -n chem -c conda-forge fpocket`または`mamba env update -f environment.yml`でインストール)
- `gemmi`は既存環境に(他パッケージ経由で)入っているが、今回は使わずBio.PDBのみで実装する
- `pandas`は`notebook`extrasに既に含まれる(identityマトリクスの表表示に使う)

## サンプルノートブック

`notebooks/cdk20_similar_targets.ipynb`の「CDK20's own UniProt annotation」セクションのコードセルは、`chem.protein.summary("Q8IZL9")`(または`"CDK20_HUMAN"`)を呼び、戻り値の`dict`を`pandas.DataFrame(list(props.items()), columns=["Property", "Value"])`で表(`.style.hide(axis="index")`、`Value`列は`white-space: pre-wrap`で長文を折り返し)として表示するだけにする。同ノートブックのBLASTPセクション直前、配列を取得するセルは`requests.get(".../Q8IZL9.fasta")`の代わりに`chem.protein.get_fasta("Q8IZL9")`を使う(`email`はデフォルト値のまま渡さない — `chem.blast.blastp`側の`email`引数にもデフォルト値が実装されたため、以前ノートブック内に直接記述していた`EBI_EMAIL`変数はリポジトリ所有者の指示で削除した)。

`notebooks/sequence_alignment.ipynb`(新規)は、特定の1蛋白専用ではなく`sequence_align`を汎用的に使うためのテンプレートノートブックとして作る。全セルはStep 1の設定セル(`UNIPROT_ACCESSION`・`CANONICAL_FEATURE_TYPE`/`CANONICAL_FEATURE_DESCRIPTION`・`MARKER_FEATURE_TYPES`・`PDB_IDS`・`DATA_DIR`)を書き換えるだけで別の蛋白・別の構造セットに差し替えられるようにする(当初はHIV-1プロテアーゼ(`P03366`のGag-Polポリプロテインから`Chain`feature`"Protease"`で切り出し)、続けてTEM-1 β-ラクタマーゼ(`P62593`、`Chain`feature`"Beta-lactamase TEM"`でシグナルペプチドを除く)、BRD4のBromo1ドメイン(`O60885`、`Domain`feature`"Bromo 1"`)と対象を差し替えながら育てたが、いずれも構造が良く折り畳まれており電子密度の欠損(ギャップ)がほとんど出ない題材だったため、リポジトリ所有者の判断で最終的にβ2アドレナリン受容体(`ADRB2_HUMAN`、`P07550`、`CANONICAL_FEATURE_TYPE=None`でfull-length使用、`MARKER_FEATURE_TYPES=("Binding site",)`でオルソステリックポケット残基をマーカーに)に差し替えた — T4リゾチーム融合・ナノボディ・ヘテロ三量体Gタンパク質複合体など構築物の多様性が高く、チェーン自動選択やギャップ処理を実地で検証できる題材のため)。最終的なPDB_IDSは18構造、インバースアゴニスト(`2RH1`のcarazolol等)・中性アンタゴニスト(`3NYA`のalprenolol)・共有結合性アゴニスト(`3PDS`)・フル/部分アゴニスト(`4LDO`のアドレナリン、`7DHI`のsalbutamol等)・正/負のアロステリックモジュレーター(`6N48`/`6OBA`)・細胞内アロステリック部位アンタゴニスト(`5X7D`)・Gs蛋白複合体(`3SN6`、`7DHR`等cryo-EM)を薬理学的に幅広くカバーする。

セル構成: (1) タイトル+概要md(「アラインメント自体は`sequence_align`に任せ、このノートブックは表示に専念する」設計方針を明記)、(2) Step 1設定md+コード、(3) Step 2 PDBダウンロードmd+コード(`chem.rcsb.download_structures(PDB_IDS, outdir=DATA_DIR, filetype="pdb")`後、`structure_paths = [os.path.join(DATA_DIR, f"{pdb_id}.pdb") for pdb_id in PDB_IDS]`でパスリストを作る)、(4) Step 3 `sequence_align()`呼び出しmd+コード: `result = chem.protein.sequence_align(...)`を実行し、`result["sequences"]`をファイルパスからPDB ID(`os.path.splitext(os.path.basename(p))[0]`)にre-keyした`observed`辞書を作り、蛋白名・organism・マーカー位置・canonical配列を`print`した後、`sequences_by_label = {"Canonical": canonical_seq, **observed}`として**canonicalも含めた**全蛋白の総当たりペアワイズidentityマトリクス(gap位置は分母分子どちらからも除外、`chem.protein.structural_align.identity_matrix`と同じ定義)を`pandas.DataFrame`で作り、対角(自明な自己一致1.0)は`np.nan`にして`na_rep="-"`で`"-"`表示、非対角でidentity`1.0`(=完全に同一の構築物系統)のセルだけ`.style.map(...)`で黄色背景にする(`identity_matrix.to_numpy(dtype=float, copy=True)`を使う点に注意 — pandasのcopy-on-write下では`.values`への直接書き込みが読み取り専用エラーになるため、素の`.copy()`では不十分)、(5) Step 4 色付きHTML表示md+コード: `render_alignment_html(observed, canonical_seq, marker_pos, width=100, reference=None)`。60→100文字幅への変更、ブロック先頭に10残基おきの位置ルーラー(`_ruler_lines`ヘルパー、右詰め数字行+`|`行)、相違点は赤太字、canonicalに無い挿入(gapped_canonicalの`-`列)は落として詰めるためHTML上には出ない、canonical側に無い(未observed)位置は灰色イタリック、マーカー位置は太字+下線(`<b><u>`)。`reference`引数で基準をcanonical以外の構造(例: `2RH1`)に切り替えられ、その場合基準構造自身が観測していない位置は比較不能として赤字にせずそのまま表示する、(6) Step 5 2RH1(最初の高分解能構造)を基準にした表示md+`render_alignment_html(observed, canonical_seq, marker_pos, reference="2RH1")`のコール。

このノートブック用のアラインメントロジック(タイブレーク、`mismatch_score=-3`)は`sequence_align`本体に実装されているため、ノートブック側は`render_alignment_html`・`pairwise_identity`など表示専用のヘルパーのみを持つ。実行して確認できた具体的な生物学的シグナル: `N187E`(N-グリコシル化部位除去、結晶化構築物の定番)が18構造全てに共通、`M96T`/`M98T`が新しめの構築物系統(cryo-EM Gs複合体含む)にのみ共通、`H93C`が共有結合性アゴニスト構造`3PDS`だけに出現(タイトル通りアゴニストをシステインに共有結合させるための変異と整合)、`6KR8`(フルアゴニスト結合状態)は`A59C`/`T136C`/`N148C`等の新規Cys導入と`C77V`/`C265A`/`C327S`等の天然Cys除去を組み合わせた「システイン再設計」構築物であることが変異点から読み取れる。identityマトリクスでは`3NY8`/`3NYA`/`6PS3`/`6PS5`同士、`7DHR`/`7BZ2`/`7DHI`同士、`5JQH`/`4LDO`/`6N48`同士がそれぞれペアワイズ`1.000`(完全一致)になり、構築物の世代・系統がそのまま可視化される。

`notebooks/alphafold_pocket_thrb_human.ipynb`のRCSBダウンロードセルの直後に、ダウンロードした複数のトロンビン構造をAlphaFold予測構造を参照にアラインする独立セクション(タイトルmd + コード。コードは`align()`実行と`align_df`(rmsd/identity列)の表示のみ)を置く。重ね書きpy3Dmolビューア(トグルボタンで表示構造を選ぶウィジェット)はさらに別の独立セクション(独自のH2タイトルmd + コード)として、アラインメントセクションの直後に続ける -- 1セルに両方を詰め込まない。fpocketのサンプルは別途、リガンド入り構造(例: トリプシン+ベンザミジン `3PTB`)で`find_pocket`を実行し、選ばれたポケットの残基をハイライト表示するセルを追加する。`list_pockets`のサンプルは、AlphaFold予測構造セクション(リガンドが存在しない)に以下3セルを追加する: (1) `pandas.DataFrame`で「pocket_id / score / druggability_score / volume / n_residues」の表として全候補ポケットを表示、(2) `druggability_score >= 0.2`の候補ポケットを、半透明cartoonの上にポケットごとに異なる色の構成残基stickでハイライトして可視化(色とpocket_id/druggability_scoreの凡例付き)、(3) 同じ候補ポケットを、構成残基のstickの代わりに`spheres`フィールド(fpocketのアルファ球)を`py3Dmol.addSphere`で球ごとに描画し、空洞を充填された体積として可視化(こちらもポケットごとに色分け・凡例付き)。既存セルは書き換えず、新規セルとして追記する。

`notebooks/cdk9_cdk2_split.ipynb`(新規)は、同一の2-amino-4-heteroaryl-pyrimidine系
阻害剤シリーズをCDK9-サイクリンT複合体4構造(`4BCF`/`4BCH`/`4BCI`/`4BCJ`、キナーゼ+サイクリンの
2チェーン)とCDK2-サイクリンA複合体5構造(`4BCK`/`4BCM`/`4BCN`/`4BCO`/`4BCQ`、非対称単位に
2コピー含む4チェーン)にそれぞれ結合させた計9構造を`chem.rcsb.download_structures`
(PDB idのリストを渡す形)でダウンロードし、`chem.protein.split`を9構造それぞれに適用する
(当初はネプリライシン3構造`1R1H`/`1R1I`/`1R1J`の単一チェーン例だったが、リポジトリ所有者が
「もっと良い例」としてこの9構造を提示したため、そちらの`notebooks/neprilysin_split.ipynb`は
削除しこのノートブックに差し替えた。複数チェーン構成の違い(2チェーン vs 4チェーン)を
横断的に確認できる点、および活性化ループのリン酸化トレオニン`TPO`という、糖鎖修飾`NAG`とは
別種の「結合次数テンプレートと一致しない」ケースが一貫して現れる点で、より充実した題材になっている)。
セル構成: (1) タイトル+概要md、(2) 9構造ダウンロード、(3) `split`の説明md(`remove_water=True`
を使う旨も明記)、(4) 9構造をループして`split(..., remove_water=True, outdir="cdk9_cdk2_split")`
実行(デフォルトの`all_chains=False`のまま、`"protein"`はCDK9セットが`{"A", "B"}`、
CDK2セットが`{"A", "B", "C", "D"}`の辞書になる)、(5) チェーン構成確認md、(6) `entry_id`/
`n_chains`/`chain_ids`の`pandas.DataFrame`でチェーン数の違いを表示、(7) 分割結果の説明md、
(8) 各構造の`split`結果(`"ligands"`)をそのまま`entry_id`列付きで結合した`pandas.DataFrame`で
一覧表示、(9) コードごとに`bond_orders_restored`が常に同じ値になること(`TPO`は常に`False`、
それ以外は常に`True`)を`groupby`で集計表示、(10) 網羅性確認md、(11) `chem.ligand.list_ligand_instances`
(`exclude=WATER`)で構造ファイル自身から数えた水以外のHETATM残基インスタンス数と、`split`が
実際に書き出したSDF数が一致すること(`assert`)を9構造それぞれで確認、(12) リガンドフリー
蛋白質PDBの確認md、(13) `"protein"`辞書の全チェーンファイル(計28ファイル)に対して再度
`list_ligand_instances`を実行し水以外のHETATM残基が0件であることを`assert`、
(14) `remove_water=True`確認md、(15) `Bio.PDB`で28ファイル全てを直接パースし`res.id[0] == "W"`
(水のhetero flag)の残基が1件も無いことを`assert`(`list_ligand_instances`は水自体を対象外に
しているため、この確認だけは別途必要)、(16) `TPO`/阻害剤比較md、(17) `4BCF`の`TPO`
(`bond_orders_restored=False`)と`T6Q`(`True`)のSDFを読み込み直し、ボンドタイプ集合
(前者は`{"SINGLE"}`のみ、後者は`AROMATIC`を含む)を比較、(18) 阻害剤SAR比較md、
(19) `T6Q`/`T7Z`/`T3E`/`T9N`/`TJF`5化合物についてコードごとに1回だけSDFを読み込み直し、
`bound_to_entries`(その化合物が結合していた構造id、重複除去してソート)・原子数・
`chem.ligand.molecular_weight`/`qed`・SMILESを`pandas.DataFrame`で比較(`T6Q`/`T7Z`/`T3E`/`T9N`は
CDK9・CDK2の両方に、`TJF`はCDK2のみに結合していることが`bound_to_entries`列に現れる)、
(20) 可視化md、(21) 4チェーン構造`4BCK`の蛋白質PDB4ファイルをチェーンごとに色分けした
cartoonとして、両キナーゼコピー(チェーンA/C)に結合した阻害剤`T3E`のSDF(stick)を
py3Dmolで重ねて表示、(22) `identity_matrix`の説明md(リポジトリ所有者の「splitされた蛋白の
identityのマトリクスを作りたい」という要望を受けて追加。キナーゼ同士・サイクリン同士は
~1.0、キナーゼ vs サイクリンは無関係なフォールドなので~0.15、CDK9 vs CDK2のキナーゼドメイン
同士は同ファミリーのパラログなので~0.4前後になる、という期待値の説明を含む)、
(23) 28ファイル全ての`"protein"`パスと`{entry_id}_{chain_id}`ラベルを集め、
`chem.protein.identity_matrix`を1回呼んで28×28の`numpy`配列に変換、(24) `matplotlib`の
`imshow`でヒートマップとして可視化(カラーバー付き、軸ラベルに`{entry_id}_{chain_id}`)。

同ノートブックは続けて`chem.protein.compute_transform`/`apply_transform`/`chem.ligand.apply_transform`/
`chem.protein.residues_near`のデモを追加する: (25) 説明md(「片方のチェーンだけアラインし、
同じ変換を他方に使い回す」旨。同一PDBエントリ内のチェーンは元々同じ結晶の非対称単位に
収まっているため、キナーゼ側1本の変換をそのエントリの全チェーン・全リガンドSDFへそのまま
適用してよい、という理屈も明記)、(26) 各エントリのキナーゼ(chain A)を基準構造`4BCF`へ
`compute_transform`で変換計算(`4BCF`自身は単位回転・ゼロ並進を手動セット)、
(27) 変換適用md、(28) 各エントリの全蛋白質チェーンと全リガンドSDFに、キナーゼ1本分の
変換をそのまま`protein.apply_transform`/`ligand.apply_transform`で適用し`cdk9_cdk2_aligned/`
に書き出し(非対称単位内2コピー目のチェーンC/Dも含め、個別にアラインし直さない点がポイント)、
(29) 検算md、(30) 変換後のキナーゼchain Aを改めてreferenceに対し`compute_transform`し直し、
回転行列がほぼ単位行列・並進がほぼゼロ・rmsdが元の値とほぼ一致することを確認、
(31) 可視化md、(32) 9構造全ての全チェーン+主要阻害剤(`TPO`除く)を共通座標系に重ねてpy3Dmol
表示(CDK9セットとCDK2セットで色分け)、(33) アクティブサイト判定md、(34) 同一chainに
結合した`TPO`以外のリガンド候補のうち最大重原子数のものを本命とし、本命から半径8Å以内の
候補だけを残す2段階アルゴリズム(`active_site_ligands`)で結晶化添加剤(`4BCK`の`SGM`など)
を除外、(35) 可視化md、(36) キナーゼ(chain A)のみ+その活性部位リガンドのみを重ねて表示、
(37) `chem.protein.residues_near`の説明md、(38)-(39) `residues_near`(半径8Å)でアクティブ
サイト残基を抽出しラベル付け、`ipywidgets.SelectMultiple`で表示する構造(entry_id)を選べる
インタラクティブなpy3Dmolビュー(選んだ構造すべてを同じ座標系に重ねて比較できる)。

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- テストは`tests/test_protein_sequence_align.py`に、ネットワーク不要な範囲で追加する: `_marker_positions`が`marker_feature_types=None`/空タプルで空集合を返すこと・指定typeでfeature位置をcanonical相対座標に正しく変換すること、`_fetch_canonical`(`requests.get`を`monkeypatch`)がfeature指定時にスライスすること・`canonical_feature_type=None`でfull-lengthをそのまま使うこと、`_align_to_canonical`が長い異種挿入(融合パートナー相当)を1つの連続したギャップに畳み込むこと(canonicalの置換対象範囲が十分長い場合)・同点最適解群からマッチブロック数最小のものを選ぶこと(canonical末尾が繰り返し残基のケース)、置換対象範囲が短い(2残基)場合はなおスミアリングし得るという既知の限界そのものを固定するテスト、`_load_ca_sequences`が全チェーンのCA配列を返すこと、`_best_matching_chain_sequence`が(`chain=None`なら)最もアラインメントスコアの高いチェーンを選ぶこと・`chain`明示指定時はそれを尊重すること、`sequence_align`をエンドツーエンドで(`requests.get`を`monkeypatch`、合成PDBテキストを`tmp_path`に書き出して)実行し、`protein_name`/`organism`/`canonical_seq`/`feature_start`/`feature_end`/`marker_positions`/`raw_sequences`/`sequences`が期待通りの内容になること(1残基欠損・1残基点変異を含む合成構造で確認)。実データでの動作確認は`notebooks/sequence_alignment.ipynb`のβ2アドレナリン受容体18構造で行い、`chem.protein.__init__`が`from .sequence_align import sequence_align`で`sequence_align`という関数名をサブモジュール名と同名で再エクスポートしているため(既存の「実装ファイル名を公開関数名と同じにしない」原則の裏返しの注意点として)、`import chem.protein.sequence_align as x`のような`from ... import ... as`形式でサブモジュール自体にアクセスしようとすると`x`は関数オブジェクトになってしまう(`chem.protein`パッケージの`sequence_align`属性が関数で上書きされているため)。テストでモジュール内のプライベートヘルパー(`_fetch_canonical`等)に直接アクセスする場合は`importlib.import_module("chem.protein.sequence_align")`を使うこと
- テストは`tests/test_protein_split.py`に、`split`のロジック(合成PDBテキストに対する、
  デフォルト(`all_chains=False`)でのチェーンごとの蛋白質PDB書き出し・ファイル名への
  チェーンID埋め込みと水以外のHETATM残基除去、`outdir`未存在時の自動作成、
  `all_chains=True`での単一蛋白質PDBへの統合書き出し、デフォルト(`remove_water=False`)では
  蛋白質PDBに水が残り`remove_water=True`では水も取り除かれること(いずれの場合も水自体が
  SDF化されることはない)、リガンドSDFの書き出しと
  `AssignBondOrdersFromTemplate`によるベンゼン環の芳香族性復元(`bond_orders_restored=True`)、
  テンプレートの取得に失敗するインスタンス(`requests.get`を`monkeypatch`で偽装し、該当コードに
  ディスクリプタを一切返さない)がスキップされずに`_pick_ligand_residue`+`_residue_to_raw_mol`
  による生の単結合SDFとして書き出され(`bond_orders_restored=False`)、かつ蛋白質PDB側からは
  正しく取り除かれること)を、`chem.ligand.extract`の`requests.get`を`monkeypatch`することで
  ネットワーク不要なオフラインテストとして追加する。実データでの動作確認は
  `4BCF`/`4BCH`/`4BCI`/`4BCJ`(CDK9-サイクリンT複合体、2チェーン)・
  `4BCK`/`4BCM`/`4BCN`/`4BCO`/`4BCQ`(CDK2-サイクリンA複合体、4チェーン)の計9構造で
  `split`を実行し、活性化ループのリン酸化トレオニン`TPO`(主鎖内の修飾残基)・低分子阻害剤
  (`T6Q`/`T7Z`/`T3E`/`T9N`/`TJF`)・結晶化添加剤(`GOL`/`SGM`/`SO4`)の全インスタンスが
  (構造ファイル自身から`list_ligand_instances`で数えた件数と一致する形で)1つも欠けずSDF化され、
  `TPO`は9構造・全14インスタンスにわたって一貫して`bond_orders_restored=False`(生の単結合)、
  それ以外は全て`True`(結合次数復元済み)になること、2/4チェーンいずれの構成でも蛋白質PDBが
  正しくチェーンごとに分割され水以外のHETATM残基が残らないこと、`T6Q`/`T7Z`/`T3E`/`T9N`が
  CDK9・CDK2の両方に、`TJF`がCDK2のみに結合しているという既知の実験事実が
  `bound_to_entries`集計に正しく現れること、`identity_matrix`を全28チェーンPDBに対して
  実行すると同一蛋白質の複数コピー同士(例: 4構造にわたるCDK9キナーゼチェーン)が~0.997-1.0、
  キナーゼ vs サイクリン(無関係なフォールド)が~0.15、CDK9 vs CDK2のキナーゼドメイン同士
  (同ファミリーのパラログ)が~0.4前後になり、ヒートマップ上で4つの高一致度ブロックとして
  はっきり分かれることを`notebooks/cdk9_cdk2_split.ipynb`の実行で確認した
- テストは`tests/test_protein_align.py`・`tests/test_protein_pocket.py`にネットワーク・外部バイナリ(fpocket)不要な範囲(`identity_matrix`が単一構造で自己一致`1.0`のみの辞書を返すこと、同一配列同士のペアが`1.0`・異なる配列同士のペアが対称(`matrix[a][b] == matrix[b][a]`)かつ`1.0`未満になること、標準アミノ酸残基を持たない構造(HETATMのみ)がスキップされ行・列どちらにも現れないこと、`chain`引数が実際に尊重され(サイズの大きい方のチェーンを自動選択するのではなく)指定したチェーン同士が比較されること、シーケンスマッチングロジック、`_matched_ca_pairs`が返す`identity`について同一配列で`1.0`・ギャップを含む場合はギャップ位置を分母/分子どちらからも除外・ギャップなしミスマッチを含む場合は分母に数えて分子には数えないこと、`align()`の戻り値が`{"rmsd":..., "identity":...}`の形でreferenceは`{"rmsd": 0.0, "identity": 1.0}`になること、`_best_matching_chain_alignment`が「マッチ位置数は多いが一致度が低い大きなチェーン」より「短くても一致度が高いチェーン」を選ぶこと(`3B9F`相当の合成データで再現)、`align()`をエンドツーエンドで実行してもサイズ最大の無関係なチェーンではなく正しいチェーンが選ばれ`identity`が高くなること、リガンド自動検出・HETコード判定・ファイル判定の分岐、ポケット選択の距離計算、`_info.txt`パーサ、fpocket未インストール時のエラーメッセージ、`_pocket_atm_paths`/`_pocket_result`の単体動作(`pocket{N}_vert.pqr`が存在する/しない両方のケース)、`_parse_pocket_spheres`が合成PQRテキストから`x`/`y`/`z`/`radius`を正しく取り出すこと、`list_pockets`を`_run_fpocket`を`monkeypatch`で偽の出力ディレクトリに差し替えて実行し`druggability_thres=None`なら全ポケットが`druggability_score`降順(値なしは最後)で返ること・デフォルト(`0.1`)では値なし/閾値未満のポケットが除外されること・`druggability_thres`を明示指定すればその閾値で絞り込まれること、`WATER`が`SOLVENT_AND_IONS`の真部分集合であること)のみ追加する。実際のBio.PDB構造アラインメントとfpocket実行は、簡易的な合成PDBテキスト(固定カラム位置で手書きしたATOM/HETATMレコード)を使ったオフラインテストと、実データでの手動実行確認(3PTB+BENでPocket 1が選ばれAsp189・Ser195が残基リストに含まれること、AlphaFold予測構造(リガンド無し)で`list_pockets`が複数候補を返し、そのうち`druggability_score >= 0.2`の3件を`spheres`経由でHTMLに書き出しブラウザで表示確認したところ、cartoon上にポケットごとに色分けされた充填体積として描画されることを確認済み)の組み合わせで検証する
