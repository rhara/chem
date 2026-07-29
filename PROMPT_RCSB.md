# chem.rcsb.download_structures の再現プロンプト

`chem`リポジトリの`chem`パッケージ内に、RCSB PDB構造ダウンロード用のサブパッケージ`rcsb`(`chem.rcsb`)を追加し、ターゲット別のPDB構造ファイルをダウンロードする関数を実装するための指示。`chem.chembl`と対になる機能で、ID解決ロジックを共有する。

## 要件

`src/chem/rcsb/fetch.py` に関数`download_structures`を実装し、`src/chem/rcsb/__init__.py`で`from .fetch import download_structures`として再エクスポートする。以下の形で呼び出せること:

```python
from chem import rcsb
rcsb.download_structures(id, resolution_thres=2.0, outdir="data", filetype="cif")
```

- `id`: ターゲット蛋白の識別子、またはPDBエントリIDのリスト
  - ターゲット識別子として渡す場合、`chem.chembl.download_activities`と全く同じ3種類を受け付ける
    - ChEMBLターゲットID (例: `CHEMBL204`)
    - UniProt accession (例: `P00734`)
    - UniProt Entry名/mnemonic (例: `THRB_HUMAN`, `BRAF_HUMAN`)
  - list/tuple/setとして渡す場合は、ダウンロードしたいPDBエントリID(4文字、例: `["6LU7", "7BQY"]`)そのものとみなし、ターゲット解決・RCSB Search APIによる検索を一切行わずそれらのエントリを直接ダウンロード対象にする(大文字小文字は区別せず内部で大文字に正規化する)
- `resolution_thres`: 省略可(デフォルト`None`)。指定時は分解能(Å)がこの値以下(inclusive)の構造のみを対象にする。分解能情報を持たない構造(NMRなど)は、`resolution_thres`が指定されているときは常に除外し、`None`のときは(分解能の有無によらず)全て含める。`id`にPDBエントリIDのリストを渡した場合も同様にフィルタリングされる
- `outdir`: 出力先ディレクトリ。存在しなければ作成する
- `filetype`: `"cif"`(デフォルト)、`"pdb"`、`"both"`のいずれか。ダウンロードするファイル形式を指定する
- 戻り値: 少なくとも1ファイルが揃った(新規ダウンロードまたは既存ファイル)PDBエントリ数

## ID解決ロジックの共通化

`chem.chembl.fetch`にあったID解決関連の関数を`src/chem/ids.py`に切り出し、`chem.chembl.fetch`と`chem.rcsb.fetch`の両方から`from ..ids import ...`で利用する:

- `resolve_uniprot_accession(id_)`: UniProt accession/entry名 → UniProt accession(既存の`_resolve_uniprot_accession`をそのまま移動)
- `resolve_target_chembl_id(id_)`: ChEMBLターゲットID/UniProt accession/entry名 → ChEMBLターゲットID(既存の`_resolve_target_chembl_id`をそのまま移動。`chem.chembl.download_activities`が使用)
- `resolve_uniprot_accession_any(id_)`: ChEMBLターゲットID/UniProt accession/entry名 → UniProt accession(新規追加。`chem.rcsb.download_structures`が使用)。ChEMBLターゲットIDが渡された場合は、ChEMBL API `target/{id}.json`を叩いて`target_components`から`component_type == "PROTEIN"`の`accession`を逆引きする。それ以外は`resolve_uniprot_accession`にそのまま委譲する

`chem.chembl.fetch`の`_CHEMBL_ID_RE`/`_UNIPROT_ACCESSION_RE`/`_resolve_uniprot_accession`/`_resolve_target_chembl_id`は削除し、`ids.py`のものを使う(既存のテスト`test_resolve_chembl_id_passthrough`も`tests/test_ids.py`に移動し、`chem.ids`を直接テストする形にする)。

## データ取得方法

1. `id`がlist/tuple/setの場合は、各要素を大文字に正規化し4文字(先頭が数字、残り3文字が英数字)のPDBエントリIDとして妥当性検証する(空リストや形式不正な要素があれば`ValueError`)。それ以外の場合は`resolve_uniprot_accession_any(id)`でUniProt accessionを得た上で、**RCSB Search API** (`https://search.rcsb.org/rcsbsearch/v2/query`、POST)でそのaccessionに紐づく全PDBエントリIDを検索する。クエリはAND条件の2ノード:
   - `rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession` の`exact_match`がaccession
   - `rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_name` の`exact_match`が`"UniProt"`
   - `return_type: "entry"`、`request_options.paginate`でページネーション(1回1000件、`total_count`に達するまで`start`を進める)
   - ヒットが0件の場合(HTTP 204または空の`result_set`)は`ValueError`
2. 得られたエントリID群(list指定の場合はそれ自身)について、**RCSB GraphQL API** (`https://data.rcsb.org/graphql`、POST)で分解能を一括取得する。200件程度ずつバッチ化し、クエリは`entries(entry_ids: $ids) { rcsb_id rcsb_entry_info { resolution_combined } }`。`resolution_combined`が空/nullなら分解能なし(`None`)として扱う(NMR構造など)。バッチ処理は`tqdm`で進捗表示(`CHEM_QUIETNESS`が非quietのときのみ)
3. `resolution_thres`でフィルタリングした上で、各エントリを`https://files.rcsb.org/download/{entry_id}.{ext}`からダウンロードする(`ext`は`filetype`に応じて`cif`/`pdb`、`"both"`なら両方)。ダウンロードも`tqdm`で進捗表示する
4. `outdir`に`{entry_id}.{ext}`として保存する。**保存先に同名ファイルが既に存在する場合はネットワークアクセスせずスキップする**(既存ファイルの中身が最新かどうかの検証は行わない。同じentry_id+拡張子の内容はRCSB側で不変なため)
5. 一部の大きな構造(cryo-EMアセンブリなど)は legacy `.pdb`形式が提供されておらず404になることがある。その場合はその形式だけスキップし(quiet時以外は`stderr`に警告)、処理は継続する

## 呼び出しログと進捗出力

`chem.chembl.download_activities`と同様に、`download_structures`に`@logged`デコレータ(`chem.verbosity.logged`)を適用する。終了時の`"wrote N structures to {outdir}"`メッセージも`is_quiet()`のとき出力しない。

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- 新規の外部依存追加は不要(既存の`requests`、`tqdm`のみ使用)
- `src/chem/rcsb/`として`chem`パッケージのサブパッケージに配置(`pyproject.toml`の`[tool.setuptools.packages.find]`が`chem.rcsb`を自動検出)

## サンプルノートブック

`notebooks/example_proteins.ipynb`に、`THRB_HUMAN`(トロンビン)を対象としたサンプルセルを追加する:
1. `rcsb.download_structures("THRB_HUMAN", resolution_thres=1.2, outdir="data")`を実行(`from chem import rcsb`。デモを軽くするため分解能閾値は厳しめにして数件程度に絞る)
2. ダウンロードされた`.cif`ファイル一覧を表示
3. `py3Dmol`で最初の構造をカートゥーン表示

生成される`data/`ディレクトリは`.gitignore`に追加し、コミットしない。

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- テストは`tests/test_rcsb_fetch.py`にネットワーク不要な範囲(`_select_entries`の分解能フィルタリングロジック、`_download_one`の既存ファイルスキップ、不正な`filetype`のバリデーション、`_validate_pdb_ids`の正規化・バリデーション、PDBエントリIDリストを渡したときにターゲット解決・検索が呼ばれないことの確認)のみ追加する。ネットワークを伴う関数(UniProt/ChEMBL解決、RCSB検索・ダウンロード)はユニットテスト対象外とし、手動での実行確認(実際にAPIを叩く)で動作検証する
