# chem.chembl.download_activities の再現プロンプト

`chem`リポジトリの`chem`パッケージ内に、ChEMBLデータ取得用のサブパッケージ`chembl`(`chem.chembl`)を追加し、ChEMBLのターゲット別バイオアクティビティデータをダウンロードする関数を実装するための指示。

## 要件

`src/chem/chembl/fetch.py` に関数`download_activities`を実装し、`src/chem/chembl/__init__.py`で`from .fetch import download_activities`として再エクスポートする。以下の形で呼び出せること:

```python
from chem import chembl
chembl.download_activities(id, mw=[250, 650], normalize_smiles=True, output="test.tsv")
```

- `id`: ターゲット蛋白の識別子。以下3種類すべてを受け付ける
  - ChEMBLターゲットID (例: `CHEMBL204`)
  - UniProt accession (例: `P00734`)
  - UniProt Entry名/mnemonic (例: `THRB_HUMAN`, `BRAF_HUMAN`)
- `mw`: 省略可。`[下限, 上限]` の分子量範囲(両端inclusive)。指定時はこの範囲外の化合物を除外する
- `normalize_smiles`: 省略可、デフォルト`False`。`True`の場合、化合物をChEMBL Structure Pipelineで標準化・脱塩し、さらに重複化合物(正規化後SMILESが同一)を1行に集約する(下記参照)
- `output`: 出力ファイルパス。拡張子`.csv`ならカンマ区切り、それ以外(`.tsv`など)はタブ区切りで書き出す
- 戻り値: 書き出した行数

## データ取得方法

- **ChEMBL Web REST API** (`https://www.ebi.ac.uk/chembl/api/data`) を`requests`で直接叩く。`chembl_webresource_client`パッケージは使わない
- ID解決:
  - `^CHEMBL\d+$` にマッチすればそのままターゲットChEMBL IDとして使う
  - それ以外はUniProt REST API (`https://rest.uniprot.org/uniprotkb/search`) でaccessionに解決してから、ChEMBLの`target.json?target_components__accession=<accession>`でターゲットChEMBL IDを引く
  - UniProtクエリでは、入力が正式なUniProtKB accessionの形式(正規表現 `^([A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[OPQ][0-9][A-Z0-9]{3}[0-9])$`)にマッチする場合のみ`accession:`フィールドを使う。マッチしない場合(Entry名)は`id:`フィールドを使う。`accession:`フィールドは形式が不正だと400エラーを返すため、フィールドの出し分けが必須
  - 複数ターゲットがヒットした場合は`target_type == "SINGLE PROTEIN"`を優先する
- Activity取得: `activity.json?target_chembl_id=<id>&pchembl_value__isnull=false&limit=1000` をページネーション(`page_meta.next`)で全件取得する(pChEMBL値を持つレコードのみをAPI側で絞り込む)。取得件数が多い(数千件規模)ため、`tqdm`で進捗を表示する(`CHEM_QUIETNESS`が非quietのときのみ、下記参照)
- pChEMBL値を持たないレコードは除外する(`pchembl_value is None`のものはスキップ)

## 化合物処理

- 各activityレコードの`canonical_smiles`をRDKitでパース
- `normalize_smiles=True`の場合: [ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline)(`pip`/`conda-forge`パッケージ`chembl_structure_pipeline`)の`standardizer.standardize_mol()` → `standardizer.get_parent_mol()`で標準化・脱塩し、`Chem.MolToSmiles`で正規化する
- 分子量は`rdkit.Chem.Descriptors.MolWt`で計算し、`mw`が指定されていれば範囲外の行をスキップする

## 出力列

`normalize_smiles=False`(1activity=1行):
```
molecule_chembl_id, assay_chembl_id, target_chembl_id, document_chembl_id, pchembl_value, smiles, mw
```

`normalize_smiles=True`(正規化後SMILESで重複化合物を集約、1化合物=1行):
```
parent_chembl_id, target_chembl_id, smiles, mw, n, pchembl_mean, pchembl_median, pchembl_std
```
- 集計キーは正規化後のSMILES(異なる`molecule_chembl_id`でも脱塩後に同一構造になれば同一化合物として扱う)
- `parent_chembl_id`はグループ内で最初に出現した`molecule_chembl_id`を代表値として採用
- `pchembl_mean`/`pchembl_median`/`pchembl_std`は`statistics.fmean`/`median`/`pstdev`(母標準偏差、n=1のときは0)で計算し、小数点3桁に丸める

## 呼び出しログと進捗出力 (CHEM_QUIETNESS)

`src/chem/verbosity.py` に共通デコレータを実装し、`chem.chembl.download_activities`(実体は`chem/chembl/fetch.py`)に適用する:

- `is_quiet()`: 環境変数`CHEM_QUIETNESS`が未設定なら`False`(非quiet)。設定されていて値が`"0"`/`"N"`/`"FALSE"`(大文字小文字不問)のいずれでもなければ`True`(quiet)
- `@logged`デコレータ: `is_quiet()`が`False`のとき、呼び出された関数名と実引数(デフォルト値も含めてbindしたもの)を`関数名(引数名=値, ...)`の形式で標準エラー出力に書き出す
- `download_activities`のtqdm進捗バーは`disable=is_quiet()`を渡し、quiet時は表示しない
- `download_activities`終了時の`"wrote N rows to ..."`メッセージも`is_quiet()`のとき出力しない
- 今後`chem`/`chembl`に追加する公開関数にも同じ`@logged`を使い回す想定

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- `pyproject.toml`の`dependencies`: `rdkit`, `py3dmol`, `tqdm`, `requests`, `chembl_structure_pipeline`
- `environment.yml`(conda-forge)にも`chembl_structure_pipeline`を追加
- `src/chem/chembl/`として`chem`パッケージのサブパッケージに配置(`pyproject.toml`の`[tool.setuptools.packages.find]`が`chem.chembl`を自動検出)
- `chem/chembl/fetch.py`は同じ`chem`パッケージ内の`chem.verbosity`を相対import(`from ..verbosity import ...`)で使う
- `chem/chembl/__init__.py`は`from .fetch import download_activities`で公開関数を再エクスポートし、利用側は`chem.chembl.fetch`ではなく`chem.chembl`から直接呼び出す

## サンプルノートブック

`notebooks/00_getting_started.ipynb`に、`BRAF_HUMAN`を対象としたサンプルセルを追加する:
1. `chembl.download_activities("BRAF_HUMAN", mw=[250, 650], normalize_smiles=True, output="braf_activities.tsv")`を実行(`from chem import chembl`)
2. `pandas`で読み込み、`pchembl_mean`降順で上位を表示
3. 最も活性の高い化合物をRDKitの`MolToImage`で描画

`pandas`は`[notebook]` extraに追加する(コア依存ではない)。生成される`.tsv`/`.csv`は`.gitignore`に追加し、コミットしない。

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない(ただしChEMBL Structure Pipelineのようなpublicな公式OSSツールの利用は問題ない)
- テストは`tests/test_fetch.py`・`tests/test_verbosity.py`にネットワーク不要な範囲(ID直接一致のパススルー、正規化ロジック、集計ロジック、ログ抑制ロジック)のみ追加する。ネットワークを伴う関数(UniProt解決・ChEMBL取得)はユニットテスト対象外とし、手動での実行確認(実際にAPIを叩く)で動作検証する
