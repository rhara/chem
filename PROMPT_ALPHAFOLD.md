# chem.alphafold.download_structures の再現プロンプト

`chem`リポジトリの`chem`パッケージ内に、AlphaFold DB構造ダウンロード用のサブパッケージ`alphafold`(`chem.alphafold`)を追加し、ターゲット別のAlphaFold予測構造ファイルをダウンロードする関数を実装するための指示。`chem.chembl`・`chem.rcsb`と対になる機能で、ID解決ロジック(`chem.ids`)を共有する。

## 要件

`src/chem/alphafold/fetch.py` に関数`download_structures`を実装し、`src/chem/alphafold/__init__.py`で`from .fetch import download_structures`として再エクスポートする。以下の形で呼び出せること:

```python
from chem import alphafold
alphafold.download_structures(id, plddt_thres=70.0, outdir="data", filetype="cif")
```

- `id`: ターゲット蛋白の識別子。`chem.chembl.download_activities`・`chem.rcsb.download_structures`と全く同じ3種類を受け付ける
  - ChEMBLターゲットID (例: `CHEMBL204`)
  - UniProt accession (例: `P00734`)
  - UniProt Entry名/mnemonic (例: `THRB_HUMAN`, `BRAF_HUMAN`)
  - 解決には`chem.ids.resolve_uniprot_accession_any`をそのまま使う(`chem.rcsb`と共用)
- `plddt_thres`: 省略可(デフォルト`None`)。指定時は平均pLDDT信頼度(`globalMetricValue`、0-100)がこの値以上(inclusive)のエントリのみを対象にする。`chem.rcsb`の`resolution_thres`と違い、pLDDTは常に値が存在する(NMRのような「値なし」ケースはない)ため、値なし除外ロジックは不要
- `outdir`: 出力先ディレクトリ。存在しなければ作成する
- `filetype`: `"cif"`(デフォルト)、`"pdb"`、`"both"`のいずれか
- 戻り値: 少なくとも1ファイルが揃った(新規ダウンロードまたは既存ファイル)予測エントリ数

## データ取得方法

1. `resolve_uniprot_accession_any(id)`でUniProt accessionを得る
2. **AlphaFold DB API** (`https://alphafold.ebi.ac.uk/api/prediction/{accession}`、GET)でそのaccessionの予測エントリ一覧を取得する。レスポンスはJSON配列で、通常1件だが以下のケースがある:
   - 非常に大きい蛋白質はフラグメント分割されて複数エントリになることがある
   - 公式のAlphaFold予測がないターゲットではコミュニティ投稿の代替モデル(`providerId`が`GDM`以外)が返る場合があり、その場合`entryId`は`AF-{accession}-F{n}`形式とは限らない(例: `AF-0000000365840311`のような数値ID)
   - HTTPステータス404、または空配列の場合はエントリなしとして`ValueError`
3. 各エントリの構造ファイルURLは、自前でURLを組み立てず**APIレスポンスに含まれる`cifUrl`/`pdbUrl`をそのまま使う**(entryIdの命名規則がケースによって異なるため、URL構築の頑健性のため必須)
4. `plddt_thres`でフィルタリング(`globalMetricValue >= plddt_thres`)した上で、各エントリを`entry["cifUrl"]`/`entry["pdbUrl"]`からダウンロードする(`filetype="both"`なら両方)。`tqdm`で進捗表示(`CHEM_QUIETNESS`が非quietのときのみ)
5. `outdir`に`{entryId}.{ext}`として保存する。**保存先に同名ファイルが既に存在する場合はネットワークアクセスせずスキップする**(`chem.rcsb`と同じ挙動)
6. 該当形式のURLがレスポンスに存在しない場合はその形式だけスキップし(quiet時以外は`stderr`に警告)、処理は継続する

## 呼び出しログと進捗出力

`chem.chembl.download_activities`・`chem.rcsb.download_structures`と同様に、`download_structures`に`@logged`デコレータ(`chem.verbosity.logged`)を適用する。終了時の`"wrote N structures to {outdir}"`メッセージも`is_quiet()`のとき出力しない。

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- 新規の外部依存追加は不要(既存の`requests`、`tqdm`のみ使用)
- `src/chem/alphafold/`として`chem`パッケージのサブパッケージに配置(`pyproject.toml`の`[tool.setuptools.packages.find]`が`chem.alphafold`を自動検出)

## サンプルノートブック

`notebooks/alphafold_pocket_thrb_human.ipynb`の末尾に、`THRB_HUMAN`(トロンビン)を対象としたサンプルセルを追加する:
1. `alphafold.download_structures("THRB_HUMAN", outdir="af_data", filetype="pdb")`を実行(`from chem import alphafold`)
2. `py3Dmol`で予測構造をカートゥーン表示し、B-factor列に格納されているper-residue pLDDT信頼度でカラーリングする(`colorscheme={"prop": "b", "gradient": "roygb", "min": 50, "max": 90}`)

生成される`af_data/`ディレクトリは`.gitignore`に追加し、コミットしない。

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- テストは`tests/test_alphafold_fetch.py`にネットワーク不要な範囲(`_select_entries`のpLDDTフィルタリングロジック、`_download_one`の既存ファイルスキップ、不正な`filetype`のバリデーション)のみ追加する。ネットワークを伴う関数(UniProt/ChEMBL解決、AlphaFold DB問い合わせ・ダウンロード)はユニットテスト対象外とし、手動での実行確認(実際にAPIを叩く)で動作検証する
