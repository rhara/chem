# chem.blast.blastp の再現プロンプト

`chem`パッケージ内に、BLAST検索用のサブパッケージ`blast`(`chem.blast`)を追加し、EBI Job Dispatcher REST API経由でblastp検索を実行しヒット一覧を返す関数を実装するための指示。他のサブパッケージ(`chem.rcsb`/`chem.alphafold`)と異なり外部データベースへのダウンロード目的ではなく、配列の類似検索そのものが目的。

## 要件

`src/chem/blast/search.py` に関数`blastp`を実装し、`src/chem/blast/__init__.py`で`from .search import blastp`として再エクスポートする。以下の形で呼び出せること:

```python
from chem import blast
hits = blast.blastp(sequence, database="pdb", email="you@example.com")
```

- `sequence`: 検索する蛋白質配列(プレーン配列でもFASTA形式でもよい。EBIのAPIはどちらも受け付ける)
- `database`: 検索対象のEBIホストデータベース名(例: `"pdb"` — PDBの全チェーン、ヒットはチェーン単位、例`"1ABC_A"`。`"uniprotkb_swissprot"` — レビュー済みUniProtエントリ、蛋白質単位で1ヒット)
- `email`: EBI Job Dispatcher APIが要求する連絡先メールアドレス(必須引数、デフォルト値なし)
- `matrix`: 置換行列。デフォルト`"BLOSUM62"`
- `expect`: E-value閾値(この値以下、inclusive)。デフォルト`1e-10`
- `max_hits`: 取得する最大ヒット数。デフォルト`50`(EBI側の1検索あたりの上限でもある)
- `poll_interval`: ジョブ完了待ちのポーリング間隔(秒)。デフォルト`10`
- `timeout`: ジョブ完了待ちのタイムアウト(秒)。デフォルト`600`
- `title`: EBI側に記録されるジョブタイトル。デフォルト`"chem.blast"`
- 戻り値: ヒットごとの辞書のリスト(EBI自身のランキング順、スコア良い順)。各辞書は`accession`(ヒットのデータベースaccession、`database="pdb"`なら`"1ABC_A"`のようなチェーンID、`database="uniprotkb_swissprot"`ならUniProt accession)、`description`(ヒットの説明文全体)、`identity_pct`(アライメント領域での配列同一性%)、`align_len`(アライメント長、残基数)、`evalue`(E-value)を持つ

## データ取得方法

[EBI Job Dispatcher REST API](https://www.ebi.ac.uk/jdispatcher/docs/webservices/)(`https://www.ebi.ac.uk/Tools/services/rest/ncbiblast`)を直接叩く(専用クライアントライブラリは使わない):

1. `POST {base}/run` に `email`/`program=blastp`/`stype=protein`/`sequence`/`database`/`matrix`/`exp`/`alignments`(=`max_hits`)/`scores`(=`max_hits`)/`title` をフォームデータとして送信し、レスポンス本文(ジョブID文字列、前後の空白を`strip`)を受け取る
2. `GET {base}/status/{job_id}` を`poll_interval`秒間隔でポーリングし、状態が`FINISHED`/`FAILURE`/`ERROR`/`NOT_FOUND`のいずれか(終端状態)になるまで繰り返す。`timeout`秒を超えたら`TimeoutError`を送出する。ポーリング中は`tqdm`で進捗表示(`is_quiet()`のときは非表示) — 所要時間が事前にわからないため、単純なブロッキング待機ではなく必ず進捗を出す
3. 終端状態が`FINISHED`以外なら`RuntimeError`を送出する
4. `GET {base}/result/{job_id}/json` でJSON形式の結果を取得し、`hits`配列の各要素から`hit_acc`→`accession`、`hit_def`→`description`、`hit_hsps[0].hsp_identity`→`identity_pct`(小数点1桁に丸め)、`hit_hsps[0].hsp_align_len`→`align_len`、`hit_hsps[0].hsp_expect`→`evalue`を取り出しリストにする

## 呼び出しログについて(他サブパッケージとの違い)

このリポジトリの公開関数は原則`chem.verbosity.logged`デコレータを付与し呼び出し引数をstderrにログ出力する規約だが、`blastp`は**あえて`@logged`を付けない**。`email`引数がそのままログに出力されてしまい、呼び出し元の連絡先メールアドレスをノートブックの出力やログに無意味に残してしまうため。代わりに、`is_quiet()`でないときだけ`"{job_id}: {N} hit(s) against {database!r}"`という1行サマリーをstderrに出力する(引数の`email`/`sequence`はログに含めない)。この理由をdocstringおよびAPI.md/README.mdに明記すること。

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- 新規の外部依存追加は不要(既存の`requests`、`tqdm`のみ使用)
- `src/chem/blast/`として`chem`パッケージのサブパッケージに配置(`pyproject.toml`の`[tool.setuptools.packages.find]`が`chem.blast`を自動検出)

## サンプルノートブック

`notebooks/cdk20_similar_targets.ipynb`の「BLASTP search via EBI Job Dispatcher」節を、素のPOST/ポーリングコードのベタ書きから`chem.blast.blastp`の呼び出しに置き換える。`email`は環境変数ではなくノートブック内に直接記述する(リポジトリ所有者からの明示指示)。

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- テストは`tests/test_blast_search.py`にネットワーク不要な範囲(`_parse_hits`のJSONパース、`_submit`が送るPOSTパラメータの検証、`_status`のレスポンス文字列のstrip、`_wait`のポーリングループが終端状態で停止すること・タイムアウトで`TimeoutError`を送出すること、`blastp`本体が成功時にヒットリストを返すこと・ジョブ失敗時に`RuntimeError`を送出すること)のみ追加する。EBIへの実際のネットワークアクセスはユニットテスト対象外とし、手動での実行確認(実際にAPIを叩く)で動作検証する
