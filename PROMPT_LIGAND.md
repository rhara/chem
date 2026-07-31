# chem.ligand の再現プロンプト

`chem`リポジトリの`chem`パッケージ内に、ダウンロード済み構造からリガンドを抽出して結合次数を修正し、薬らしさ(QED)を評価するためのサブパッケージ`ligand`(`chem.ligand`)を追加する指示。`chem.protein`の`SOLVENT_AND_IONS`を再利用する。

## 背景

`chem.rcsb.download_structures`でダウンロードしたPDBファイルのHETATM座標には結合次数(単結合/二重結合/芳香環)の情報が無く、RDKitで素直に読み込むと全て単結合・非芳香族として解釈される。リガンドの3-letterコード(PDB Chemical Component Dictionary, CCD)を手がかりに、RCSBが公開しているそのコードの理想SMILES(結合次数が正しい)をテンプレートとして`AllChem.AssignBondOrdersFromTemplate`で結合次数を復元し、得られた正しい分子でQEDなどの記述子を計算できるようにする。

## 要件

`src/chem/ligand/extract.py`に以下を実装し、`src/chem/ligand/__init__.py`で`from .extract import list_ligand_codes, list_ligand_instances, load_ligand, molecular_weight, qed`として再エクスポートする。以下の形で呼び出せること:

```python
from chem import ligand

instances = ligand.list_ligand_instances("data/3RM0.pdb")
# [{"code": "S54", "chain": "H", "resnum": 1, "icode": ""},
#  {"code": "S54", "chain": "H", "resnum": 2, "icode": ""}]
codes = ligand.list_ligand_codes("data/3RM0.pdb")        # ["S54"]
mol = ligand.load_ligand("data/3RM0.pdb", "S54", chain="H", resnum=1)  # 結合次数修正済みRDKit Mol
ligand.molecular_weight(mol)                               # 499.6
ligand.qed(mol)                                             # 0.29
```

### `list_ligand_instances(structure_path, exclude=SOLVENT_AND_IONS)`

- `structure_path`の全HETATM残基(`exclude`、デフォルト`chem.protein.SOLVENT_AND_IONS`、に含まれるものを除く)を**物理的な出現ごとに1件**列挙し、`{"code", "chain", "resnum", "icode"}`の辞書のリストをファイル順で返す
- 同じコードが複数箇所(複数チェーン、あるいは同一チェーン内の異なるresnum)にあってもまとめず、出現ごとに別のエントリとして返す(例: ホモ二量体の各サブユニットに同じ阻害剤が1分子ずつ結合していれば2件)
- `chain`/`resnum`/`icode`は`load_ligand`に渡してその特定の1コピーを指定するのに使う

### `list_ligand_codes(structure_path, exclude=SOLVENT_AND_IONS)`

- `list_ligand_instances`の結果からコードだけを取り出し、重複を除いてソートして返す(同じコードの複数コピーは1つにまとまる)
- 呼び出し側で`SOLVENT_AND_IONS | {"NAG", "TYS", "MRD"}`のように追加のコードを渡せば、データセット固有の非リガンドHETATM(糖鎖修飾・修飾残基・結晶化添加剤など)も除外できる(このデータセット固有の追加除外はnotebook側の責務とし、関数のデフォルトには含めない)

### `load_ligand(structure_path, ligand, chain=None, resnum=None, icode=None)`

- `ligand`は3文字のPDB HETコード(`list_ligand_codes`/`list_ligand_instances`が返すもの)
- `chain`/`resnum`/`icode`(`list_ligand_instances`の1エントリの値をそのまま渡す想定)を指定すると、そのコードが複数出現する場合に特定の1コピーへ絞り込む。全て`None`(デフォルト)なら絞り込みを行わず、マッチした中で最も原子数の多い(=最も完全な)インスタンスを使う(化学的に同一なコピー間ではどれでもよい場合向け)
- 内部処理:
  1. Bio.PDBで`structure_path`をパースし、resnameが`ligand`(大文字小文字を無視)に一致するHETATM残基を集める。`chain`/`resnum`/`icode`が指定されていればそれぞれで絞り込む。絞り込んだ結果が空なら`ValueError`。複数残る場合は原子数最大のものを使う
  2. その残基だけを含む一時PDBファイルを`Bio.PDB.PDBIO`(+ `Select`サブクラスで対象残基のみ`accept_residue`)で書き出し、`Chem.MolFromPDBFile(path, sanitize=False, removeHs=False)`でRDKit Molとして読み込む。**注意**: RDKitの`removeHs=True`は`sanitize=True`のときしか効かないため(`sanitize=False`のままだと無視される)、読み込み後に明示的に`Chem.RemoveHs(mol, sanitize=False)`でHを除去すること
  3. RCSBの Chemical Component Dictionary REST API (`https://data.rcsb.org/rest/v1/core/chemcomp/{CODE}`) を叩き、`pdbx_chem_comp_descriptor`から理想SMILESの候補を集める。優先順位: (OpenEye OEToolkits, SMILES_CANONICAL) → (OpenEye OEToolkits, SMILES) → (CACTVS, SMILES_CANONICAL) → (CACTVS, SMILES) → その他の`SMILES`/`SMILES_CANONICAL`型全て(順不同、重複除去)
  4. 候補ごとに`Chem.MolFromSmiles`でテンプレートを作り、`Chem.RemoveHs`で明示的Hを落とす(一部のCCD SMILESは立体配置を示すためのH原子を`[H]/N=C(...)`のように明示的原子として書いており、そのままだと抽出した分子と原子数が合わない)。テンプレートの重原子数が抽出分子と一致しないものはスキップし、一致するものだけ`AllChem.AssignBondOrdersFromTemplate(template, mol)`を試す。成功したら`Chem.SanitizeMol`して返す
  5. どの候補も一致・成功しなければ`ValueError`(共有結合で繋がった複数残基からなるペプチド様リガンドの一部を単独残基として抽出した場合や、電子密度が不完全な残基などで起こりうる)
- 戻り値: 結合次数・芳香族性が修正され、sanitize済みのRDKit Mol

### `qed(mol)` / `molecular_weight(mol)`

- それぞれ`rdkit.Chem.QED.qed(mol)` / `rdkit.Chem.Descriptors.MolWt(mol)`への薄いラッパー

## 前提環境

- `~/chem`リポジトリ、`chem` conda-forge環境(Python 3.12)
- 新規の外部依存追加は不要(既存の`rdkit`、`biopython`、`requests`のみ使用)
- `src/chem/ligand/`として`chem`パッケージのサブパッケージに配置(`pyproject.toml`の`[tool.setuptools.packages.find]`が`chem.ligand`を自動検出)

## サンプルノートブックの更新

`notebooks/alphafold_pocket_thrb_human.ipynb`の「RCSB PDB structures for thrombin」ダウンロードセルの直後(構造ビューアセルの前)に新セクションを追加する:

- `chem.protein.SOLVENT_AND_IONS | {"NAG", "TYS", "MRD"}`を除外セットとして、`data/`の全PDBファイルに対し`list_ligand_instances`でリガンドの出現を列挙(同一コードの複数コピーも含め、出現ごとに1件)
- 各インスタンスを`load_ligand(path, inst["code"], chain=inst["chain"], resnum=inst["resnum"], icode=inst["icode"])`で抽出し、失敗した場合は`try`/`except`で捕捉して理由付きでスキップリストに積む(結合次数テンプレートが一致しない共有結合ペプチド様リガンドの断片や、密度不完全な残基などが該当しうる)
- 成功した分子について`molecular_weight`/`qed`を計算し、`pandas.DataFrame`で「pdb_id / ligand / chain / resnum / mw / qed」の列を持つ表にまとめて`pdb_id`→`ligand`→`chain`→`resnum`の順でソートして表示する(`chain`/`resnum`列は、同一コードの複数コピーが同じ見た目の行として区別できなくなるのを防ぐため)

## 注意

- このリポジトリはdd_*プロジェクト群、`~/lab/chembl`、`dd_chembl`とは無関係な独立プロジェクト。それらのコードやロジックを参照・流用しない
- `src/chem/view3d/render.py`の`_template_mol`/`_ligand_molblock`も同様にCCDの`SMILES_CANONICAL`を使った結合次数復元を行っているが、そちらは表示用に候補1つのみ・失敗時は`None`を返して静かにフォールバックする実装であり、本モジュールとは目的が異なるため統合しない(ただし同じ`removeHs=True`+`sanitize=False`の組み合わせによる無効化の可能性があるため、修正するなら別タスクとして扱う)
- テストは`tests/test_ligand_extract.py`にネットワーク不要な範囲(`list_ligand_codes`/`list_ligand_instances`の除外ロジック、同一コードの複数コピーが`list_ligand_codes`では1件にまとまり`list_ligand_instances`では別々のエントリのまま残ること、`_pick_ligand_residue`が見つからない場合や`chain`指定が一致しない場合に`ValueError`になること、`chain`指定で複数コピーのうち特定の1つを選べること、`_fetch_template_smiles_candidates`の優先順位・重複除去(`requests.get`を`monkeypatch`で差し替え)、ベンゼン環をRDKitで生成して結合次数無しのPDBとして書き出し・CCD応答を`monkeypatch`で偽装した上で`load_ligand`が芳香族性を正しく復元すること、原子数が合わないテンプレートでは`ValueError`になること、`qed`/`molecular_weight`が既知の分子に対して妥当な範囲の値を返すこと)を追加する
