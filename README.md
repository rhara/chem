# chem

Chemistry utilities for interactive use in Jupyter notebooks.

Currently a minimal package skeleton — functionality is added incrementally.

## Install

AmberTools has no PyPI package (conda-forge or source build only), so the
environment — including AmberTools and RDKit — is created from
`environment.yml` via conda-forge, not from `pyproject.toml` alone:

```bash
mamba env create -f environment.yml
mamba activate chem
python -m ipykernel install --user --name chem --display-name "Python 3.12 (chem)"
```

Then select the "Python 3.12 (chem)" kernel in Jupyter.

If you only need the pure-Python parts (RDKit, no AmberTools/`tleap`/`sander`/
`antechamber`/`cpptraj`), you can instead do a plain venv + pip install:

```bash
python3.12 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,notebook]"
```

## Dependencies

- Core: `rdkit`, `py3dmol` (lightweight in-notebook 3D viewer), `tqdm`,
  `requests`, `chembl_structure_pipeline` (ChEMBL's own structure
  standardization/desalting library)
- `[notebook]` extra: `jupyter`, `jupyterlab`, `notebook`, `ipykernel`,
  `ipywidgets`, `nglview` (fuller-featured 3D/trajectory viewer — installed
  via conda-forge in `environment.yml` for reliable widget asset setup),
  `pandas`
- `[dev]` extra: `pytest`
- AmberTools (`tleap`, `sander`, `antechamber`, `cpptraj`, ...): conda-forge
  only, see `environment.yml` — no PyPI package exists

## chem.chembl.fetch — ChEMBL bioactivity download

```python
import chem.chembl.fetch as cf

cf.download_activities(
    "BRAF_HUMAN",  # ChEMBL target id, UniProt accession, or UniProt entry name
    mw=[250, 650],  # optional [lower, upper] molecular weight filter, inclusive
    normalize_smiles=True,  # standardize/desalt + aggregate duplicate compounds
    output="activities.tsv",  # ".csv" -> comma-delimited, otherwise tab-delimited
)
```

Talks to the ChEMBL and UniProt REST APIs directly (no
`chembl_webresource_client`). Only activities with a pChEMBL value are kept.
When `normalize_smiles=True`, compounds are standardized/desalted via the
[ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline)
and duplicate compounds (same resulting smiles) are collapsed into one row
with `n`/`pchembl_mean`/`pchembl_median`/`pchembl_std` plus a
`parent_chembl_id` pointing at one representative `molecule_chembl_id`; when
`False`, one row is written per raw activity record instead.

### Progress and call logging

Every call to `chem.chembl.fetch.download_activities` prints its function name and
arguments to stderr, and long-running downloads show a `tqdm` progress bar —
both controlled by the `CHEM_QUIETNESS` environment variable. Unset (or set
to `"0"`/`"N"`/`"FALSE"`, case-insensitive) means verbose (the default); any
other value suppresses both.

### A note on `conda activate` and plain `python`/`pip`

If your shell auto-activates another conda env at startup (e.g. via
`auto_activate_base` or a login script), `PATH` ordering after
`mamba activate chem` is not always reliable — `CONDA_PREFIX` will correctly
point at `chem`, but a bare `python`/`pip` call can still resolve to a
different env's binary. When in doubt, prefer explicit paths or `conda run`:

```bash
$CONDA_PREFIX/bin/python -c "import rdkit; print(rdkit.__version__)"
# or
conda run -n chem python -c "import rdkit; print(rdkit.__version__)"
```

## Layout

- `src/chem/` — core package (`verbosity.py`: the `CHEM_QUIETNESS`-aware `@logged` decorator)
- `src/chem/chembl/` — ChEMBL data access (`fetch.py`: `download_activities`)
- `notebooks/` — example notebooks
- `tests/` — pytest test suite

## Test

```bash
pytest
```
