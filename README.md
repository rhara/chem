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

- Core: `rdkit`, `py3dmol` (lightweight in-notebook 3D viewer), `tqdm`
- `[notebook]` extra: `jupyter`, `jupyterlab`, `notebook`, `ipykernel`,
  `ipywidgets`, `nglview` (fuller-featured 3D/trajectory viewer — installed
  via conda-forge in `environment.yml` for reliable widget asset setup)
- `[dev]` extra: `pytest`
- AmberTools (`tleap`, `sander`, `antechamber`, `cpptraj`, ...): conda-forge
  only, see `environment.yml` — no PyPI package exists

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

- `src/chem/` — package source
- `notebooks/` — example notebooks
- `tests/` — pytest test suite

## Test

```bash
pytest
```
