# chem

Chemistry utilities for interactive use in Jupyter notebooks.

Currently a minimal package skeleton — functionality is added incrementally.

## Install

```bash
mamba create -n chem python=3.12 -y
mamba activate chem
pip install -e ".[dev]"
python -m ipykernel install --user --name chem --display-name "Python 3.12 (chem)"
```

Then select the "Python 3.12 (chem)" kernel in Jupyter.

## Layout

- `src/chem/` — package source
- `notebooks/` — example notebooks
- `tests/` — pytest test suite

## Test

```bash
pytest
```
