# chem

Chemistry utilities for interactive use in Jupyter notebooks.

See [API.md](API.md) for the full function reference.

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

Notebooks under `notebooks/` are committed with outputs stripped. After
cloning, run this once to have `git add`/`git commit` strip a notebook's
outputs and execution counts automatically (the working-tree file on disk is
untouched -- only what's staged/committed is stripped):

```bash
nbstripout --install --attributes .gitattributes
```

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
  standardization/desalting library), `biopython` (structural alignment),
  `numpy`
- `[notebook]` extra: `jupyter`, `jupyterlab`, `notebook`, `ipykernel`,
  `ipywidgets`, `nglview` (fuller-featured 3D/trajectory viewer — installed
  via conda-forge in `environment.yml` for reliable widget asset setup),
  `pandas`
- `[dev]` extra: `pytest`
- AmberTools (`tleap`, `sander`, `antechamber`, `cpptraj`, ...) and `fpocket`
  (pocket detection, used by `chem.protein.find_pocket`): conda-forge only,
  see `environment.yml` — no PyPI package exists for either

## chem.chembl — ChEMBL bioactivity download

```python
from chem import chembl

chembl.download_activities(
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
`False`, one row is written per raw activity record instead. If `output`
already exists, it is left as-is and not re-downloaded (as with
`chem.rcsb`/`chem.alphafold`).

### Progress and call logging

Every call to `chem.chembl.download_activities` prints its function name and
arguments to stderr, and long-running downloads show a `tqdm` progress bar —
both controlled by the `CHEM_QUIETNESS` environment variable. Unset (or set
to `"0"`/`"N"`/`"FALSE"`, case-insensitive) means verbose (the default); any
other value suppresses both.

## chem.rcsb — RCSB PDB structure download

```python
from chem import rcsb

rcsb.download_structures(
    "THRB_HUMAN",  # ChEMBL target id, UniProt accession, or UniProt entry name
    resolution_thres=2.0,  # optional max resolution in Angstrom, inclusive
    outdir="data",  # destination directory, created if missing
    filetype="cif",  # "cif" (default), "pdb", or "both"
)

# or download specific PDB entries directly, skipping target resolution entirely
rcsb.download_structures(["6LU7", "7BQY"], outdir="data")
```

Resolves `id` to a UniProt accession, finds every PDB entry whose polymer entities
are annotated with that accession via the RCSB Search API, and downloads each
qualifying entry's structure file(s) from `files.rcsb.org`. When `resolution_thres`
is set, entries without a resolution (e.g. NMR structures) are excluded; when it's
left as `None`, every entry is downloaded regardless of resolution. Passing a list
(or tuple/set) of PDB entry ids as `id` instead downloads exactly those entries,
with no target resolution or search step. A file already present in `outdir` is
left as-is and not re-downloaded, so re-running only fetches what's missing.

## chem.alphafold — AlphaFold DB predicted structure download

```python
from chem import alphafold

alphafold.download_structures(
    "THRB_HUMAN",  # ChEMBL target id, UniProt accession, or UniProt entry name
    plddt_thres=70.0,  # optional min average pLDDT confidence (0-100), inclusive
    outdir="data",  # destination directory, created if missing
    filetype="cif",  # "cif" (default), "pdb", or "both"
)
```

Resolves `id` to a UniProt accession and downloads every AlphaFold DB prediction
entry for it (usually one, but very large proteins may be split into fragments,
and some targets without an official prediction have community-submitted
alternatives instead), using the download URLs the AlphaFold API itself returns.
When `plddt_thres` is set, only entries whose average pLDDT confidence
(`globalMetricValue`) meets the threshold are kept; when left as `None`, every
entry is downloaded regardless of confidence. As with `chem.rcsb`, a file already
present in `outdir` is left as-is and not re-downloaded.

## chem.protein — structural alignment and pocket detection

```python
from chem import protein

# Sequence-align and structurally superpose a set of same-target structures
# (mix PDB/CIF, RCSB/AlphaFold freely). Writes one PDB file per input into outdir.
rmsd = protein.align(
    ["data/1PPB.cif", "data/1BTH.cif", "af_data/AF-P00734-F1.pdb"],
    reference=None,  # index into the list, or any path (need not be in the list);
                     # defaults to the first entry
    chain=None,  # override auto chain selection (see below)
    outdir="aligned",
)

# Run fpocket on a structure and identify the pocket nearest a ligand.
pocket = protein.find_pocket(
    "aligned/1PPB.pdb",
    ligand=None,  # None = auto-detect from HETATM; or a HET code; or an external file
    outdir="pocket_out",  # keep fpocket's raw output; omit to use a discarded temp dir
)
```

`align` selects each structure's primary polymer chain (the one with the most
standard amino acid residues -- e.g. thrombin's catalytic heavy chain rather than
its short light chain; pass `chain="A"` etc. to override; HETATM residues are
excluded even when their resname matches a standard amino acid, so a covalently-
linked peptidomimetic ligand sharing the protein's chain id doesn't get pulled
into its sequence), sequence-aligns it
against the reference chain, and superposes the whole structure (Kabsch fit on the
sequence-matched CA atoms, applied to every atom including ligands and waters) onto
the reference's frame. Every input, including the reference, is written out as a
PDB file in `outdir` so each can be loaded and overlaid individually (e.g. one
py3Dmol `addModel` call per file) -- the reference does not need to be a member of
the input list; either way it's written to `outdir` exactly once. Returns
`{path: rmsd}`; a structure with too few
residues in common with the reference is skipped with a warning rather than raising.

`find_pocket` runs [fpocket](https://github.com/Discngine/fpocket) on a PDB file
(fpocket requires legacy PDB format, which is what `align` always writes) and picks
the fpocket pocket whose lining atoms are closest to a ligand's 3D coordinates. The
ligand can be auto-detected (the largest non-solvent/ion HETATM group in the file,
e.g. a co-crystallized inhibitor from an RCSB download), a 3-letter PDB HET code to
disambiguate when several ligand-like groups are present, or a path to an external
ligand file (`.pdb`/`.sdf`/`.mol`/`.mol2`, e.g. a docking pose) for structures that
don't contain the ligand themselves. Returns a dict with `pocket_id`,
`score`/`druggability_score`/`volume` convenience fields, the full raw fpocket score
dict under `info`, and `residues` (`[{"chain", "resnum", "resname"}, ...]`) lining
the selected pocket.

## chem.view3d — interactive py3Dmol structure viewer

```python
from chem import view3d
from chem.protein import SOLVENT_AND_IONS

# Rainbow (N -> C) cartoon backbone, plus any HETATM ligand group as sticks,
# followed by a caption. Displays directly -- just call it.
view3d.render_protein(
    "data/1PPB.pdb",
    exclude=SOLVENT_AND_IONS | {"NAG"},  # default: SOLVENT_AND_IONS
    width=600,
    height=500,
)

# Or color by the file's B-factor column instead (e.g. AlphaFold's per-residue
# pLDDT confidence, which is what AlphaFold DB downloads store there).
view3d.render_protein("af_data/AF-P00734-F1.pdb", coloring="bfactor")
```

Cartoon-only styles don't draw ligands, so any HETATM group not in `exclude` is
added as magenta sticks. `coloring="spectrum"` (default) colors the cartoon
N -> C by residue position; `coloring="bfactor"` colors it by the file's
per-atom B-factor column instead, scaled over `bfactor_range` (default `(50,
90)`, AlphaFold's pLDDT confidence convention -- pass the structure's own
B-factor range for crystallographic temperature factors). Below the view, a
caption shows the PDB id, chain ids, ligand HET codes, and experimental
resolution (parsed from the file's `REMARK 2 RESOLUTION` record, or `"N/A"` if
absent -- e.g. NMR/AlphaFold/`protein.align` output). `render_protein`
displays the view and caption directly and returns nothing, so there's no need
to chain `.show()` or use it as a cell's last expression.

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

- `src/chem/` — core package (`verbosity.py`: the `CHEM_QUIETNESS`-aware `@logged` decorator; `ids.py`: shared ChEMBL target id / UniProt accession / UniProt entry name resolution)
- `src/chem/chembl/` — ChEMBL data access (`fetch.py`: `download_activities`, re-exported at package level)
- `src/chem/rcsb/` — RCSB PDB structure download (`fetch.py`: `download_structures`, re-exported at package level)
- `src/chem/alphafold/` — AlphaFold DB structure download (`fetch.py`: `download_structures`, re-exported at package level)
- `src/chem/protein/` — structural tools (`structural_align.py`: `align`; `pocket.py`: `find_pocket`; both re-exported at package level)
- `src/chem/view3d/` — interactive structure viewing (`render.py`: `render_protein`, re-exported at package level)
- `notebooks/` — example notebooks
- `tests/` — pytest test suite

## Test

```bash
pytest
```
