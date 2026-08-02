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

## chem.blast — BLAST search via EBI's Job Dispatcher

```python
from chem import blast

hits = blast.blastp(
    sequence,  # a protein sequence, plain or FASTA-formatted
    database="pdb",  # e.g. "pdb" (chain-level hits with a solved structure) or "uniprotkb_swissprot"
    email="you@example.com",  # required by EBI's Job Dispatcher API
)
# [{"accession": "1ABC_A", "description": "...", "identity_pct": 44.8, "align_len": 288, "evalue": 7e-82}, ...]
```

Submits a `blastp` search to [EBI's Job Dispatcher](https://www.ebi.ac.uk/jdispatcher/docs/webservices/),
polls until it finishes (progress shown via `tqdm`), and returns the ranked
hit list — handy for finding a data-rich homolog (more PDB structures, more
ChEMBL activity data) when the actual target of interest has little of
either. Unlike every other `chem` function, `blastp` is deliberately **not**
decorated with the call-logging described below, since that would echo the
caller's `email` to stderr; `CHEM_QUIETNESS` still controls its progress bar
and summary line.

## chem.protein — structural alignment and pocket detection

```python
from chem import protein

# Sequence-align and structurally superpose a set of same-target structures
# (mix PDB/CIF, RCSB/AlphaFold freely). Writes one PDB file per input into outdir.
# Returns {path: {"rmsd": ..., "identity": ...}}.
align_results = protein.align(
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

# Or, for a structure with no bound ligand at all (e.g. an AlphaFold
# prediction), list every candidate pocket fpocket finds instead.
pockets = protein.list_pockets("af_data/AF-P00734-F1.pdb")

# Split a structure into a ligand-free protein PDB (one per chain) and one SDF
# per non-water HETATM ligand instance (ions/sugars/additives included) -- e.g.
# to prep a receptor/ligand pair for docking.
split_result = protein.split(
    "data/1R1H.pdb",
    all_chains=False,  # True writes a single protein PDB with every chain together instead
    remove_water=False,  # True also strips water out of the protein PDB
    outdir="split",
)
# split_result["protein"] -> {"A": "split/1R1H_protein_A.pdb"}
# split_result["ligands"] -> [{"path": "split/1R1H_ligand_BIR_A2001.sdf", "code": "BIR",
#                               "bond_orders_restored": True, ...}, ...]
```

`align` selects each non-reference structure's chain by whichever best matches the
reference (most identical residues at gap-free matched positions) -- not simply the
chain with the most residues, since a structure can contain a larger bound partner
protein (e.g. thrombin co-crystallized with a serpin inhibitor) that a size-only
heuristic would wrongly prefer (pass `chain="A"` etc. to override; the reference's
own chain is still picked by size alone, so use an unambiguous single-chain
reference, e.g. an AlphaFold prediction, if in doubt). HETATM residues are excluded
even when their resname matches a standard amino acid, so a covalently-linked
peptidomimetic ligand sharing the protein's chain id doesn't get pulled into its
sequence. It then sequence-aligns the selected chain
against the reference chain, and superposes the whole structure (Kabsch fit on the
sequence-matched CA atoms, applied to every atom including ligands and waters) onto
the reference's frame. Every input, including the reference, is written out as a
PDB file in `outdir` so each can be loaded and overlaid individually (e.g. one
py3Dmol `addModel` call per file) -- the reference does not need to be a member of
the input list; either way it's written to `outdir` exactly once. Returns
`{path: {"rmsd": ..., "identity": ...}}` -- `identity` is the fraction of
matched (gap-free) sequence positions with an identical residue, so a
mismatch still counts as "matched" (identity can be `< 1.0` even with no
gaps at all), but a gapped position (e.g. a loop present in one structure but
not the other) counts in neither the numerator nor the denominator. A
structure with too few residues in common with the reference is skipped with
a warning rather than raising.

`find_pocket` runs [fpocket](https://github.com/Discngine/fpocket) on a PDB file
(fpocket requires legacy PDB format, which is what `align` always writes) and picks
the fpocket pocket whose lining atoms are closest to a ligand's 3D coordinates. The
ligand can be auto-detected (the largest non-solvent/ion HETATM group in the file,
e.g. a co-crystallized inhibitor from an RCSB download), a 3-letter PDB HET code to
disambiguate when several ligand-like groups are present, or a path to an external
ligand file (`.pdb`/`.sdf`/`.mol`/`.mol2`, e.g. a docking pose) for structures that
don't contain the ligand themselves. Returns a dict with `pocket_id`,
`score`/`druggability_score`/`volume` convenience fields, the full raw fpocket score
dict under `info`, `residues` (`[{"chain", "resnum", "resname"}, ...]`) lining
the selected pocket, and `spheres` (`[{"x", "y", "z", "radius"}, ...]`) -- fpocket's
own alpha spheres approximating the pocket cavity's shape, handy for rendering the
pocket as a filled volume (e.g. one py3Dmol `addSphere` per entry) instead of sticks
on the lining residues.

`list_pockets` runs the same fpocket analysis but, instead of picking the one
pocket nearest a ligand, returns every pocket fpocket detected -- for structures
with no bound ligand to anchor on. Each entry is shaped exactly like a single
`find_pocket` result, sorted by `druggability_score` descending. fpocket routinely
reports dozens of low-quality cavities on a typical structure, so `druggability_thres`
(default `0.1`) drops any pocket scoring below it, or with no score at all;
pass `None` to keep everything unfiltered.

`split` decomposes a structure file into a ligand-free protein (water kept by
default, everything else HETATM stripped) and one SDF molecule per non-water HETATM
residue instance -- real ligands, ions (`ZN`), glycosylation sugars (`NAG`),
crystallization additives, all of it, every single one. Each is written out
via `chem.ligand.load_ligand` when possible (bond orders/aromaticity restored
against the PDB Chemical Component Dictionary, `bond_orders_restored=True`
in the returned entry); when `load_ligand` can't resolve a template for it
(e.g. a covalently-linked glycosylation sugar or peptidomimetic ligand
missing the atom(s) involved in that link, relative to the free/standalone
template, or incomplete crystallographic density), a warning is printed and
the raw single-bond connectivity RDKit guesses from atomic distances is
written instead (`bond_orders_restored=False`) -- every instance still gets
an SDF either way, just not always with corrected bonds. By default the
protein PDB is written one file per chain, the chain id folded into the
filename; `all_chains=True` writes a single file with every chain together
instead. `remove_water=True` additionally strips water out of the protein
PDB (off by default, since crystallographic waters are routinely useful
downstream).

## chem.ligand — extract a ligand and score its drug-likeness

```python
from chem import ligand

# Every distinct non-solvent/ion HETATM *occurrence* in the file -- two copies
# of the same ligand code (e.g. one per chain) come back as two entries.
instances = ligand.list_ligand_instances("data/3RM0.pdb")
# [{"code": "S54", "chain": "H", "resnum": 1, "icode": ""},
#  {"code": "S54", "chain": "H", "resnum": 2, "icode": ""}]

# Extract one as a proper RDKit molecule: real 3D coordinates, plus bond
# orders/aromaticity fixed up against the PDB Chemical Component Dictionary's
# ideal SMILES for that code (PDB files carry no bond-order information, so
# RDKit's raw read is otherwise all single bonds). chain/resnum/icode pin
# down a specific copy when the code isn't unique in the file.
inst = instances[0]
mol = ligand.load_ligand("data/3RM0.pdb", inst["code"], chain=inst["chain"], resnum=inst["resnum"])

ligand.molecular_weight(mol)  # 499.6
ligand.qed(mol)  # 0.29 -- Quantitative Estimate of Drug-likeness, 0-1
```

A residue that's really one piece of a covalently-linked multi-residue
ligand (e.g. a peptidomimetic inhibitor built from linked D-amino acid
HETATM groups), or has incomplete crystallographic density, can't be matched
to a standalone template -- `load_ligand` raises `ValueError` rather than
guessing, so batch over `list_ligand_instances` with a `try`/`except` to skip
and report those.

## chem.view3d — interactive py3Dmol structure viewer

```python
from chem import view3d
from chem.protein import SOLVENT_AND_IONS

# Rainbow (N -> C) cartoon backbone, plus any HETATM ligand group as sticks,
# followed by a caption. Displays directly -- just call it. Default exclude is
# chem.protein.WATER (just water), so bound ions/additives show up too.
view3d.render_protein(
    "data/1PPB.pdb",
    exclude=SOLVENT_AND_IONS | {"NAG"},  # go narrower than the WATER default if wanted
    width=600,
    height=500,
)

# Or color by the file's B-factor column instead (e.g. AlphaFold's per-residue
# pLDDT confidence, which is what AlphaFold DB downloads store there).
view3d.render_protein("af_data/AF-P00734-F1.pdb", coloring="bfactor")

# Or render the backbone as a solid volume instead of a cartoon ribbon --
# a translucent van der Waals surface (union of a smooth blob at every atom,
# not individual spheres), so a bound ligand can still show through it.
view3d.render_protein("af_data/AF-P00734-F1.pdb", coloring="bfactor", style="surface")
```

Cartoon/surface styles don't draw ligands, so any HETATM group not in
`exclude` is added as magenta sticks. `coloring="spectrum"` (default) colors
the backbone N -> C by residue position; `coloring="bfactor"` colors it by
the file's per-atom B-factor column instead, scaled over `bfactor_range`
(default `(50, 90)`, AlphaFold's pLDDT confidence convention -- pass the
structure's own B-factor range for crystallographic temperature factors).
`style="cartoon"` (default) draws a ribbon; `style="surface"` draws a solid,
translucent (opacity 0.85) van der Waals volume computed by 3Dmol.js via
marching cubes, restricted to non-HETATM atoms so a ligand keeps its own
stick rendering rather than being swallowed by the volume. The view sits in a
light-gray border -- exactly the area where 3Dmol.js's mouse controls
(rotate/zoom/pan) take over -- with a caption beside it on the right, one
field per line: the PDB id, chain ids, ligand HET codes, and experimental
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
- `src/chem/blast/` — BLAST search via EBI's Job Dispatcher (`search.py`: `blastp`, re-exported at package level)
- `src/chem/protein/` — structural tools (`structural_align.py`: `align`; `pocket.py`: `find_pocket`; both re-exported at package level)
- `src/chem/ligand/` — ligand extraction and scoring (`extract.py`: `list_ligand_codes`, `list_ligand_instances`, `load_ligand`, `qed`, `molecular_weight`, all re-exported at package level)
- `src/chem/view3d/` — interactive structure viewing (`render.py`: `render_protein`, re-exported at package level)
- `notebooks/` — example notebooks
- `tests/` — pytest test suite

## Test

```bash
pytest
```
