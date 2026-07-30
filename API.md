# chem API reference

Public functions in the `chem` package, grouped by subpackage. For install
instructions and quick-start examples see [README.md](README.md); for the
original implementation spec of each function see the matching
`PROMPT_*.md` file.

Every function below is decorated with `chem.verbosity.logged`: unless
`CHEM_QUIETNESS` is set to something other than `"0"`/`"N"`/`"FALSE"`
(case-insensitive), each call prints its function name and fully-bound
arguments to stderr, and any `tqdm` progress bars are shown. Setting
`CHEM_QUIETNESS` to any other value (e.g. `"1"`) silences both.

## chem.chembl

### `chembl.download_activities(id, mw=None, normalize_smiles=False, output="activities.tsv")`

Download ChEMBL bioactivity records for a target into a tsv/csv file. Talks to
the ChEMBL and UniProt REST APIs directly (no `chembl_webresource_client`).

- `id` — ChEMBL target id (e.g. `"CHEMBL204"`), UniProt accession (e.g.
  `"P00734"`), or UniProt entry name (e.g. `"THRB_HUMAN"`).
- `mw` — optional `[lower, upper]` molecular weight range, inclusive on both
  ends.
- `normalize_smiles` — if `True`, each compound is standardized/desalted via
  the [ChEMBL Structure Pipeline](https://github.com/chembl/ChEMBL_Structure_Pipeline)
  (its parent compound is used for the `mw` filter), and duplicate compounds
  (same resulting smiles) are aggregated into one row with
  `n`/`pchembl_mean`/`pchembl_median`/`pchembl_std` plus a `parent_chembl_id`
  pointing at one representative `molecule_chembl_id`. If `False`, one row is
  written per raw activity record instead.
- `output` — destination file path; delimiter is chosen from the extension
  (`.csv` → comma, otherwise tab). **If this file already exists, it is left
  as-is and not re-downloaded.**

Only activities with a pChEMBL value are kept. Returns the number of rows
written (or already present, if `output` existed).

## chem.rcsb

### `rcsb.download_structures(id, resolution_thres=None, outdir="data", filetype="cif")`

Download RCSB PDB structure files for a target, or an explicit list of PDB
entries, into a directory. When `id` is a ChEMBL target id, UniProt accession,
or UniProt entry name, it's resolved to a UniProt accession, and every PDB
entry whose polymer entities are annotated with that accession via the RCSB
Search API is downloaded. When `id` is instead a list (or tuple/set) of PDB
entry ids (e.g. `["6LU7", "7BQY"]`), those entries are downloaded directly,
skipping target resolution and search entirely. Either way, downloads come
from `files.rcsb.org`, and the resulting entries can optionally be filtered
by `resolution_thres` regardless of which form `id` took.

- `id` — ChEMBL target id, UniProt accession, or UniProt entry name (same
  three forms as `chembl.download_activities`), or a list (or tuple/set) of
  PDB entry ids to download directly. Explicit PDB ids are case-insensitive
  (normalized to uppercase); a `ValueError` is raised if the list is empty or
  any entry isn't shaped like a PDB id (4 characters: a digit followed by
  three alphanumerics).
- `resolution_thres` — optional maximum resolution in Å, inclusive. When set,
  entries without a resolution (e.g. NMR structures) are excluded; when
  `None` (default), every entry is kept regardless of resolution.
- `outdir` — destination directory; created if missing.
- `filetype` — `"cif"` (default), `"pdb"`, or `"both"`.

**A file already present in `outdir` (same entry id and extension) is left
as-is and not re-downloaded**, so repeated calls only fetch what's missing.
Returns the number of entries for which at least one file is present (newly
downloaded or already on disk).

## chem.alphafold

### `alphafold.download_structures(id, plddt_thres=None, outdir="data", filetype="cif")`

Download AlphaFold DB predicted structure files for a target into a
directory. Resolves `id` to a UniProt accession and downloads every
AlphaFold DB prediction entry for it — usually one, but very large proteins
may be split into fragments, and some targets without an official prediction
have community-submitted alternatives instead — using the download URLs the
AlphaFold API itself returns (robust to non-standard entry ids).

- `id` — same three forms as above.
- `plddt_thres` — optional minimum average pLDDT confidence (0-100),
  inclusive (`globalMetricValue`). `None` (default) keeps every entry
  regardless of confidence.
- `outdir` — destination directory; created if missing.
- `filetype` — `"cif"` (default), `"pdb"`, or `"both"`.

Same skip-if-exists behavior as `chem.rcsb`. Returns the number of entries
for which at least one file is present.

## chem.protein

### `protein.align(structures, reference=None, chain=None, outdir="aligned")`

Sequence-align and structurally superpose a set of same-target structures
(PDB/CIF freely mixed — combine `chem.rcsb`/`chem.alphafold` downloads).

- `structures` — list of PDB/CIF file paths for the same protein.
- `reference` — which structure to align everything onto: an index into
  `structures`, or a path. Defaults to `structures[0]`. Does not need to be a
  member of `structures`; either way it's written to `outdir` exactly once.
- `chain` — optional chain id to use in every structure, overriding
  auto-selection. Default: each structure's chain with the most standard
  amino acid residues (its primary polymer chain — e.g. thrombin's catalytic
  heavy chain rather than its short light chain). HETATM residues are excluded
  even when their resname matches a standard amino acid, so a covalently-linked
  peptidomimetic ligand sharing the protein's chain id (e.g. a D-amino-acid-
  containing inhibitor) doesn't get pulled into the sequence.
- `outdir` — destination directory; created if missing.

For each non-reference structure, its selected chain is sequence-aligned
against the reference chain (`Bio.Align.PairwiseAligner`, global mode), and
the whole structure is superposed via a Kabsch fit (`Bio.PDB.Superimposer`)
on the sequence-matched CA atoms, applied to every atom including ligands and
waters. Every structure, including the reference, is written out as a PDB
file in `outdir` (regardless of input format) — ready to be loaded and
overlaid, e.g. one py3Dmol `addModel` call per file, or as input to
`chem.protein.find_pocket`.

Returns `{path: rmsd}` over the sequence-matched CA atoms for every structure
that could be aligned (the reference maps to `0.0`). A structure with no
usable chain, or too few residues in common with the reference, is skipped
with a warning (unless quiet) rather than raising.

Large structures (e.g. cryo-EM assemblies with >26 chains or >99999 atoms)
are not supported by the legacy PDB writer used here.

### `protein.find_pocket(structure, ligand=None, outdir=None)`

Run [fpocket](https://github.com/Discngine/fpocket) on a PDB file and
identify the pocket nearest a ligand, along with its lining residues.

- `structure` — path to a PDB file. fpocket requires legacy PDB format;
  `chem.protein.align`'s output is a natural fit.
- `ligand` — how to locate the ligand:
  - `None` (default) — auto-detect the largest non-solvent/ion HETATM group
    in `structure` (e.g. a co-crystallized inhibitor from an RCSB download).
    This is a size heuristic, not a drug-likeness check — a glycosylation
    sugar (`NAG`) or a modified peptide residue (`TYS`) can outsize a small
    inhibitor and win instead; pass an explicit HET code to disambiguate.
  - a 1-3 character PDB HET code (e.g. `"STI"`) — use that HETATM group in
    `structure`.
  - a path to an external ligand file (`.pdb`/`.sdf`/`.mol`/`.mol2`, e.g. a
    docking pose) with 3D coordinates already in `structure`'s frame, for
    structures that don't contain the ligand themselves.
- `outdir` — optional directory to keep fpocket's full raw output
  (`pockets/`, `*_info.txt`, ...); if `None`, a temporary directory is used
  and discarded.

The fpocket pocket whose lining atoms have the smallest minimum distance to
the ligand's atoms is selected. Returns a dict:

- `pocket_id` — fpocket's pocket number.
- `score` / `druggability_score` / `volume` — convenience fields pulled from
  fpocket's info file for the selected pocket.
- `residues` — list of `{"chain", "resnum", "icode", "resname"}` lining the
  pocket, in first-seen order, deduplicated. `icode` is `""` when the residue
  has no PDB insertion code — chymotrypsin-numbered serine proteases (e.g.
  trypsin, thrombin, factor Xa) commonly have insertion-code residues like
  `60A`/`60B` sharing a resnum with plain residue `60`, so `icode` is kept
  distinct from `resnum` rather than folded into it.
- `info` — the full raw fpocket score dict for the selected pocket (all
  fields from `*_info.txt`, values converted to `float` where possible).

### `protein.SOLVENT_AND_IONS`

A `frozenset` of common crystallization solvent/ion/additive HET codes (not
exhaustive) — the same exclusion list `find_pocket`'s ligand auto-detection
uses. Handy for deciding which HETATM groups in a structure are worth
displaying as a ligand (e.g. adding a `stick` style for everything in a
structure's HETATM records except this set, since cartoon-only styles don't
draw ligands at all).

## chem.ligand

### `ligand.list_ligand_instances(structure_path, exclude=SOLVENT_AND_IONS)`

Every non-excluded HETATM residue *instance* (physical occurrence) in a
structure file, as a list of `{"code", "chain", "resnum", "icode"}` dicts, in
file order. Two copies of the same ligand code (e.g. one per chain in a
dimer) produce two separate entries here, not one — use these dicts'
`chain`/`resnum`/`icode` to pin down a specific copy in `load_ligand`.
Solvent/ions (`chem.protein.SOLVENT_AND_IONS`) are excluded by default; pass
a wider `exclude` (e.g. adding a dataset-specific non-ligand HETATM code like
a glycosylation sugar) to filter those out too.

### `ligand.list_ligand_codes(structure_path, exclude=SOLVENT_AND_IONS)`

Every distinct non-excluded HETATM residue code (3-letter PDB chemical
component id) in a structure file, e.g. `["S54"]`. Multiple copies of the
same code collapse to a single entry here — use `list_ligand_instances` to
enumerate every physical occurrence instead.

### `ligand.load_ligand(structure_path, ligand, chain=None, resnum=None, icode=None)`

Extract a ligand from a structure file as a proper RDKit molecule.

- `structure_path` — path to a PDB file.
- `ligand` — its 3-letter PDB chemical component code (see
  `list_ligand_codes`/`list_ligand_instances`).
- `chain`/`resnum`/`icode` — pin down one specific instance when `ligand`'s
  code occurs more than once (values as found in a `list_ligand_instances`
  entry). Left as `None` (default), the most complete matching instance is
  used — fine when the code is known to be unique, or when any copy will do
  since they're chemically identical.

The residue's atoms and 3D coordinates come straight from the structure file,
but PDB format has no bond-order information, so RDKit's initial guess is all
single bonds with no aromaticity. This is corrected by fetching the PDB
Chemical Component Dictionary's ideal SMILES for `ligand` (trying a few of
its recorded descriptors, since some spell out a stereo-defining hydrogen as
an explicit atom that would otherwise break the atom-count match) and using
it as a bond-order template (`rdkit.Chem.AllChem.AssignBondOrdersFromTemplate`).

Raises `ValueError` if no matching residue is found (bad code, or bad
chain/resnum/icode), or if no candidate template matches the extracted atoms
— e.g. one piece of a covalently-linked multi-residue ligand (a
peptidomimetic inhibitor built from linked amino-acid HETATM groups
extracted on its own no longer has the same atoms as the free amino acid), or
a residue with incomplete crystallographic density.

### `ligand.qed(mol)`

Quantitative Estimate of Drug-likeness (Bickerton et al., 2012) for an RDKit
molecule, 0-1. Thin wrapper around `rdkit.Chem.QED.qed`.

### `ligand.molecular_weight(mol)`

Average molecular weight (g/mol) for an RDKit molecule. Thin wrapper around
`rdkit.Chem.Descriptors.MolWt`.

## chem.view3d

### `view3d.render_protein(path, exclude=SOLVENT_AND_IONS, width=600, height=500, coloring="spectrum", bfactor_range=(50, 90), style="cartoon")`

Display a PDB structure file as an interactive py3Dmol view in a light-gray
bordered frame, with a caption beside it on the right: the protein backbone
-- as a cartoon or a solid volume per `style`, colored per `coloring` -- plus
any HETATM ligand group not in `exclude` as magenta sticks (neither backbone
style draws ligands on its own). The border marks exactly the area where
3Dmol.js's mouse controls (rotate/zoom/pan) take over.

- `path` — path to a PDB file.
- `exclude` — HET codes to leave off the ligand sticks. Defaults to
  `chem.protein.SOLVENT_AND_IONS`; pass a superset (e.g.
  `SOLVENT_AND_IONS | {"NAG", "TYS"}`) to also exclude structure-specific
  non-ligand HETATM groups such as glycosylation sugars or modified residues.
- `width` / `height` — viewer size in pixels.
- `coloring` — `"spectrum"` (default): rainbow N -> C by residue position; or
  `"bfactor"`: rainbow by the file's per-atom B-factor column (e.g.
  AlphaFold's per-residue pLDDT confidence).
- `bfactor_range` — `(min, max)` the `"bfactor"` gradient is scaled over;
  ignored for `"spectrum"`. Defaults to AlphaFold's pLDDT confidence
  convention (`50, 90`); pass the structure's own B-factor range for
  crystallographic temperature factors.
- `style` — `"cartoon"` (default): ribbon backbone; or `"surface"`: a solid
  van der Waals volume (the union of a smooth blob at every backbone atom,
  computed by 3Dmol.js via marching cubes over a grid of atomic radii, rather
  than a ribbon or individual per-atom spheres), translucent (opacity 0.85)
  so a ligand bound underneath still shows through. Restricted to
  non-HETATM atoms either way, so the ligand keeps its own separate stick
  rendering instead of being enveloped by the volume.

The caption lists, one per line: the PDB id (`path`'s filename stem), the
distinct chain ids found in `ATOM` records, the ligand HET codes shown as
sticks (`"none"` if empty), and the experimental resolution parsed from the
file's legacy-PDB `REMARK 2 RESOLUTION` record (`"N/A"` if absent — e.g. NMR
structures, AlphaFold predictions, or files written by `chem.protein.align`,
which doesn't preserve header/REMARK records).

Displays the view and caption directly as a side effect and returns nothing
— just call it, no need to chain `.show()` or use it as a cell's last
expression.

## chem.ids (shared identifier resolution)

Internal to `chem.chembl`/`chem.rcsb`/`chem.alphafold`, not typically called
directly, but usable standalone:

- `ids.resolve_uniprot_accession(id_)` — UniProt accession or entry name →
  canonical UniProt accession.
- `ids.resolve_target_chembl_id(id_)` — ChEMBL target id / UniProt accession
  / UniProt entry name → ChEMBL target id.
- `ids.resolve_uniprot_accession_any(id_)` — ChEMBL target id / UniProt
  accession / UniProt entry name → UniProt accession (reverse-resolves via
  the ChEMBL API when given a ChEMBL target id).
