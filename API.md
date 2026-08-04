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

## chem.blast

### `blast.blastp(sequence, database, email="user@example.com", matrix="BLOSUM62", expect=1e-10, max_hits=50, poll_interval=10, timeout=600, title="chem.blast")`

Run a protein BLAST (`blastp`) search via
[EBI's Job Dispatcher REST API](https://www.ebi.ac.uk/jdispatcher/docs/webservices/)
and return the ranked hit list. Submits the search, polls until it finishes,
then fetches and parses the JSON result.

- `sequence` — a protein sequence, plain or FASTA-formatted (EBI accepts
  either).
- `database` — which EBI-hosted database to search, e.g. `"pdb"` (every PDB
  chain — hits are chain-level, e.g. `"1ABC_A"`, so a homolog with an
  experimentally solved structure appears once per chain) or
  `"uniprotkb_swissprot"` (reviewed UniProt entries, one hit per protein).
  See the webservices docs above for the full list EBI hosts.
- `email` — contact email required by EBI's Job Dispatcher API (their
  abuse-prevention/contact policy — not stored or used for anything else).
  Defaults to the placeholder `"user@example.com"`; pass your own if you're
  making many calls.
- `matrix` — substitution matrix, default `"BLOSUM62"`.
- `expect` — E-value threshold (upper bound, inclusive), default `1e-10` —
  notably stricter than EBI's own tool default (`10`). Not a free-form
  float: EBI's API only accepts one of a fixed set of values (`ValueError`
  otherwise) — `1e-200`, `1e-100`, `1e-50`, `1e-10`, `1e-5`, `1e-4`, `1e-3`,
  `1e-2`, `1e-1`, `1.0`, `10`, `100`, `1000`. Pick a larger one (e.g. `1e-3`,
  `1.0`) to surface more distant/lower-identity homologs.
- `max_hits` — maximum number of hits to request, default `50`. Also a fixed
  EBI enum (`ValueError` otherwise) — `0`, `5`, `10`, `20`, `50`, `100`,
  `150`, `200`, `250`, `500`, `750`, `1000`. EBI's own default is also `50`,
  but a search often has more hits than that within the `expect` cutoff —
  raise this to see further down the ranked list rather than just the
  closest matches.
- `poll_interval` — seconds between status polls while the job runs, default
  `10`.
- `timeout` — seconds to wait for the job to reach a terminal state before
  raising `TimeoutError`, default `600`.
- `title` — job title EBI records for the search, default `"chem.blast"`.

Returns a list of dicts, one per hit, in EBI's own ranking order (best
first): `accession`, `description`, `identity_pct`, `align_len`, `evalue`.
Raises `RuntimeError` if the search job itself ends in a failure state, and
`TimeoutError` if it doesn't reach a terminal state within `timeout` seconds.

**Not decorated with `chem.verbosity.logged`**, unlike every other function
in this reference (see the note at the top of this file) — that decorator
logs every bound argument, including `email`, to stderr on each call, which
would echo the caller's contact email for no benefit. `CHEM_QUIETNESS` still
silences its `tqdm` progress bar and the one-line hit-count summary printed
on success.

## chem.protein

### `protein.summary(id)`

Fetch a UniProt entry — by accession (e.g. `"Q8IZL9"`), entry
name/mnemonic (e.g. `"CDK20_HUMAN"`), or ChEMBL target id (e.g.
`"CHEMBL3559690"`), resolved via `chem.ids.resolve_uniprot_accession_any`
— and return a flat `dict` of properties useful for drug-discovery target
triage:

- `entry_name`, `accession` — UniProt entry name (mnemonic) and accession.
- `protein_name`, `gene_name` (with synonyms, if any), `organism`,
  `sequence_length`, `ec_number`.
- `family` — the `SIMILARITY` comment (protein family classification).
- `function` — the `FUNCTION` comment.
- `subcellular_location` — comma-joined `SUBCELLULAR LOCATION` values.
- `kinase_domain_range` — `"{start}-{end}"` residue range of the first
  `Domain` feature whose description mentions "kinase" (`None` for
  non-kinases or entries without a domain annotation).
- `active_site_residue` — residue number of the first `Active site`
  feature (`None` if unannotated).
- `n_pdb_xrefs` — count of `PDB` cross-references in the UniProt entry
  (not RCSB's own search API — use `chem.rcsb` for a live/authoritative
  count).
- `has_alphafold_model`, `has_bindingdb_entry` — whether an `AlphaFoldDB` /
  `BindingDB` cross-reference exists.
- `chembl_target_id` — the `ChEMBL` cross-reference id, if any.
- `pharos_development_level` — Pharos target development level plus a short
  gloss, e.g. `"Tbio (biology characterized, no known drug/chemical
  probe)"` (`Tclin`/`Tchem`/`Tbio`/`Tdark`; `None` if UniProt has no Pharos
  cross-reference).
- `protein_existence`, `annotation_score` — UniProt's own evidence-level and
  annotation-completeness scores.

Any field UniProt doesn't have data for comes back as `None` (or `False`
for the `has_*` flags) rather than raising.

### `protein.get_fasta(id, email="user@example.com")`

Fetch a UniProt entry's sequence as a FASTA string. `id` accepts the same
three forms as `protein.summary` — UniProt accession, entry name/mnemonic,
or ChEMBL target id — resolved the same way
(`chem.ids.resolve_uniprot_accession_any`).

`email` isn't required by UniProt's REST API, but is sent as a contact
address in the request's `User-Agent` header per [UniProt's own API usage
guidelines](https://www.uniprot.org/help/api); pass your own if you're
making many calls.

### `protein.sequence_align(accession, structures, canonical_feature_type=None, canonical_feature_description=None, marker_feature_types=None, chain=None)`

Align a set of PDB/CIF structures' observed sequences against a UniProt
canonical sequence, one residue position at a time — e.g. to compare WT and
mutant structures, or many ligand-bound states of the same receptor, on a
shared, canonical-numbered coordinate system. Unlike `protein.align`, this
doesn't superpose or write any structure files; it only produces aligned
sequence data, leaving rendering (text blocks, colored HTML, mutation lists,
...) to the caller.

- `accession` — UniProt accession (e.g. `"P07550"`).
- `structures` — list of PDB/CIF file paths (e.g. from
  `chem.rcsb.download_structures`).
- `canonical_feature_type` / `canonical_feature_description` — optional
  UniProt feature to slice canonical down to, e.g. `("Chain", "Beta-lactamase
  TEM")` to drop a cleaved signal peptide, or `("Domain", "Bromo 1")` for a
  single domain of a larger protein. `canonical_feature_type=None` (default)
  uses the full-length UniProt sequence as-is.
- `marker_feature_types` — optional tuple of UniProt feature types to collect
  as positions of interest, e.g. `("Active site",)` for an enzyme's catalytic
  residues, `("Binding site",)` for a receptor's ligand pocket, `("Site",)`
  for a non-catalytic functional residue. `None` (default) returns an empty
  marker set.
- `chain` — optional chain id to use in every structure, overriding
  auto-selection. Default: whichever chain in each structure has the
  highest-scoring global alignment to canonical — not simply "chain A", since
  chain lettering isn't consistent across depositions, and a structure can
  contain chains that aren't the target protein at all (a fusion partner
  spliced into a loop, e.g. T4 lysozyme in many GPCR structures; a
  stabilizing nanobody; other complex partners such as a heterotrimeric G
  protein).

For each structure, the selected chain's sequence is extracted from its ATOM
records via ProDy — not SEQRES, so this is exactly what's resolved in the
deposited coordinates, no more and no less — then globally aligned
(`Bio.Align.PairwiseAligner`) onto canonical's own numbering. Ties between
equally-optimal alignments (which can arise near a long unmodeled stretch
when a residue at one end coincidentally matches canonical at the other end)
are broken in favor of the alignment with the fewest separate matched blocks,
so a real gap comes out as one contiguous gap rather than a coincidental
cross-gap match plus a shorter gap elsewhere.

Returns a dict:
- `protein_name`, `organism` — from the UniProt entry.
- `canonical_seq` — the (possibly feature-sliced) canonical sequence.
- `feature_start`, `feature_end` — `canonical_seq`'s 1-based bounds within
  the full UniProt sequence (`1`, `len(full_seq)` when
  `canonical_feature_type` is `None`).
- `marker_positions` — set of 1-based `canonical_seq` positions collected
  from `marker_feature_types` (empty when `marker_feature_types` is `None`).
- `raw_sequences` — `{path: sequence}`, the selected chain's CA-derived
  sequence exactly as observed, before alignment.
- `sequences` — `{path: sequence}`, each raw sequence collapsed onto
  canonical's coordinates — for every `path`, `len(sequences[path]) ==
  len(canonical_seq)`, with `-` marking any canonical position not observed
  in that structure (missing density, or outside the modeled construct; a
  residue belonging to a fusion partner or other foreign chain is dropped
  entirely rather than shown, since it has no canonical position of its own).

### `protein.align(structures, reference=None, chain=None, outdir="aligned")`

Sequence-align and structurally superpose a set of same-target structures
(PDB/CIF freely mixed — combine `chem.rcsb`/`chem.alphafold` downloads).

- `structures` — list of PDB/CIF file paths for the same protein.
- `reference` — which structure to align everything onto: an index into
  `structures`, or a path. Defaults to `structures[0]`. Does not need to be a
  member of `structures`; either way it's written to `outdir` exactly once.
- `chain` — optional chain id to use in every non-reference structure,
  overriding auto-selection. Default: whichever chain has the most residues
  identical to the reference at matched (gap-free) sequence positions — not
  simply the chain with the most residues overall, since a structure can
  contain a larger bound partner protein (e.g. thrombin co-crystallized with
  a serpin inhibitor next to its own, smaller heavy chain) that a size-only
  heuristic would wrongly prefer; gap-averse global alignment will also
  happily align most of an unrelated chain's length at near-zero identity, so
  raw matched-position count doesn't discriminate either. The reference's own
  chain is still picked by size alone (its primary polymer chain), since
  there's nothing yet to compare it against — pass an unambiguous
  single-chain reference (e.g. an AlphaFold prediction) if in doubt. HETATM
  residues are excluded even when their resname matches a standard amino
  acid, so a covalently-linked peptidomimetic ligand sharing the protein's
  chain id (e.g. a D-amino-acid-containing inhibitor) doesn't get pulled into
  the sequence.
- `outdir` — destination directory; created if missing.

For each non-reference structure, its selected chain is sequence-aligned
against the reference chain (`Bio.Align.PairwiseAligner`, global mode), and
the whole structure is superposed via a Kabsch fit (`Bio.PDB.Superimposer`)
on the sequence-matched CA atoms, applied to every atom including ligands and
waters. Every structure, including the reference, is written out as a PDB
file in `outdir` (regardless of input format) — ready to be loaded and
overlaid, e.g. one py3Dmol `addModel` call per file, or as input to
`chem.protein.find_pocket`.

Returns `{path: {"rmsd": ..., "identity": ...}}` over the sequence-matched CA
atoms for every structure that could be aligned (the reference maps to
`{"rmsd": 0.0, "identity": 1.0}`). `rmsd` is in Ångströms; `identity` is the
fraction of matched (gap-free) positions with an identical residue — a
mismatch (substitution) still counts as "matched" as long as no gap was
opened there, so `identity` can be `< 1.0` even when every residue found a
counterpart; gapped positions (e.g. a loop present in one structure but not
the other) count in neither the numerator nor the denominator. Both are plain
floats rounded to 3 decimal places. A structure with no usable chain, or too
few residues in common with the reference, is skipped with a warning (unless
quiet) rather than raising.

Large structures (e.g. cryo-EM assemblies with >26 chains or >99999 atoms)
are not supported by the legacy PDB writer used here.

### `protein.compute_transform(mobile, reference, mobile_chain=None, reference_chain=None)`

Compute the rigid-body superposition of `mobile`'s chain onto `reference`'s
chain, without moving or writing anything — the same underlying alignment
`align` uses, but handed back as raw numbers instead of applied to a whole
structure file. Lets you compute the transform from one chain of a
multi-chain complex (e.g. a kinase chain written out by `chem.protein.split`)
and reuse the exact same rotation/translation, via `apply_transform`, on that
entry's other chains (a bound partner) and ligand SDFs
(`chem.ligand.apply_transform`) — reconstructing a whole complex in a common
frame, including pieces (partner chains, ligands) that have no sequence of
their own to align on.

- `mobile`, `reference` — PDB/CIF file paths.
- `mobile_chain`, `reference_chain` — optional chain ids to use (default: for
  `reference`, its primary polymer chain; for `mobile`, whichever chain
  actually matches the reference best — same selection `align` uses for its
  reference/`chain` handling).

Returns `{"rotation": 3x3 nested list, "translation": 3-element list, "rmsd":
float, "identity": float}`, in Biopython's `Superimposer`/`Atom.transform`
convention (`new_coord = old_coord @ rotation + translation`, row-vector
coordinates) — pass `rotation`/`translation` straight to `apply_transform`.
`rmsd` is in Ångströms over the sequence-matched CA atoms; `identity` is the
same definition `align` reports. Raises `ValueError` if either structure has
too few residues, or too few of them match the other, to superpose on.

### `protein.apply_transform(structure_path, rotation, translation, outpath)`

Apply a rotation+translation — as returned by `compute_transform` — to every
atom of a PDB/CIF structure file, and write the result as PDB.

- `structure_path` — PDB/CIF file to transform (e.g. a `chem.protein.split`
  chain PDB that wasn't itself used to compute the transform, such as a bound
  partner chain).
- `rotation`, `translation` — as returned by `compute_transform`'s
  `"rotation"`/`"translation"`.
- `outpath` — destination PDB file path; parent directory created if
  missing.

Returns `outpath`.

### `protein.identity_matrix(structures, chain=None)`

Pairwise sequence identity across a set of structures — e.g. the per-chain
protein PDB files `chem.protein.split` writes out, to see at a glance which
chains are the same protein (~1.0) versus unrelated (near 0). Same identity
definition `align` reports, but computed for every pair among `structures`
directly, with no reference or structural superposition step.

- `structures` — list of PDB/CIF file paths.
- `chain` — optional chain id to use in every structure (default: each
  structure's own primary polymer chain — the one with the most standard
  amino acid residues, same auto-selection `align` uses for its reference).
  Most `chem.protein.split` outputs already contain a single chain, so this
  rarely needs to be passed.

Returns a dict of dicts, `{path_i: {path_j: identity, ...}, ...}`, one entry
for every pair among the structures that had a usable chain (including
`path_i == path_i` → `1.0`), symmetric (`identity[a][b] == identity[b][a]`).
`identity` is a plain float rounded to 3 decimal places, same definition as
`align`'s. A structure with no chain usable for sequence comparison (no
polymer residues, or the requested `chain` id not found) is skipped with a
warning (unless quiet) and omitted from the matrix entirely, rather than
raising.

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
- `spheres` — list of `{"x", "y", "z", "radius"}`: fpocket's alpha spheres
  approximating the pocket cavity's shape/volume (as opposed to `residues`'
  lining protein atoms), in the same frame as `structure`. Handy for a
  "filled space" view, e.g. one py3Dmol `addSphere` per entry, rather than
  sticks on the lining residues.
- `info` — the full raw fpocket score dict for the selected pocket (all
  fields from `*_info.txt`, values converted to `float` where possible).

### `protein.list_pockets(structure, outdir=None, druggability_thres=0.1)`

Run fpocket on a PDB file and return every pocket it detects — no reference
ligand needed, unlike `find_pocket`. Useful for blind pocket detection on
structures with no bound ligand at all, e.g. an AlphaFold prediction, to see
every candidate binding site rather than just the one nearest a known ligand.

- `structure` — path to a PDB file.
- `outdir` — same as `find_pocket`'s.
- `druggability_thres` — minimum `druggability_score` (inclusive) to keep.
  fpocket routinely reports dozens of low-quality cavities with near-zero
  scores on a typical structure; the default (`0.1`) drops those, along with
  any pocket fpocket didn't assign a `druggability_score` to at all. Pass
  `None` to keep every detected pocket, unfiltered.

Returns a list of dicts, one per kept pocket, each shaped exactly like a
single `find_pocket` result (`pocket_id`/`score`/`druggability_score`/
`volume`/`residues`/`spheres`/`info`), sorted by `druggability_score`
descending.

### `protein.residues_near(structure_path, ligands, radius=8.0)`

Protein residues in `structure_path` with at least one atom within `radius`
Å of any atom in `ligands` — a simple distance-based active-site definition,
unlike `find_pocket`/`list_pockets` (which run fpocket's cavity detection).
Handy when the ligand(s) already have their own file(s) in the same
coordinate frame as `structure_path` — e.g. the ligand-free protein PDB and
per-ligand SDF(s) `chem.protein.split` writes out, optionally after
`apply_transform`/`chem.ligand.apply_transform` moved both into a shared
reference frame.

- `structure_path` — PDB/CIF file. HETATM residues in it (water, ligands,
  ions — if any are still present) are never returned, even if within
  `radius`: this only ever selects polymer (ATOM) residues.
- `ligands` — one ligand file path, or a list of them (`.pdb`/`.sdf`/`.mol`/
  `.mol2`, 3D coordinates already in the same frame as `structure_path`).
- `radius` — distance threshold in Å, inclusive. Default `8.0`.

Returns a list of `{"chain", "resnum", "icode", "resname"}` (`icode` is `""`
when the residue has no PDB insertion code), sorted by `(chain, resnum,
icode)` — same shape as `find_pocket`'s `"residues"`.

### `protein.split(structure_path, all_chains=False, remove_water=False, outdir="split")`

Split a structure file into a ligand-free protein PDB and one SDF file per
non-water HETATM ligand instance — e.g. to prep a receptor/ligand pair for
docking.

- `structure_path` — path to a PDB/CIF structure file.
- `all_chains` — if `True`, write a single protein PDB with every chain
  together instead of one PDB per chain. Default `False` (split by chain).
- `remove_water` — if `True`, also strip water (HETATM `HOH`/`WAT`/`DOD`) out
  of the protein PDB. Default `False` (water is kept).
- `outdir` — destination directory; created if missing. Default `"split"`.

Returns a dict:

- `"protein"` — a `{chain_id: path}` dict, one entry per chain
  (`all_chains=False`, the default), or the single protein PDB path
  (`all_chains=True`). By default water is kept — pass `remove_water=True` to
  strip it out too. Every other HETATM residue — real ligands, ions,
  crystallization additives, glycosylation sugars, alike — is always
  stripped, since those are exactly what end up in `"ligands"` below instead.
- `"ligands"` — list of `{"path", "code", "chain", "resnum", "icode",
  "bond_orders_restored"}`, one entry per non-water HETATM residue *instance*
  (see `ligand.list_ligand_instances`) — e.g. two copies of the same ligand
  code bound to different chains produce two entries, each its own SDF file.
  Every such instance gets an SDF, no matter what — 3D coordinates always come
  straight from `structure_path`. When `ligand.load_ligand` can resolve a
  bond-order template for it against the PDB Chemical Component Dictionary,
  that's what's written (`bond_orders_restored=True`: proper bond
  orders/aromaticity). When it can't (e.g. a covalently-linked glycosylation
  sugar or peptidomimetic ligand missing the atom(s) involved in that link,
  relative to the free/standalone template; or incomplete crystallographic
  density), a warning is printed (unless quiet) and the raw connectivity
  RDKit guesses from atomic distances is written instead
  (`bond_orders_restored=False`: single bonds only, no aromaticity/stereo).

### `protein.SOLVENT_AND_IONS`

A `frozenset` of common crystallization solvent/ion/additive HET codes (not
exhaustive) — the same exclusion list `find_pocket`'s ligand auto-detection
and `chem.ligand`'s default `exclude` use, so bare ions/crystallization junk
don't get treated as "the" ligand. For deciding what's worth *displaying*,
see `WATER` below instead — an ion or additive is usually still worth
seeing, even though it isn't drug-like.

### `protein.WATER`

A `frozenset` of just water HET codes (`HOH`, `WAT`, `DOD`) — a strict subset
of `SOLVENT_AND_IONS`. This is `view3d.render_protein`'s default `exclude`:
everything except water is drawn as a ligand, since bound ions and
crystallization additives are usually worth seeing even though they aren't
drug-like.

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

### `ligand.apply_transform(sdf_path, rotation, translation, outpath)`

Apply a rotation+translation — as returned by `chem.protein.compute_transform`
— to an SDF file's 3D coordinates and write the result to `outpath`.
Complements `chem.protein.apply_transform` (for protein PDB chains): a
transform computed from one chain of a complex (e.g. via
`chem.protein.compute_transform` on the kinase chain that `chem.protein.split`
wrote out) can be reapplied here to that same entry's ligand SDF(s), so the
ligand ends up in the same coordinate frame as the aligned protein —
reconstituting the bound complex without needing a sequence to align the
ligand on.

- `sdf_path` — SDF file to transform (e.g. a `chem.protein.split` ligand
  path). Read with `sanitize=False`, since some SDFs written by `split`
  (`bond_orders_restored=False`) can't survive sanitization.
- `rotation`, `translation` — as returned by `chem.protein.compute_transform`'s
  `"rotation"`/`"translation"` — same convention as `chem.protein.apply_transform`,
  so the same pair of matrices moves both the protein and the ligand
  consistently.
- `outpath` — destination SDF file path; parent directory created if
  missing.

Returns `outpath`.

### `ligand.qed(mol)`

Quantitative Estimate of Drug-likeness (Bickerton et al., 2012) for an RDKit
molecule, 0-1. Thin wrapper around `rdkit.Chem.QED.qed`.

### `ligand.molecular_weight(mol)`

Average molecular weight (g/mol) for an RDKit molecule. Thin wrapper around
`rdkit.Chem.Descriptors.MolWt`.

## chem.view3d

### `view3d.render_protein(path, exclude=WATER, width=600, height=500, coloring="spectrum", bfactor_range=(50, 90), style="cartoon")`

Display a PDB structure file as an interactive py3Dmol view in a light-gray
bordered frame, with a caption beside it on the right: the protein backbone
-- as a cartoon or a solid volume per `style`, colored per `coloring` -- plus
any HETATM ligand group not in `exclude` as magenta sticks (neither backbone
style draws ligands on its own). The border marks exactly the area where
3Dmol.js's mouse controls (rotate/zoom/pan) take over.

- `path` — path to a PDB file.
- `exclude` — HET codes to leave off the ligand sticks. Defaults to
  `chem.protein.WATER` (just water) — bound ions and crystallization
  additives are shown, since they're usually worth seeing even though they
  aren't drug-like. Pass `chem.protein.SOLVENT_AND_IONS` for the old broader
  default, or a superset (e.g. `SOLVENT_AND_IONS | {"NAG", "TYS"}`) to also
  exclude structure-specific non-ligand HETATM groups such as glycosylation
  sugars or modified residues.
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
