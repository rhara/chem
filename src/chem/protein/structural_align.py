import os
import sys

from Bio.Align import PairwiseAligner
from Bio.PDB import MMCIFParser, PDBIO, PDBParser
from Bio.PDB.Polypeptide import index_to_one, is_aa, three_to_index
from Bio.PDB.Superimposer import Superimposer
from tqdm import tqdm

from ..verbosity import is_quiet, logged

# The minimum number of matched CA atoms a rigid-body superposition is well-posed
# with; in practice legitimate same-target alignments match far more than this.
_MIN_MATCHED_RESIDUES = 3


def _load_structure(path):
    ext = os.path.splitext(path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    return parser.get_structure(os.path.splitext(os.path.basename(path))[0], path)


def _is_polymer_residue(r):
    """A real polymer (ATOM) residue -- excludes HETATM residues even when their
    resname matches a standard amino acid, e.g. a D-amino acid or proline in a
    covalently-linked peptidomimetic ligand that shares the protein's chain id.
    Bio.PDB's is_aa() checks resname only, not the hetero flag, so it alone would
    wrongly pull such ligand residues into the sequence/coordinate list.
    """
    return r.id[0] == " " and is_aa(r, standard=True)


def _select_chain(model, chain_id=None):
    """Return chain_id if given, otherwise the model's chain with the most standard
    amino acid residues (its primary polymer chain).
    """
    chains = list(model)
    if chain_id is not None:
        matches = [c for c in chains if c.id == chain_id]
        if not matches:
            raise ValueError(f"chain '{chain_id}' not found")
        return matches[0]
    return max(chains, key=lambda c: sum(1 for r in c if _is_polymer_residue(r)))


def _chain_seq_and_ca(chain):
    residues = [r for r in chain if _is_polymer_residue(r) and "CA" in r]
    seq = "".join(index_to_one(three_to_index(r.get_resname())) for r in residues)
    ca_atoms = [r["CA"] for r in residues]
    return seq, ca_atoms


def _matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca):
    """Sequence-align ref_seq/mob_seq and return the CA atoms at every matched
    (non-gap) position, as parallel (ref_points, mob_points) lists, plus the
    fraction of those matched positions where the residue is identical
    (gapped/unmatched positions -- residues with nothing on the other side --
    aren't counted in either the numerator or the denominator).
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2
    aligner.mismatch_score = -1
    alignment = aligner.align(ref_seq, mob_seq)[0]

    ref_pts, mob_pts = [], []
    n_identical = 0
    for (r0, r1), (m0, m1) in zip(*alignment.aligned):
        for k in range(r1 - r0):
            ref_pts.append(ref_ca[r0 + k])
            mob_pts.append(mob_ca[m0 + k])
            if ref_seq[r0 + k] == mob_seq[m0 + k]:
                n_identical += 1
    identity = n_identical / len(ref_pts) if ref_pts else 0.0
    return ref_pts, mob_pts, identity


def _best_matching_chain_alignment(model, ref_seq, ref_ca):
    """Align every polymer chain in model against (ref_seq, ref_ca) and return
    the (ref_pts, mob_pts, identity) for whichever has the most identical
    matched residues (identity * matched count).

    A plain "most residues" chain-size heuristic picks the wrong chain when the
    structure also contains a larger bound partner protein -- e.g. thrombin's
    real structure files sometimes include a co-crystallized serpin inhibitor
    (a full ~350-400 residue fold) alongside thrombin's own ~250-residue heavy
    chain, and gap-averse global alignment (open_gap_score=-10) will happily
    "match" most of an unrelated chain's length at near-zero identity rather
    than open gaps, so matched-count alone doesn't discriminate either --
    identity is what actually tells the same-protein chain apart.
    """
    best = None  # (n_identical, ref_pts, mob_pts, identity)
    for chain in model:
        seq, ca = _chain_seq_and_ca(chain)
        if not seq:
            continue
        ref_pts, mob_pts, identity = _matched_ca_pairs(ref_seq, ref_ca, seq, ca)
        n_identical = identity * len(ref_pts)
        if best is None or n_identical > best[0]:
            best = (n_identical, ref_pts, mob_pts, identity)
    if best is None:
        raise ValueError("no chain with polymer residues found")
    return best[1], best[2], best[3]


@logged
def align(structures, reference=None, chain=None, outdir="aligned"):
    """Sequence-align and structurally superpose a set of same-target structures.

    For each non-reference structure: selects its chain -- `chain` if given,
    otherwise whichever chain actually matches the reference best (see below) --
    sequence-aligns it against the reference chain, and superposes the whole
    structure -- Kabsch fit on the sequence-matched CA atoms, applied to every
    atom including ligands and waters -- onto the reference's coordinate frame.
    Every structure, including the reference itself, is written out as a PDB file in
    `outdir` (regardless of input format), ready to be loaded and overlaid, e.g. one
    py3Dmol addModel call per file, or as input to chem.protein.find_pocket.

    structures: list of PDB/CIF file paths for the same protein (e.g. downloaded via
        chem.rcsb/chem.alphafold).
    reference: which structure to align everything onto -- an index into
        `structures`, or a path. Defaults to structures[0]. The reference does not
        need to be a member of `structures`; either way it is written to `outdir`
        exactly once (it is skipped in the alignment loop if it also appears in
        `structures`). The reference's own chain is picked by size alone (its
        primary polymer chain, the one with the most standard amino acid
        residues) since there's nothing yet to compare it against -- pass an
        unambiguous single-chain reference (e.g. an AlphaFold prediction) if in
        doubt.
    chain: optional chain id to use in every non-reference structure, overriding
        auto-selection. Default: pick whichever chain has the most residues
        identical to the reference at matched (gap-free) sequence positions --
        not simply the chain with the most residues overall, since a structure
        can contain a larger bound partner protein (e.g. thrombin co-crystallized
        with a ~350-400-residue serpin inhibitor, next to its own ~250-residue
        heavy chain) that a size-only heuristic would wrongly prefer. Gap-averse
        global alignment (a single substitution is far cheaper than opening a
        gap) will also happily align most of an unrelated chain's length at
        near-zero identity rather than open gaps, so raw matched-position count
        doesn't discriminate either -- identity is what actually distinguishes
        the same-protein chain from an unrelated bound partner.
    outdir: destination directory; created if missing.

    Returns {path: {"rmsd": ..., "identity": ...}} over the sequence-matched CA
    atoms for every structure that could be aligned (the reference maps to
    {"rmsd": 0.0, "identity": 1.0}). rmsd is in Angstroms; identity is the
    fraction of matched (gap-free) positions with an identical residue -- gapped
    positions (residues with nothing on the other side, e.g. a loop present in
    one structure but not the other) don't count in either the numerator or the
    denominator. Both are plain floats rounded to 3 decimal places. A structure
    with no usable chain, or too few residues in common with the reference, is
    skipped with a warning (unless quiet) rather than raising.

    Very large structures (e.g. cryo-EM assemblies with >26 chains or >99999 atoms)
    are not supported by the legacy PDB writer used here.
    """
    if not structures:
        raise ValueError("structures must be non-empty")

    if reference is None:
        ref_path = structures[0]
    elif isinstance(reference, int):
        ref_path = structures[reference]
    else:
        ref_path = reference

    ref_structure = _load_structure(ref_path)
    ref_chain = _select_chain(next(ref_structure.get_models()), chain)
    ref_seq, ref_ca = _chain_seq_and_ca(ref_chain)
    if len(ref_ca) < _MIN_MATCHED_RESIDUES:
        raise ValueError(f"reference structure {ref_path} has too few residues to align on")

    os.makedirs(outdir, exist_ok=True)
    io = PDBIO()
    results = {}

    ref_out_path = os.path.join(outdir, os.path.splitext(os.path.basename(ref_path))[0] + ".pdb")
    io.set_structure(ref_structure)
    io.save(ref_out_path)
    results[ref_path] = {"rmsd": 0.0, "identity": 1.0}

    for path in tqdm(
        structures, desc=f"aligning structures to {ref_path}", unit="structure", disable=is_quiet()
    ):
        if path == ref_path:
            continue

        out_path = os.path.join(outdir, os.path.splitext(os.path.basename(path))[0] + ".pdb")

        try:
            mob_structure = _load_structure(path)
            mob_model = next(mob_structure.get_models())
            if chain is not None:
                mob_chain = _select_chain(mob_model, chain)
                mob_seq, mob_ca = _chain_seq_and_ca(mob_chain)
                ref_pts, mob_pts, identity = _matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca)
            else:
                # Reference-aware: picks whichever chain is actually the same
                # protein as the reference, not just the largest chain in the
                # file (see _best_matching_chain_alignment).
                ref_pts, mob_pts, identity = _best_matching_chain_alignment(mob_model, ref_seq, ref_ca)
            if len(ref_pts) < _MIN_MATCHED_RESIDUES:
                raise ValueError(f"only {len(ref_pts)} residue(s) matched the reference sequence")
        except ValueError as e:
            if not is_quiet():
                print(f"skipping {path}: {e}", file=sys.stderr)
            continue

        sup = Superimposer()
        sup.set_atoms(ref_pts, mob_pts)
        sup.apply(list(mob_structure.get_atoms()))

        io.set_structure(mob_structure)
        io.save(out_path)
        results[path] = {"rmsd": round(float(sup.rms), 3), "identity": round(identity, 3)}

    return results


@logged
def identity_matrix(structures, chain=None):
    """Pairwise sequence identity matrix across a set of structures -- e.g. the
    per-chain protein PDB files chem.protein.split writes out, to see at a glance
    which chains are the same protein (~1.0) versus unrelated (near 0).

    structures: list of PDB/CIF file paths.
    chain: optional chain id to use in every structure (default: each structure's
        own primary polymer chain, the one with the most standard amino acid
        residues -- same auto-selection align() uses for its reference). Most
        chem.protein.split outputs already contain a single chain, so this
        rarely needs to be passed.

    Returns a dict of dicts, `{path_i: {path_j: identity, ...}, ...}`, one entry
    for every pair among the structures that had a usable chain (including
    path_i == path_i -> 1.0), symmetric (identity[a][b] == identity[b][a]).
    identity is the fraction of matched (gap-free) sequence positions with an
    identical residue -- same definition as align()'s identity, a plain float
    rounded to 3 decimal places. A structure with no chain usable for sequence
    comparison (no polymer residues, or the requested chain id not found) is
    skipped with a warning (unless quiet) and omitted from the matrix entirely,
    rather than raising.
    """
    seqs = {}
    for path in structures:
        try:
            structure = _load_structure(path)
            sel_chain = _select_chain(next(structure.get_models()), chain)
            seq, ca = _chain_seq_and_ca(sel_chain)
            if not seq:
                raise ValueError("no standard amino acid residues with a CA atom found")
        except ValueError as e:
            if not is_quiet():
                print(f"skipping {path}: {e}", file=sys.stderr)
            continue
        seqs[path] = (seq, ca)

    paths = list(seqs)
    matrix = {p: {p: 1.0} for p in paths}
    for i, path_i in enumerate(paths):
        seq_i, ca_i = seqs[path_i]
        for path_j in paths[i + 1 :]:
            seq_j, ca_j = seqs[path_j]
            _, _, identity = _matched_ca_pairs(seq_i, ca_i, seq_j, ca_j)
            identity = round(identity, 3)
            matrix[path_i][path_j] = identity
            matrix[path_j][path_i] = identity

    return matrix
