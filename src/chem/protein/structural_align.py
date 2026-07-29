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
    return max(chains, key=lambda c: sum(1 for r in c if is_aa(r, standard=True)))


def _chain_seq_and_ca(chain):
    residues = [r for r in chain if is_aa(r, standard=True) and "CA" in r]
    seq = "".join(index_to_one(three_to_index(r.get_resname())) for r in residues)
    ca_atoms = [r["CA"] for r in residues]
    return seq, ca_atoms


def _matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca):
    """Sequence-align ref_seq/mob_seq and return the CA atoms at every matched
    (non-gap) position, as parallel (ref_points, mob_points) lists.
    """
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2
    aligner.mismatch_score = -1
    alignment = aligner.align(ref_seq, mob_seq)[0]

    ref_pts, mob_pts = [], []
    for (r0, r1), (m0, m1) in zip(*alignment.aligned):
        for k in range(r1 - r0):
            ref_pts.append(ref_ca[r0 + k])
            mob_pts.append(mob_ca[m0 + k])
    return ref_pts, mob_pts


@logged
def align(structures, reference=None, chain=None, outdir="aligned"):
    """Sequence-align and structurally superpose a set of same-target structures.

    For each non-reference structure: selects its primary polymer chain (or `chain`
    if given), sequence-aligns it against the reference chain, and superposes the
    whole structure -- Kabsch fit on the sequence-matched CA atoms, applied to every
    atom including ligands and waters -- onto the reference's coordinate frame.
    Every structure, including the reference itself, is written out as a PDB file in
    `outdir` (regardless of input format), ready to be loaded and overlaid, e.g. one
    py3Dmol addModel call per file, or as input to chem.protein.find_pocket.

    structures: list of PDB/CIF file paths for the same protein (e.g. downloaded via
        chem.rcsb/chem.alphafold).
    reference: which structure to align everything onto -- an index into
        `structures`, or one of its paths. Defaults to structures[0].
    chain: optional chain id to use in every structure, overriding auto-selection.
        Default: auto-select each structure's chain with the most standard amino
        acid residues (its primary polymer chain) -- e.g. thrombin's catalytic
        heavy chain rather than its short light chain.
    outdir: destination directory; created if missing.

    Returns {path: rmsd} over the sequence-matched CA atoms for every structure that
    could be aligned (the reference maps to 0.0). A structure with no usable chain,
    or too few residues in common with the reference, is skipped with a warning
    (unless quiet) rather than raising.

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

    for path in tqdm(
        structures, desc=f"aligning structures to {ref_path}", unit="structure", disable=is_quiet()
    ):
        out_path = os.path.join(outdir, os.path.splitext(os.path.basename(path))[0] + ".pdb")

        if path == ref_path:
            io.set_structure(ref_structure)
            io.save(out_path)
            results[path] = 0.0
            continue

        try:
            mob_structure = _load_structure(path)
            mob_chain = _select_chain(next(mob_structure.get_models()), chain)
            mob_seq, mob_ca = _chain_seq_and_ca(mob_chain)
            ref_pts, mob_pts = _matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca)
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
        results[path] = sup.rms

    return results
