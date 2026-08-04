"""Compare a set of PDB/CIF structures' observed sequences against a UniProt
canonical sequence, one residue position at a time. Rendering (text blocks,
colored HTML, mutation summaries, ...) is deliberately left to the caller --
this module only produces the aligned data.
"""

import itertools
import os

import prody
import requests
from Bio.Align import PairwiseAligner
from tqdm import tqdm

from ..verbosity import is_quiet, logged

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"

# Cap on how many equally-optimal alignments to compare when breaking ties (see
# _align_to_canonical) -- large enough for the ties seen in practice (usually a
# handful, at most a couple dozen, from short runs of identical residues near a
# large gap), small enough to never be slow even for a pathologically
# repetitive sequence.
_MAX_TIED_ALIGNMENTS = 500


def _make_aligner():
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2
    # Mismatch is penalized more harshly than a plain -1 would: with a mild
    # mismatch penalty, a long foreign insertion (e.g. T4 lysozyme fused into a
    # GPCR's loop, ~160 residues) can come out *cheaper* as a smear of many
    # small mismatches at the junction than as one clean gap, since a single
    # gap-open is a steep one-time cost. -3 makes each mismatch expensive
    # enough that this is consistently represented as a clean gap instead --
    # for the large insertions actually seen in practice (dozens of residues
    # or more). Going harsher still (e.g. matching open_gap_score) overcorrects
    # a different way: extending an *already-open* adjacent gap only costs
    # extend_gap_score (-0.5), so any mismatch penalty above that will tempt a
    # real point mutation sitting right next to a real gap to get silently
    # absorbed into the gap instead of reported -- observed in practice as
    # e.g. a genuine S262D or C265F disappearing from the mutation list. -3
    # is harsh enough to fix the large-insertion smearing actually seen in
    # real structures while staying gentle enough not to swallow adjacent
    # point mutations; it can still smear a couple of residues in a
    # small-scale, low-stakes edge case (a short foreign block only a residue
    # or two longer than the canonical span it replaces), which is an
    # accepted, irreducible tradeoff of a single affine-gap scoring scheme
    # covering both regimes.
    aligner.mismatch_score = -3
    return aligner


def _fetch_canonical(accession, feature_type, feature_description):
    resp = requests.get(f"{UNIPROT_API}/{accession}.json", timeout=30)
    resp.raise_for_status()
    entry = resp.json()
    full_seq = entry["sequence"]["value"]

    if feature_type is None:
        start, end = 1, len(full_seq)
    else:
        feature = next(
            f
            for f in entry["features"]
            if f["type"] == feature_type
            and (feature_description is None or f["description"] == feature_description)
        )
        start = feature["location"]["start"]["value"]
        end = feature["location"]["end"]["value"]

    return entry, full_seq[start - 1 : end], start, end


def _marker_positions(entry, start, end, marker_feature_types):
    if not marker_feature_types:
        return set()
    return {
        f["location"]["start"]["value"] - start + 1
        for f in entry["features"]
        if f["type"] in marker_feature_types and start <= f["location"]["start"]["value"] <= end
    }


def _load_ca_sequences(path):
    """Return {chain_id: CA-derived one-letter sequence} for every polymer chain
    in a PDB or mmCIF file -- from the ATOM records actually present, not
    SEQRES, so this reflects exactly what's resolved in the deposited
    coordinates, no more and no less.
    """
    ext = os.path.splitext(path)[1].lower()
    structure = prody.parseMMCIF(path) if ext in (".cif", ".mmcif") else prody.parsePDB(path)
    sequences = {}
    for chain in prody.HierView(structure):
        ca_atoms = chain.select("protein and name CA")
        if ca_atoms is not None:
            sequences[chain.getChid()] = ca_atoms.getSequence()
    return sequences


def _best_matching_chain_sequence(aligner, canonical_seq, path, chain):
    """Return the CA-derived sequence of whichever polymer chain in path's
    structure best matches canonical_seq (highest global-alignment score), or
    of `chain` specifically if given.

    Not simply "chain A": chain lettering isn't consistent across depositions,
    and a structure can contain chains that aren't the target protein at all --
    a fusion partner spliced into a loop (e.g. T4 lysozyme in many GPCR
    structures), a stabilizing nanobody, or other complex partners (e.g. a
    heterotrimeric G protein).
    """
    chain_sequences = _load_ca_sequences(path)
    if chain is not None:
        if chain not in chain_sequences:
            raise ValueError(f"chain '{chain}' not found in {path}")
        return chain_sequences[chain]

    best_seq, best_score = None, None
    for seq in chain_sequences.values():
        score = aligner.align(canonical_seq, seq).score
        if best_score is None or score > best_score:
            best_seq, best_score = seq, score
    if best_seq is None:
        raise ValueError(f"no polymer chain with CA atoms found in {path}")
    return best_seq


def _align_to_canonical(aligner, canonical_seq, query_seq):
    """Align query_seq to canonical_seq and collapse it onto canonical's own
    coordinates: returns a string of exactly len(canonical_seq) characters, one
    per canonical position, holding query_seq's residue there or '-' if that
    canonical position isn't observed in query_seq (missing density, or a
    residue outside the modeled construct). Any query residue that aligns to a
    gap *in canonical* (e.g. a fusion partner such as T4 lysozyme spliced into
    a loop for crystallization, or an expression tag) is dropped -- it doesn't
    correspond to any canonical position, so there's no single shared column to
    place it in when comparing many structures that carry different
    insertions of different lengths.

    When a residue near one end of query_seq happens to also match canonical
    at the *other* end of a long unmodeled stretch (e.g. a Leu right after the
    last resolved residue, with canonical's very last residue also a Leu, and
    70+ residues of unresolved tail in between), a plain global alignment can
    score that spurious cross-gap match exactly the same as correctly placing
    it next to the rest of the matched block -- Biopython then returns
    whichever of the tied-best alignments it happens to find first, which is
    not necessarily the biologically real one. Comparing the tied-best
    alignments (up to _MAX_TIED_ALIGNMENTS of them) and keeping the one with
    the fewest separate matched blocks resolves this in favor of one
    contiguous gap over "gap, one coincidental match, another gap" -- though a
    tie can still remain when canonical itself ends (or starts) in a run of
    2+ identical residues, which no scoring scheme can resolve from sequence
    alone.
    """
    alignments = aligner.align(canonical_seq, query_seq)
    best_n_blocks, best_alignment = None, None
    for alignment in itertools.islice(alignments, _MAX_TIED_ALIGNMENTS):
        n_blocks = len(alignment.aligned[0])
        if best_n_blocks is None or n_blocks < best_n_blocks:
            best_n_blocks, best_alignment = n_blocks, alignment
    gapped_canonical, gapped_query = str(best_alignment[0]), str(best_alignment[1])
    return "".join(q for c, q in zip(gapped_canonical, gapped_query) if c != "-")


@logged
def sequence_align(
    accession,
    structures,
    canonical_feature_type=None,
    canonical_feature_description=None,
    marker_feature_types=None,
    chain=None,
):
    """Align a set of PDB/CIF structures' observed sequences against a UniProt
    canonical sequence, one residue position at a time.

    For each structure, the chain that best matches canonical (or `chain` if
    given) is extracted from its ATOM records via ProDy -- not SEQRES, so this
    is exactly what's resolved in the deposited coordinates, no more and no
    less -- and aligned onto canonical's own numbering. The result is a set of
    equal-length, canonical-indexed sequences directly comparable position by
    position across every structure, with '-' marking any canonical position
    not observed in that structure (missing density, or outside the modeled
    construct). Rendering (text blocks, colored HTML, mutation lists, ...) is
    left to the caller; this function only produces the aligned data.

    accession: UniProt accession (e.g. "P07550").
    structures: list of PDB/CIF file paths (e.g. from
        chem.rcsb.download_structures).
    canonical_feature_type / canonical_feature_description: optional UniProt
        feature to slice canonical down to -- e.g. ("Chain", "Beta-lactamase
        TEM") to drop a cleaved signal peptide, or ("Domain", "Bromo 1") for a
        single domain of a larger protein. canonical_feature_type=None (the
        default) uses the full-length UniProt sequence as-is.
    marker_feature_types: optional tuple of UniProt feature types to collect as
        positions of interest -- e.g. ("Active site",) for an enzyme's
        catalytic residues, ("Binding site",) for a receptor's ligand pocket,
        ("Site",) for a non-catalytic functional residue. None (the default,
        matching a protein with no such annotation, or when markers just
        aren't wanted) returns an empty marker set.
    chain: optional chain id to use in every structure, overriding
        auto-selection. Default: whichever chain in each structure has the
        highest-scoring global alignment to canonical.

    Returns a dict:
        "protein_name", "organism": from the UniProt entry.
        "canonical_seq": the (possibly feature-sliced) canonical sequence.
        "feature_start", "feature_end": canonical_seq's 1-based bounds within
            the full UniProt sequence (1, len(full_seq) when
            canonical_feature_type is None).
        "marker_positions": set of 1-based canonical_seq positions collected
            from marker_feature_types (empty when marker_feature_types is
            None).
        "raw_sequences": {path: sequence}, the selected chain's CA-derived
            sequence exactly as observed, before alignment.
        "sequences": {path: sequence}, each raw sequence collapsed onto
            canonical's coordinates (see _align_to_canonical) -- for every
            path, len(sequences[path]) == len(canonical_seq).
    """
    entry, canonical_seq, feat_start, feat_end = _fetch_canonical(
        accession, canonical_feature_type, canonical_feature_description
    )
    marker_positions = _marker_positions(entry, feat_start, feat_end, marker_feature_types)

    aligner = _make_aligner()
    raw_sequences = {}
    sequences = {}
    for path in tqdm(structures, desc="aligning structures", unit="structure", disable=is_quiet()):
        raw_seq = _best_matching_chain_sequence(aligner, canonical_seq, path, chain)
        raw_sequences[path] = raw_seq
        sequences[path] = _align_to_canonical(aligner, canonical_seq, raw_seq)

    return {
        "protein_name": entry["proteinDescription"]["recommendedName"]["fullName"]["value"],
        "organism": entry["organism"]["scientificName"],
        "canonical_seq": canonical_seq,
        "feature_start": feat_start,
        "feature_end": feat_end,
        "marker_positions": marker_positions,
        "raw_sequences": raw_sequences,
        "sequences": sequences,
    }
