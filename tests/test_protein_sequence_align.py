import importlib

# chem/protein/__init__.py does `from .sequence_align import sequence_align`,
# which rebinds the chem.protein.sequence_align *attribute* to the function
# (same name as the submodule) -- so `import chem.protein.sequence_align as sa`
# would resolve to the function, not the module. Go through importlib instead
# to get the actual module and its private helpers.
sa = importlib.import_module("chem.protein.sequence_align")

# Trimmed-down but structurally faithful stand-in for a UniProt entry with a
# cleaved signal peptide (canonical_feature_type slicing) and a couple of
# "Active site" features (marker_feature_types).
_FAKE_ENTRY = {
    "proteinDescription": {"recommendedName": {"fullName": {"value": "Fake Enzyme"}}},
    "organism": {"scientificName": "Homo sapiens"},
    "sequence": {"value": "MKKAGLVPSTCNQKR"},  # 15 aa; mature chain is residues 3-15
    "features": [
        {"type": "Chain", "description": "Fake Enzyme", "location": {"start": {"value": 3}, "end": {"value": 15}}},
        {"type": "Active site", "description": "Nucleophile", "location": {"start": {"value": 5}, "end": {"value": 5}}},
        {"type": "Active site", "description": "Proton donor", "location": {"start": {"value": 12}, "end": {"value": 12}}},
    ],
}


def _atom_line(serial, name, resname, chain, resseq, x, y, z):
    element = name.strip()[0]
    return (
        f"ATOM  {serial:>5} {name:<4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element:>2}"
    )


_THREE_LETTER = {
    "A": "ALA", "G": "GLY", "L": "LEU", "V": "VAL", "P": "PRO",
    "S": "SER", "T": "THR", "C": "CYS", "N": "ASN", "Q": "GLN",
    "M": "MET", "I": "ILE", "F": "PHE", "W": "TRP", "H": "HIS",
    "K": "LYS", "R": "ARG", "D": "ASP", "E": "GLU", "Y": "TYR",
}


def _write_chain(f, chain, one_letter_seq, x_offset=0.0):
    for i, c in enumerate(one_letter_seq):
        f.write(_atom_line(i + 1, "CA", _THREE_LETTER[c], chain, i + 1, i * 3.8 + x_offset, 0.0, 0.0) + "\n")


def _write_pdb(path, chains):
    """chains: list of (chain_id, one_letter_seq) tuples."""
    with open(path, "w") as f:
        for i, (chain_id, seq) in enumerate(chains):
            _write_chain(f, chain_id, seq, x_offset=100.0 * i)
        f.write("END\n")


def test_marker_positions_none_returns_empty_set():
    assert sa._marker_positions(_FAKE_ENTRY, 3, 15, None) == set()
    assert sa._marker_positions(_FAKE_ENTRY, 3, 15, ()) == set()


def test_marker_positions_filters_by_type_and_rebases_to_slice():
    # precursor positions 5, 12 rebased to the mature-chain slice starting at 3
    # -> mature-relative positions 3, 10
    positions = sa._marker_positions(_FAKE_ENTRY, 3, 15, ("Active site",))
    assert positions == {3, 10}


def test_fetch_canonical_slices_by_feature(monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return _FAKE_ENTRY

    monkeypatch.setattr(sa.requests, "get", lambda url, timeout: FakeResponse())

    entry, canonical_seq, start, end = sa._fetch_canonical("Q00000", "Chain", "Fake Enzyme")

    assert canonical_seq == "KAGLVPSTCNQKR"  # full_seq[2:15]
    assert (start, end) == (3, 15)


def test_fetch_canonical_full_length_when_no_feature_type(monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return _FAKE_ENTRY

    monkeypatch.setattr(sa.requests, "get", lambda url, timeout: FakeResponse())

    _entry, canonical_seq, start, end = sa._fetch_canonical("Q00000", None, None)

    assert canonical_seq == _FAKE_ENTRY["sequence"]["value"]
    assert (start, end) == (1, 15)


def test_align_to_canonical_drops_foreign_insertion_as_a_clean_gap():
    aligner = sa._make_aligner()
    canonical = "AGLV" + "PSTCNQ" + "KRHWY"  # prefix(4) | replaced(6) | suffix(5)
    # A foreign, unrelated block spliced in where canonical has "PSTCNQ", like
    # a fusion partner -- should collapse to a gap at exactly those six
    # canonical positions, not be scattered across nearby mismatches. (The
    # replaced span needs to be long enough for this: see
    # test_align_to_canonical_can_smear_a_short_replaced_span for the known
    # edge case where a *short* replaced span doesn't resolve this cleanly.)
    query = "AGLV" + "WYHDEKRQNCTMIFVDE" + "KRHWY"
    observed = sa._align_to_canonical(aligner, canonical, query)
    assert observed == "AGLV------KRHWY"


def test_align_to_canonical_can_smear_a_short_replaced_span():
    # Known, accepted limitation (see the comment on mismatch_score in
    # _make_aligner): when the canonical span being replaced by a foreign
    # block is short (here, 2 residues), absorbing it into a couple of
    # mismatches is cheaper than paying for a second gap-open, so it does NOT
    # collapse to a clean gap the way a longer replaced span does. This test
    # documents/locks in that known behavior rather than treating it as a
    # regression.
    aligner = sa._make_aligner()
    canonical = "AGLV" + "PS" + "TCNQKR"  # prefix(4) | replaced(2) | suffix(6)
    query = "AGLV" + "WYHDEKRQNCTMIFV" + "TCNQKR"
    observed = sa._align_to_canonical(aligner, canonical, query)
    assert "-" not in observed  # smeared instead of gapped
    assert observed != "AGLV--TCNQKR"


def test_align_to_canonical_prefers_fewest_blocks_on_tie():
    # Canonical ends in a repeated residue ("...NQKK"); a query that only
    # covers the first half should collapse to one contiguous trailing gap
    # (not a coincidental match to the far-away final "K").
    aligner = sa._make_aligner()
    canonical = "AGLVPSTCNQKK"
    query = "AGLVPSTC"  # matches canonical[:8] exactly, nothing beyond
    observed = sa._align_to_canonical(aligner, canonical, query)
    assert observed == "AGLVPSTC----"


def test_load_ca_sequences_reads_all_chains(tmp_path):
    path = tmp_path / "multi.pdb"
    _write_pdb(path, [("A", "AGLVP"), ("B", "MIFVHKR")])

    sequences = sa._load_ca_sequences(str(path))

    assert sequences == {"A": "AGLVP", "B": "MIFVHKR"}


def test_best_matching_chain_sequence_picks_highest_scoring_chain(tmp_path):
    # Chain B is an unrelated (larger) fusion partner; chain A is the real match.
    path = tmp_path / "fusion.pdb"
    _write_pdb(path, [("B", "MIFVHKRDEYQNCT"), ("A", "AGLVPSTCNQKR")])
    aligner = sa._make_aligner()

    best = sa._best_matching_chain_sequence(aligner, "AGLVPSTCNQKR", str(path), chain=None)

    assert best == "AGLVPSTCNQKR"


def test_best_matching_chain_sequence_honors_explicit_chain(tmp_path):
    path = tmp_path / "fusion.pdb"
    _write_pdb(path, [("A", "AGLVPSTCNQKR"), ("B", "MIFVHKRDEYQNCT")])
    aligner = sa._make_aligner()

    best = sa._best_matching_chain_sequence(aligner, "AGLVPSTCNQKR", str(path), chain="B")

    assert best == "MIFVHKRDEYQNCT"


def test_sequence_align_end_to_end(tmp_path, monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return _FAKE_ENTRY

    monkeypatch.setattr(sa.requests, "get", lambda url, timeout: FakeResponse())

    # canonical mature chain (residues 3-15) is "KAGLVPSTCNQKR"; this structure
    # is missing the leading "K" (residue 3) and has a point mutation C->W
    # (canonical position 8, 0-based index 7).
    mature = _FAKE_ENTRY["sequence"]["value"][2:]  # "KAGLVPSTCNQKR"
    observed_construct = mature[1:].replace("C", "W", 1)  # drop res 3, mutate the C
    path = tmp_path / "struct.pdb"
    _write_pdb(path, [("A", observed_construct)])

    result = sa.sequence_align(
        "Q00000",
        [str(path)],
        canonical_feature_type="Chain",
        canonical_feature_description="Fake Enzyme",
        marker_feature_types=("Active site",),
    )

    assert result["protein_name"] == "Fake Enzyme"
    assert result["organism"] == "Homo sapiens"
    assert result["canonical_seq"] == mature
    assert (result["feature_start"], result["feature_end"]) == (3, 15)
    assert result["marker_positions"] == {3, 10}
    assert result["raw_sequences"][str(path)] == observed_construct
    aligned = result["sequences"][str(path)]
    assert len(aligned) == len(mature)
    assert aligned[0] == "-"  # the dropped leading K
    assert aligned[1:] == observed_construct
