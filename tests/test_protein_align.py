import os

import numpy as np
import pytest

import chem.protein.structural_align as pa


def _atom_line(serial, name, resname, chain, resseq, x, y, z, record="ATOM"):
    element = name.strip()[0]
    return (
        f"{record:<6}{serial:>5} {name:<4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element:>2}"
    )


def _write_three_residue_chain(path, chain, x_offset=0.0):
    lines = [
        _atom_line(1, "CA", "ALA", chain, 1, 0.0 + x_offset, 0.0, 0.0),
        _atom_line(2, "CA", "GLY", chain, 2, 3.8 + x_offset, 0.0, 0.0),
        _atom_line(3, "CA", "LEU", chain, 3, 7.6 + x_offset, 0.0, 0.0),
    ]
    path.write_text("\n".join(lines) + "\nEND\n")


def test_align_writes_reference_even_when_not_in_structures(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    mobile_path = tmp_path / "mobile.pdb"
    _write_three_residue_chain(ref_path, "A")
    _write_three_residue_chain(mobile_path, "A", x_offset=10.0)

    outdir = tmp_path / "aligned"
    results = pa.align([str(mobile_path)], reference=str(ref_path), outdir=str(outdir))

    assert sorted(os.listdir(outdir)) == ["mobile.pdb", "ref.pdb"]
    assert results[str(ref_path)] == {"rmsd": 0.0, "identity": 1.0}
    assert results[str(mobile_path)]["rmsd"] == pytest.approx(0.0, abs=1e-6)
    assert results[str(mobile_path)]["identity"] == 1.0


def test_align_rmsd_is_plain_float_rounded_to_3dp(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    mobile_path = tmp_path / "mobile.pdb"
    _write_three_residue_chain(ref_path, "A")
    # A small perturbation (not just a rigid translation) so the best fit leaves
    # a nonzero, not-conveniently-round residual to check the rounding on.
    lines = [
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "GLY", "A", 2, 3.8, 0.1234, 0.0),
        _atom_line(3, "CA", "LEU", "A", 3, 7.6, 0.0, 0.0),
    ]
    mobile_path.write_text("\n".join(lines) + "\nEND\n")

    results = pa.align([str(mobile_path)], reference=str(ref_path), outdir=str(tmp_path / "aligned"))

    rmsd = results[str(mobile_path)]["rmsd"]
    assert type(rmsd) is float  # not numpy.float64
    assert rmsd != 0.0
    assert rmsd == round(rmsd, 3)
    identity = results[str(mobile_path)]["identity"]
    assert type(identity) is float
    assert identity == 1.0  # same residue types, only coordinates perturbed


def test_chain_seq_and_ca_excludes_hetatm_with_standard_resname(tmp_path):
    # A covalently-linked peptidomimetic ligand can share the protein's chain id
    # and use a standard amino acid resname (e.g. "PRO") but as a HETATM -- as
    # with real thrombin-inhibitor complex 6YHG, whose bound ligand's "PRO H 307"
    # was wrongly pulled into chain H's sequence, corrupting the alignment.
    lines = [
        _atom_line(1, "CA", "ALA", "H", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "GLY", "H", 2, 3.8, 0.0, 0.0),
        _atom_line(3, "CA", "LEU", "H", 3, 7.6, 0.0, 0.0),
        _atom_line(4, "CA", "PRO", "H", 307, 500.0, 500.0, 500.0, record="HETATM"),
    ]
    path = tmp_path / "ligand_hetatm.pdb"
    path.write_text("\n".join(lines) + "\nEND\n")

    structure = pa._load_structure(str(path))
    chain = pa._select_chain(next(structure.get_models()))
    seq, ca_atoms = pa._chain_seq_and_ca(chain)
    assert seq == "AGL"
    assert len(ca_atoms) == 3


def test_matched_ca_pairs_identical_sequences():
    ref_seq = "ACDEFGHIK"
    mob_seq = "ACDEFGHIK"
    ref_ca = list(range(len(ref_seq)))
    mob_ca = list(range(100, 100 + len(mob_seq)))
    ref_pts, mob_pts, identity = pa._matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca)
    assert ref_pts == ref_ca
    assert mob_pts == mob_ca
    assert identity == 1.0


def test_matched_ca_pairs_with_gap():
    ref_seq = "ACDEFGHIK"
    mob_seq = "ACDEHIK"  # FG deleted
    ref_ca = list(range(len(ref_seq)))
    mob_ca = list(range(100, 100 + len(mob_seq)))
    ref_pts, mob_pts, identity = pa._matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca)
    assert len(ref_pts) == len(mob_pts) == 7
    # the two deleted positions (F, G at ref indices 4,5) must not appear
    assert 4 not in ref_pts
    assert 5 not in ref_pts
    # every remaining (gap-free) position is still an exact match
    assert identity == 1.0


def test_matched_ca_pairs_identity_excludes_gaps_counts_mismatches():
    ref_seq = "ACDEFGHIK"
    mob_seq = "ACDXFGHIK"  # one substitution (E->X), no length change/gap
    ref_ca = list(range(len(ref_seq)))
    mob_ca = list(range(100, 100 + len(mob_seq)))
    ref_pts, mob_pts, identity = pa._matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca)
    # no gap opened (mismatch is cheaper than a gap here), so all 9 positions match
    assert len(ref_pts) == 9
    # ...but only 8 of those 9 matched positions are an identical residue
    assert identity == pytest.approx(8 / 9)


def test_align_rejects_empty_structures():
    with pytest.raises(ValueError):
        pa.align([])


_THREE_LETTER = {
    "A": "ALA", "G": "GLY", "L": "LEU", "V": "VAL", "P": "PRO",
    "S": "SER", "T": "THR", "C": "CYS", "N": "ASN", "Q": "GLN",
    "M": "MET", "I": "ILE", "F": "PHE", "W": "TRP", "H": "HIS",
    "K": "LYS", "R": "ARG", "D": "ASP", "E": "GLU", "Y": "TYR",
}


def _write_chain(path, chain, one_letter_seq, x_offset=0.0, mode="w"):
    lines = [
        _atom_line(i + 1, "CA", _THREE_LETTER[c], chain, i + 1, i * 3.8 + x_offset, 0.0, 0.0)
        for i, c in enumerate(one_letter_seq)
    ]
    with open(path, mode) as f:
        f.write("\n".join(lines) + "\nEND\n" if mode == "w" else "\n" + "\n".join(lines) + "\nEND\n")


def test_align_picks_matching_chain_over_larger_unrelated_chain(tmp_path):
    # Reproduces a real bug seen with RCSB entry 3B9F (thrombin co-crystallized
    # with protein C inhibitor, a serpin): a chain-size-only heuristic picked
    # the larger bound partner protein (356 residues) over thrombin's own
    # heavy chain (253 residues), giving ~0.22 identity instead of ~0.996.
    ref_path = tmp_path / "ref.pdb"
    _write_chain(ref_path, "A", "AGLVPSTCNQ")  # 10 residues

    mobile_path = tmp_path / "mobile.pdb"
    # Chain S ("small"): 5 residues, exact prefix match -- the true counterpart.
    _write_chain(mobile_path, "S", "AGLVP", mode="w")
    # Chain B ("big"): 10 residues (same length as ref, so no gap is needed),
    # only one incidental match (index 3, both "V") -- an unrelated chain that
    # a matched-length-only heuristic would wrongly prefer.
    _write_chain(mobile_path, "B", "MIFVHKRDEY", mode="a")

    results = pa.align([str(mobile_path)], reference=str(ref_path), outdir=str(tmp_path / "aligned"))

    assert results[str(mobile_path)]["identity"] == 1.0
    assert results[str(mobile_path)]["rmsd"] == pytest.approx(0.0, abs=1e-6)


def test_best_matching_chain_alignment_prefers_identity_over_raw_matched_count(tmp_path):
    ref_seq = "AGLVPSTCNQ"
    ref_ca = list(range(len(ref_seq)))
    small_seq, big_seq = "AGLVP", "MIFVHKRDEY"

    path = tmp_path / "mobile.pdb"
    _write_chain(path, "S", small_seq, mode="w")
    _write_chain(path, "B", big_seq, mode="a")
    model = next(pa._load_structure(str(path)).get_models())

    ref_pts, mob_pts, identity = pa._best_matching_chain_alignment(model, ref_seq, ref_ca)

    assert identity == 1.0
    assert len(ref_pts) == len(small_seq)


def test_identity_matrix_self_identity_is_one(tmp_path):
    path = tmp_path / "a.pdb"
    _write_chain(path, "A", "AGLVPSTCNQ")

    matrix = pa.identity_matrix([str(path)])

    assert matrix == {str(path): {str(path): 1.0}}


def test_identity_matrix_symmetric_and_values(tmp_path):
    path_a = tmp_path / "a.pdb"
    path_b = tmp_path / "b.pdb"  # same sequence as a -- should be identity 1.0
    path_c = tmp_path / "c.pdb"  # unrelated sequence
    _write_chain(path_a, "A", "AGLVPSTCNQ")
    _write_chain(path_b, "A", "AGLVPSTCNQ")
    _write_chain(path_c, "A", "MIFVHKRDEY")

    matrix = pa.identity_matrix([str(path_a), str(path_b), str(path_c)])

    a, b, c = str(path_a), str(path_b), str(path_c)
    assert matrix[a][a] == matrix[b][b] == matrix[c][c] == 1.0
    assert matrix[a][b] == matrix[b][a] == 1.0
    assert matrix[a][c] == matrix[c][a]
    assert matrix[b][c] == matrix[c][b]
    assert matrix[a][c] < 1.0


def test_identity_matrix_skips_structure_with_no_polymer_residues(tmp_path):
    good_path = tmp_path / "good.pdb"
    _write_chain(good_path, "A", "AGLVPSTCNQ")

    bad_path = tmp_path / "no_polymer.pdb"
    lines = [_atom_line(1, "C1", "LIG", "A", 1, 0.0, 0.0, 0.0, record="HETATM")]
    bad_path.write_text("\n".join(lines) + "\nEND\n")

    matrix = pa.identity_matrix([str(good_path), str(bad_path)])

    assert list(matrix) == [str(good_path)]
    assert matrix == {str(good_path): {str(good_path): 1.0}}


def test_identity_matrix_respects_explicit_chain(tmp_path):
    # Chain "S" differs between the two files (and is the longer chain, so a
    # size-based default heuristic would pick it); chain "B" is identical between
    # them. Only if chain="B" is actually honored (not silently ignored in favor
    # of auto-selecting the larger chain) does the pair come out at identity 1.0.
    path_1 = tmp_path / "one.pdb"
    _write_chain(path_1, "S", "AGLVPSTCNQKR", mode="w")
    _write_chain(path_1, "B", "MIFVHKRDEY", mode="a")

    path_2 = tmp_path / "two.pdb"
    _write_chain(path_2, "S", "WYHDEKRQNCT", mode="w")
    _write_chain(path_2, "B", "MIFVHKRDEY", mode="a")

    matrix = pa.identity_matrix([str(path_1), str(path_2)], chain="B")

    assert matrix[str(path_1)][str(path_2)] == 1.0


def test_compute_transform_identity_when_already_superposed(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    mobile_path = tmp_path / "mobile.pdb"
    _write_three_residue_chain(ref_path, "A")
    _write_three_residue_chain(mobile_path, "A")

    transform = pa.compute_transform(str(mobile_path), str(ref_path))

    assert transform["rmsd"] == pytest.approx(0.0, abs=1e-6)
    assert transform["identity"] == 1.0
    np.testing.assert_allclose(transform["rotation"], np.eye(3), atol=1e-5)
    np.testing.assert_allclose(transform["translation"], [0.0, 0.0, 0.0], atol=1e-5)


def test_compute_transform_recovers_pure_translation(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    mobile_path = tmp_path / "mobile.pdb"
    _write_three_residue_chain(ref_path, "A")
    _write_three_residue_chain(mobile_path, "A", x_offset=10.0)

    transform = pa.compute_transform(str(mobile_path), str(ref_path))

    assert transform["rmsd"] == pytest.approx(0.0, abs=1e-6)
    np.testing.assert_allclose(transform["rotation"], np.eye(3), atol=1e-5)
    np.testing.assert_allclose(transform["translation"], [-10.0, 0.0, 0.0], atol=1e-5)


def test_compute_transform_rejects_too_few_matched_residues(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    mobile_path = tmp_path / "mobile.pdb"
    _write_three_residue_chain(ref_path, "A")
    _write_chain(mobile_path, "A", "AG")  # only 2 residues, below _MIN_MATCHED_RESIDUES

    with pytest.raises(ValueError):
        pa.compute_transform(str(mobile_path), str(ref_path))


def test_compute_transform_honors_explicit_chain_ids(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    _write_chain(ref_path, "A", "AGLVPSTCNQ")

    mobile_path = tmp_path / "mobile.pdb"
    _write_chain(mobile_path, "S", "MIFVHKRDEY", mode="w")  # unrelated, larger index
    _write_chain(mobile_path, "A", "AGLVPSTCNQ", mode="a", x_offset=5.0)

    transform = pa.compute_transform(
        str(mobile_path), str(ref_path), mobile_chain="A", reference_chain="A"
    )

    assert transform["identity"] == 1.0
    np.testing.assert_allclose(transform["translation"], [-5.0, 0.0, 0.0], atol=1e-5)


def test_apply_transform_moves_every_atom_and_matches_compute_transform(tmp_path):
    ref_path = tmp_path / "ref.pdb"
    mobile_path = tmp_path / "mobile.pdb"
    _write_three_residue_chain(ref_path, "A")
    _write_three_residue_chain(mobile_path, "A", x_offset=10.0)

    transform = pa.compute_transform(str(mobile_path), str(ref_path))
    out_path = tmp_path / "moved.pdb"
    result_path = pa.apply_transform(
        str(mobile_path), transform["rotation"], transform["translation"], str(out_path)
    )

    assert result_path == str(out_path)
    structure = pa._load_structure(str(out_path))
    coords = np.array([atom.coord for atom in structure.get_atoms()])
    ref_structure = pa._load_structure(str(ref_path))
    ref_coords = np.array([atom.coord for atom in ref_structure.get_atoms()])
    np.testing.assert_allclose(coords, ref_coords, atol=1e-4)


def test_apply_transform_creates_missing_outdir(tmp_path):
    path = tmp_path / "a.pdb"
    _write_chain(path, "A", "AGL")

    out_path = tmp_path / "nested" / "out.pdb"
    pa.apply_transform(str(path), np.eye(3).tolist(), [0.0, 0.0, 0.0], str(out_path))

    assert out_path.exists()
