import os

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
