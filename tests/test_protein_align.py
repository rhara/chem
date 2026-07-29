import pytest

import chem.protein.structural_align as pa


def test_matched_ca_pairs_identical_sequences():
    ref_seq = "ACDEFGHIK"
    mob_seq = "ACDEFGHIK"
    ref_ca = list(range(len(ref_seq)))
    mob_ca = list(range(100, 100 + len(mob_seq)))
    ref_pts, mob_pts = pa._matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca)
    assert ref_pts == ref_ca
    assert mob_pts == mob_ca


def test_matched_ca_pairs_with_gap():
    ref_seq = "ACDEFGHIK"
    mob_seq = "ACDEHIK"  # FG deleted
    ref_ca = list(range(len(ref_seq)))
    mob_ca = list(range(100, 100 + len(mob_seq)))
    ref_pts, mob_pts = pa._matched_ca_pairs(ref_seq, ref_ca, mob_seq, mob_ca)
    assert len(ref_pts) == len(mob_pts) == 7
    # the two deleted positions (F, G at ref indices 4,5) must not appear
    assert 4 not in ref_pts
    assert 5 not in ref_pts


def test_align_rejects_empty_structures():
    with pytest.raises(ValueError):
        pa.align([])
