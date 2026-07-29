import pytest

import chem.rcsb.fetch as rf


def test_download_one_skips_existing_file(tmp_path, monkeypatch):
    (tmp_path / "1ABC.cif").write_bytes(b"already here")

    def fail_if_called(*args, **kwargs):
        raise AssertionError("should not hit the network for an existing file")

    monkeypatch.setattr(rf.requests, "get", fail_if_called)
    assert rf._download_one("1ABC", str(tmp_path), "cif") is True
    assert (tmp_path / "1ABC.cif").read_bytes() == b"already here"


def test_select_entries_no_threshold_keeps_everything():
    resolutions = {"1ABC": 2.1, "2XYZ": None, "3DEF": 4.5}
    assert set(rf._select_entries(resolutions, None)) == {"1ABC", "2XYZ", "3DEF"}


def test_select_entries_threshold_excludes_resolutionless_and_coarse():
    resolutions = {"1ABC": 2.1, "2XYZ": None, "3DEF": 4.5}
    assert rf._select_entries(resolutions, 3.0) == ["1ABC"]


def test_select_entries_threshold_is_inclusive():
    resolutions = {"1ABC": 2.0}
    assert rf._select_entries(resolutions, 2.0) == ["1ABC"]


def test_download_structures_rejects_bad_filetype():
    with pytest.raises(ValueError):
        rf.download_structures("P00734", filetype="mol2")


def test_validate_pdb_ids_normalizes_case():
    assert rf._validate_pdb_ids(["1abc", "2XYZ"]) == ["1ABC", "2XYZ"]


def test_validate_pdb_ids_rejects_bad_shape():
    with pytest.raises(ValueError):
        rf._validate_pdb_ids(["1ABC", "TOOLONG"])


def test_validate_pdb_ids_rejects_empty_list():
    with pytest.raises(ValueError):
        rf._validate_pdb_ids([])


def test_download_structures_with_explicit_pdb_ids_skips_resolution(tmp_path, monkeypatch):
    def fail_if_called(*args, **kwargs):
        raise AssertionError("should not resolve a target when explicit PDB ids are given")

    monkeypatch.setattr(rf, "resolve_uniprot_accession_any", fail_if_called)
    monkeypatch.setattr(rf, "_search_entry_ids", fail_if_called)
    monkeypatch.setattr(rf, "_fetch_resolutions", lambda ids: {i: None for i in ids})
    monkeypatch.setattr(rf, "_download_one", lambda entry_id, outdir, filetype: True)

    n = rf.download_structures(["1abc", "2xyz"], outdir=str(tmp_path))
    assert n == 2
