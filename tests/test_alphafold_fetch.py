import pytest

import chem.alphafold.fetch as af


def test_select_entries_no_threshold_keeps_everything():
    entries = [{"entryId": "AF-1", "globalMetricValue": 40.0}, {"entryId": "AF-2", "globalMetricValue": 90.0}]
    assert af._select_entries(entries, None) == entries


def test_select_entries_threshold_filters_low_confidence():
    entries = [{"entryId": "AF-1", "globalMetricValue": 40.0}, {"entryId": "AF-2", "globalMetricValue": 90.0}]
    assert af._select_entries(entries, 70.0) == [entries[1]]


def test_select_entries_threshold_is_inclusive():
    entries = [{"entryId": "AF-1", "globalMetricValue": 70.0}]
    assert af._select_entries(entries, 70.0) == entries


def test_download_one_skips_existing_file(tmp_path, monkeypatch):
    (tmp_path / "AF-1.cif").write_bytes(b"already here")

    def fail_if_called(*args, **kwargs):
        raise AssertionError("should not hit the network for an existing file")

    monkeypatch.setattr(af.requests, "get", fail_if_called)
    entry = {"entryId": "AF-1", "cifUrl": "https://example.org/AF-1.cif"}
    assert af._download_one(entry, str(tmp_path), "cif") is True
    assert (tmp_path / "AF-1.cif").read_bytes() == b"already here"


def test_download_structures_rejects_bad_filetype():
    with pytest.raises(ValueError):
        af.download_structures("P00734", filetype="mol2")
