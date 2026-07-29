import chem.ids as ids


def test_resolve_chembl_id_passthrough():
    assert ids.resolve_target_chembl_id("CHEMBL204") == "CHEMBL204"
    assert ids.resolve_target_chembl_id("chembl204") == "CHEMBL204"


def test_uniprot_accession_regex():
    assert ids.UNIPROT_ACCESSION_RE.match("P00734")
    assert not ids.UNIPROT_ACCESSION_RE.match("THRB_HUMAN")
