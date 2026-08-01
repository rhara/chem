import chem.ids as ids


def test_resolve_chembl_id_passthrough():
    assert ids.resolve_target_chembl_id("CHEMBL204") == "CHEMBL204"
    assert ids.resolve_target_chembl_id("chembl204") == "CHEMBL204"


def test_uniprot_accession_regex():
    assert ids.UNIPROT_ACCESSION_RE.match("P00734")
    assert not ids.UNIPROT_ACCESSION_RE.match("THRB_HUMAN")


def test_resolve_uniprot_accession_prefers_exact_id_match_over_top_ranked_hit(monkeypatch):
    # Regression test: UniProt's "id" query is a relevance-ranked text search, not an
    # exact-match lookup. "id:CDK1_HUMAN" really does rank O14519 (CDKA1_HUMAN, an
    # unrelated protein) ahead of the actual exact match P06493 (CDK1_HUMAN) -- so the
    # exact match must be picked even when it isn't the first result.
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "results": [
                    {"primaryAccession": "O14519", "uniProtkbId": "CDKA1_HUMAN"},
                    {"primaryAccession": "P06493", "uniProtkbId": "CDK1_HUMAN"},
                ]
            }

    monkeypatch.setattr(ids.requests, "get", lambda url, params, timeout: FakeResponse())

    assert ids.resolve_uniprot_accession("CDK1_HUMAN") == "P06493"


def test_resolve_uniprot_accession_falls_back_to_top_hit_without_exact_match(monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {"results": [{"primaryAccession": "P00734", "uniProtkbId": "THRB_HUMAN"}]}

    monkeypatch.setattr(ids.requests, "get", lambda url, params, timeout: FakeResponse())

    assert ids.resolve_uniprot_accession("some-deprecated-alias") == "P00734"
