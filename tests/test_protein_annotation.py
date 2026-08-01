import chem.protein.annotation as ann

# Trimmed-down but structurally faithful stand-in for
# https://rest.uniprot.org/uniprotkb/Q8IZL9.json
_FAKE_ENTRY = {
    "primaryAccession": "Q8IZL9",
    "uniProtkbId": "CDK20_HUMAN",
    "organism": {"scientificName": "Homo sapiens"},
    "proteinExistence": "1: Evidence at protein level",
    "annotationScore": 5.0,
    "proteinDescription": {
        "recommendedName": {
            "fullName": {"value": "Cyclin-dependent kinase 20"},
            "ecNumbers": [{"value": "2.7.11.22"}],
        }
    },
    "genes": [{"geneName": {"value": "CDK20"}, "synonyms": [{"value": "CCRK"}, {"value": "CDCH"}]}],
    "comments": [
        {
            "commentType": "SIMILARITY",
            "texts": [{"value": "Belongs to the protein kinase superfamily."}],
        },
        {
            "commentType": "FUNCTION",
            "texts": [{"value": "Activates CDK2 by phosphorylating residue 'Thr-160'"}],
        },
        {
            "commentType": "SUBCELLULAR LOCATION",
            "subcellularLocations": [
                {"location": {"value": "Nucleus"}},
                {"location": {"value": "Cytoplasm"}},
            ],
        },
    ],
    "features": [
        {"type": "Domain", "description": "Protein kinase", "location": {"start": {"value": 4}, "end": {"value": 288}}},
        {"type": "Active site", "description": "Proton acceptor", "location": {"start": {"value": 127}, "end": {"value": 127}}},
    ],
    "sequence": {"length": 346},
    "uniProtKBCrossReferences": [
        {"database": "AlphaFoldDB", "id": "Q8IZL9"},
        {"database": "BindingDB", "id": "Q8IZL9"},
        {"database": "ChEMBL", "id": "CHEMBL3559690"},
        {"database": "Pharos", "id": "Q8IZL9", "properties": [{"key": "DevelopmentLevel", "value": "Tbio"}]},
    ],
}


def test_extract_properties_pulls_expected_fields():
    props = ann._extract_properties(_FAKE_ENTRY)

    assert props["entry_name"] == "CDK20_HUMAN"
    assert props["accession"] == "Q8IZL9"
    assert props["protein_name"] == "Cyclin-dependent kinase 20"
    assert props["gene_name"] == "CDK20 (synonyms: CCRK, CDCH)"
    assert props["organism"] == "Homo sapiens"
    assert props["sequence_length"] == 346
    assert props["ec_number"] == "2.7.11.22"
    assert props["kinase_domain_range"] == "4-288"
    assert props["active_site_residue"] == 127
    assert props["n_pdb_xrefs"] == 0
    assert props["has_alphafold_model"] is True
    assert props["has_bindingdb_entry"] is True
    assert props["chembl_target_id"] == "CHEMBL3559690"
    assert props["pharos_development_level"] == "Tbio (biology characterized, no known drug/chemical probe)"
    assert props["protein_existence"] == "1: Evidence at protein level"
    assert props["annotation_score"] == 5.0


def test_extract_properties_handles_missing_optional_fields():
    entry = dict(_FAKE_ENTRY)
    entry["genes"] = [{"geneName": {"value": "CDK20"}}]  # no synonyms
    entry["features"] = []  # no domain/active-site annotation
    entry["uniProtKBCrossReferences"] = []  # no cross-references at all

    props = ann._extract_properties(entry)

    assert props["gene_name"] == "CDK20"
    assert props["kinase_domain_range"] is None
    assert props["active_site_residue"] is None
    assert props["n_pdb_xrefs"] == 0
    assert props["has_alphafold_model"] is False
    assert props["chembl_target_id"] is None
    assert props["pharos_development_level"] is None


def test_summary_resolves_id_and_fetches_uniprot_json(monkeypatch):
    seen = {}

    def fake_resolve(id_):
        seen["resolved_with"] = id_
        return "Q8IZL9"

    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return _FAKE_ENTRY

    def fake_get(url, timeout):
        seen["url"] = url
        return FakeResponse()

    monkeypatch.setattr(ann, "resolve_uniprot_accession_any", fake_resolve)
    monkeypatch.setattr(ann.requests, "get", fake_get)

    props = ann.summary("CDK20_HUMAN")

    assert seen["resolved_with"] == "CDK20_HUMAN"
    assert seen["url"] == f"{ann.UNIPROT_API}/Q8IZL9.json"
    assert props["entry_name"] == "CDK20_HUMAN"
    assert props["accession"] == "Q8IZL9"


def test_summary_accepts_a_chembl_target_id(monkeypatch):
    seen = {}

    def fake_resolve(id_):
        seen["resolved_with"] = id_
        return "Q8IZL9"

    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return _FAKE_ENTRY

    monkeypatch.setattr(ann, "resolve_uniprot_accession_any", fake_resolve)
    monkeypatch.setattr(ann.requests, "get", lambda url, timeout: FakeResponse())

    props = ann.summary("CHEMBL3559690")

    assert seen["resolved_with"] == "CHEMBL3559690"
    assert props["accession"] == "Q8IZL9"


_FAKE_FASTA = ">sp|Q8IZL9|CDK20_HUMAN Cyclin-dependent kinase 20\nMDQYCILGRIG\n"


def test_get_fasta_resolves_id_and_fetches_uniprot_fasta(monkeypatch):
    seen = {}

    def fake_resolve(id_):
        seen["resolved_with"] = id_
        return "Q8IZL9"

    class FakeResponse:
        text = _FAKE_FASTA

        def raise_for_status(self):
            pass

    def fake_get(url, headers, timeout):
        seen["url"] = url
        seen["headers"] = headers
        return FakeResponse()

    monkeypatch.setattr(ann, "resolve_uniprot_accession_any", fake_resolve)
    monkeypatch.setattr(ann.requests, "get", fake_get)

    fasta = ann.get_fasta("CDK20_HUMAN")

    assert seen["resolved_with"] == "CDK20_HUMAN"
    assert seen["url"] == f"{ann.UNIPROT_API}/Q8IZL9.fasta"
    assert "user@example.com" in seen["headers"]["User-Agent"]
    assert fasta == _FAKE_FASTA


def test_get_fasta_accepts_a_chembl_target_id_and_custom_email(monkeypatch):
    seen = {}

    monkeypatch.setattr(ann, "resolve_uniprot_accession_any", lambda id_: "Q8IZL9")

    def fake_get(url, headers, timeout):
        seen["headers"] = headers
        return type("R", (), {"text": _FAKE_FASTA, "raise_for_status": lambda self: None})()

    monkeypatch.setattr(ann.requests, "get", fake_get)

    fasta = ann.get_fasta("CHEMBL3559690", email="me@example.org")

    assert "me@example.org" in seen["headers"]["User-Agent"]
    assert fasta == _FAKE_FASTA
