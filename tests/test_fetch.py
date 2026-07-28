import chem.chembl.fetch as cf


def test_resolve_chembl_id_passthrough():
    assert cf._resolve_target_chembl_id("CHEMBL204") == "CHEMBL204"
    assert cf._resolve_target_chembl_id("chembl204") == "CHEMBL204"


def test_normalize_smiles_strips_hcl_salt():
    # hydrochloride salt: keep the parent, drop the HCl counterion
    smi = cf._normalize_smiles("CN=C(N)N1CCC(Oc2ccccc2)CC1.Cl")
    assert smi == "CN=C(N)N1CCC(Oc2ccccc2)CC1"


def test_normalize_smiles_single_component():
    smi = cf._normalize_smiles("CCO")
    assert smi == "CCO"


def test_normalize_smiles_invalid():
    assert cf._normalize_smiles("not a smiles") is None


def test_aggregate_by_smiles():
    records = [
        {"molecule_chembl_id": "CHEMBL1", "target_chembl_id": "CHEMBL204",
         "smiles": "CCO", "mw": 46.07, "pchembl_value": 5.0},
        {"molecule_chembl_id": "CHEMBL2", "target_chembl_id": "CHEMBL204",
         "smiles": "CCO", "mw": 46.07, "pchembl_value": 7.0},
        {"molecule_chembl_id": "CHEMBL3", "target_chembl_id": "CHEMBL204",
         "smiles": "CCN", "mw": 45.08, "pchembl_value": 6.0},
    ]
    by_smiles = {row[2]: row for row in cf._aggregate_by_smiles(records)}

    cco_row = by_smiles["CCO"]
    assert cco_row[0] == "CHEMBL1"  # parent_chembl_id: first record in the group
    assert cco_row[4] == 2  # n
    assert cco_row[5] == 6.0  # mean
    assert cco_row[6] == 6.0  # median
    assert cco_row[7] == 1.0  # population std of [5.0, 7.0]

    ccn_row = by_smiles["CCN"]
    assert ccn_row[4] == 1
    assert ccn_row[7] == 0.0  # population std of a single value
