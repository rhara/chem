import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

import chem.ligand.extract as le


def _atom_line(serial, name, resname, chain, resseq, x, y, z, record="HETATM"):
    element = name.strip()[0]
    return (
        f"{record:<6}{serial:>5} {name:<4} {resname:>3} {chain:1}{resseq:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element:>2}"
    )


def _write_pdb(path, lines):
    path.write_text("\n".join(lines) + "\nEND\n")


@pytest.fixture
def structure_with_ligands(tmp_path):
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0, record="ATOM"),
        _atom_line(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0, record="ATOM"),
        _atom_line(3, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0),
        _atom_line(4, "C2", "LIG", "A", 200, 11.0, 10.0, 10.0),
        _atom_line(5, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0),
        _atom_line(6, "NA", "NA", "A", 301, 30.0, 30.0, 30.0),
        _atom_line(7, "C1", "OTH", "B", 1, 40.0, 40.0, 40.0),
    ]
    path = tmp_path / "structure.pdb"
    _write_pdb(path, lines)
    return str(path)


def test_list_ligand_codes_excludes_solvent_and_ions(structure_with_ligands):
    assert le.list_ligand_codes(structure_with_ligands) == ["LIG", "OTH"]


def test_list_ligand_codes_custom_exclude(structure_with_ligands):
    codes = le.list_ligand_codes(structure_with_ligands, exclude={"LIG"} | le.SOLVENT_AND_IONS)
    assert codes == ["OTH"]


def test_list_ligand_codes_empty_when_only_solvent(tmp_path):
    lines = [_atom_line(1, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0)]
    path = tmp_path / "no_ligand.pdb"
    _write_pdb(path, lines)
    assert le.list_ligand_codes(str(path)) == []


def test_pick_ligand_residue_not_found(structure_with_ligands):
    with pytest.raises(ValueError):
        le._pick_ligand_residue(structure_with_ligands, "ZZZ")


def test_fetch_template_smiles_candidates_prioritizes_openeye(monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "pdbx_chem_comp_descriptor": [
                    {"program": "CACTVS", "type": "SMILES_CANONICAL", "descriptor": "CACTVS_SMILES"},
                    {"program": "OpenEye OEToolkits", "type": "SMILES", "descriptor": "OE_SMILES"},
                ]
            }

    monkeypatch.setattr(le.requests, "get", lambda *a, **k: FakeResponse())
    candidates = le._fetch_template_smiles_candidates("XXX")
    assert candidates[0] == "OE_SMILES"
    assert "CACTVS_SMILES" in candidates


def test_fetch_template_smiles_candidates_raises_when_no_descriptors(monkeypatch):
    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {"pdbx_chem_comp_descriptor": []}

    monkeypatch.setattr(le.requests, "get", lambda *a, **k: FakeResponse())
    with pytest.raises(ValueError):
        le._fetch_template_smiles_candidates("XXX")


def _benzene_pdb(tmp_path):
    """A benzene ring written out with real 3D coordinates but, like an
    RCSB-downloaded ligand, no bond-order information -- RDKit's PDB reader
    will see only distance-guessed single bonds, no aromaticity."""
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    for atom in mol.GetAtoms():
        info = Chem.AtomPDBResidueInfo()
        info.SetResidueName("BNZ")
        info.SetResidueNumber(1)
        info.SetChainId("A")
        info.SetIsHeteroAtom(True)
        atom.SetMonomerInfo(info)
    path = tmp_path / "benzene.pdb"
    Chem.MolToPDBFile(mol, str(path))
    return str(path)


def test_load_ligand_restores_aromaticity(tmp_path, monkeypatch):
    path = _benzene_pdb(tmp_path)

    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "pdbx_chem_comp_descriptor": [
                    {"program": "CACTVS", "type": "SMILES_CANONICAL", "descriptor": "c1ccccc1"},
                ]
            }

    monkeypatch.setattr(le.requests, "get", lambda *a, **k: FakeResponse())

    mol = le.load_ligand(path, "BNZ")
    assert mol.GetNumAtoms() == 6
    assert all(atom.GetIsAromatic() for atom in mol.GetAtoms())
    assert Chem.MolToSmiles(mol) == Chem.CanonSmiles("c1ccccc1")


def test_load_ligand_raises_when_no_template_matches(tmp_path, monkeypatch):
    path = _benzene_pdb(tmp_path)

    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "pdbx_chem_comp_descriptor": [
                    # naphthalene has 10 heavy atoms -- can't match benzene's 6
                    {"program": "CACTVS", "type": "SMILES_CANONICAL", "descriptor": "c1ccc2ccccc2c1"},
                ]
            }

    monkeypatch.setattr(le.requests, "get", lambda *a, **k: FakeResponse())
    with pytest.raises(ValueError):
        le.load_ligand(path, "BNZ")


def test_qed_and_molecular_weight():
    mol = Chem.MolFromSmiles("c1ccccc1O")  # phenol
    assert 90 < le.molecular_weight(mol) < 96
    assert 0 < le.qed(mol) < 1
