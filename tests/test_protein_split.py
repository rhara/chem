import os

import pytest
from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem import AllChem

from chem.protein import split


def _atom_line(serial, name, resname, chain, resnum, x, y, z, record="ATOM"):
    element = name.strip()[0]
    return (
        f"{record:<6}{serial:>5} {name:<4} {resname:>3} {chain:1}{resnum:>4}    "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element:>2}"
    )


def _write_pdb(path, lines):
    path.write_text("\n".join(lines) + "\nEND\n")


def _benzene_hetatm_lines(resname, chain, resnum, offset, start_serial):
    """A benzene ring with real (embedded/MMFF-optimized) 3D coordinates but no
    bond-order info once round-tripped through PDB, like a real HETATM group --
    same trick as tests/test_ligand_extract.py's benzene fixture.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles("c1ccccc1"))
    AllChem.EmbedMolecule(mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    lines = []
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conf.GetAtomPosition(i)
        x, y, z = pos.x + offset[0], pos.y + offset[1], pos.z + offset[2]
        name = f"{atom.GetSymbol()}{i}"
        lines.append(_atom_line(start_serial + i, name, resname, chain, resnum, x, y, z, record="HETATM"))
    return lines


class _FakeResponse:
    def __init__(self, descriptors):
        self._descriptors = descriptors

    def raise_for_status(self):
        pass

    def json(self):
        return {"pdbx_chem_comp_descriptor": self._descriptors}


@pytest.fixture(autouse=True)
def mock_chemcomp_requests(monkeypatch):
    """split() calls load_ligand for every non-water HETATM residue -- mock the PDB
    Chemical Component Dictionary lookup for every test so none of them hit the
    network. "LIG"/"OTH" resolve to a benzene template; any other code (e.g. "NA")
    gets no descriptors, so load_ligand raises and split() falls back to writing
    raw (bond-order-unrestored) connectivity for it instead.
    """
    import chem.ligand.extract as le

    descriptors_by_code = {
        "LIG": [{"program": "CACTVS", "type": "SMILES_CANONICAL", "descriptor": "c1ccccc1"}],
        "OTH": [{"program": "CACTVS", "type": "SMILES_CANONICAL", "descriptor": "c1ccccc1"}],
    }

    def fake_get(url, timeout=30):
        code = url.rsplit("/", 1)[-1]
        return _FakeResponse(descriptors_by_code.get(code, []))

    monkeypatch.setattr(le.requests, "get", fake_get)


@pytest.fixture
def structure_path(tmp_path):
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0),
        *_benzene_hetatm_lines("LIG", "A", 200, (10.0, 10.0, 10.0), 3),
        _atom_line(20, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
        _atom_line(21, "NA", "NA", "A", 301, 30.0, 30.0, 30.0, record="HETATM"),
        _atom_line(30, "N", "ALA", "B", 1, 0.0, 0.0, 50.0),
        _atom_line(31, "CA", "ALA", "B", 1, 1.5, 0.0, 50.0),
        *_benzene_hetatm_lines("OTH", "B", 200, (10.0, 10.0, 60.0), 32),
        _atom_line(50, "O1", "HOH", "B", 300, 20.0, 20.0, 70.0, record="HETATM"),
    ]
    path = tmp_path / "structure.pdb"
    _write_pdb(path, lines)
    return str(path)


def _residue_hetero_codes(pdb_path, chain_id=None):
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("x", pdb_path)
    codes = set()
    for chain in next(structure.get_models()):
        if chain_id is not None and chain.id != chain_id:
            continue
        for res in chain:
            if res.id[0] != " ":
                codes.add(res.get_resname().strip())
    return codes


def test_split_writes_per_chain_protein_pdbs_by_default(structure_path, tmp_path):
    outdir = str(tmp_path / "out")
    result = split(structure_path, outdir=outdir)

    assert set(result["protein"].keys()) == {"A", "B"}
    assert result["protein"]["A"] == os.path.join(outdir, "structure_protein_A.pdb")
    assert result["protein"]["B"] == os.path.join(outdir, "structure_protein_B.pdb")

    parser = PDBParser(QUIET=True)
    chain_a_only = next(parser.get_structure("a", result["protein"]["A"]).get_models())
    assert [c.id for c in chain_a_only] == ["A"]
    chain_b_only = next(parser.get_structure("b", result["protein"]["B"]).get_models())
    assert [c.id for c in chain_b_only] == ["B"]

    codes_a = _residue_hetero_codes(result["protein"]["A"], "A")
    codes_b = _residue_hetero_codes(result["protein"]["B"], "B")
    assert codes_a == {"HOH"}  # LIG and NA (ligands) dropped, water kept
    assert codes_b == {"HOH"}  # OTH dropped


def test_split_remove_water_strips_water_from_protein_pdb(structure_path, tmp_path):
    outdir = str(tmp_path / "out")
    kept = split(structure_path, outdir=str(tmp_path / "kept"))
    removed = split(structure_path, remove_water=True, outdir=outdir)

    assert _residue_hetero_codes(kept["protein"]["A"], "A") == {"HOH"}  # default: water kept
    assert _residue_hetero_codes(removed["protein"]["A"], "A") == set()
    assert _residue_hetero_codes(removed["protein"]["B"], "B") == set()

    # water is never written as a ligand SDF either way, remove_water only affects the protein PDB
    assert "HOH" not in {lig["code"] for lig in removed["ligands"]}


def test_split_creates_outdir_if_missing(structure_path, tmp_path):
    outdir = str(tmp_path / "does" / "not" / "exist" / "yet")
    result = split(structure_path, outdir=outdir)
    assert os.path.exists(result["protein"]["A"])


def test_split_all_chains_writes_single_protein_pdb(structure_path, tmp_path):
    outdir = str(tmp_path / "out")
    result = split(structure_path, all_chains=True, outdir=outdir)

    assert result["protein"] == os.path.join(outdir, "structure_protein.pdb")
    assert os.path.exists(result["protein"])

    codes_a = _residue_hetero_codes(result["protein"], "A")
    codes_b = _residue_hetero_codes(result["protein"], "B")
    assert codes_a == {"HOH"}  # LIG and NA (ligands) dropped, water kept
    assert codes_b == {"HOH"}  # OTH dropped


def test_split_writes_ligand_sdf_with_restored_bond_orders(structure_path, tmp_path):
    outdir = str(tmp_path / "out")
    result = split(structure_path, outdir=outdir)

    ligands_by_code = {lig["code"]: lig for lig in result["ligands"]}
    # every non-water HETATM instance gets an SDF, including "NA" (see the next test) --
    # LIG/OTH are just the two that get proper bond-order restoration.
    assert set(ligands_by_code) == {"LIG", "OTH", "NA"}

    lig_path = ligands_by_code["LIG"]["path"]
    assert lig_path == os.path.join(outdir, "structure_ligand_LIG_A200.sdf")
    assert os.path.exists(lig_path)
    assert ligands_by_code["LIG"]["bond_orders_restored"] is True
    mol = next(Chem.SDMolSupplier(lig_path))
    assert mol.GetNumAtoms() == 6
    assert all(atom.GetIsAromatic() for atom in mol.GetAtoms())

    # the protein PDB has NA excluded too, even though its SDF wasn't bond-order
    # restored -- it's still a non-water ligand candidate that split() strips out
    # of the protein, same as any other.
    assert "NA" not in _residue_hetero_codes(result["protein"]["A"], "A")


def test_split_writes_raw_sdf_when_no_bond_order_template_matches(structure_path, tmp_path):
    outdir = str(tmp_path / "out")
    result = split(structure_path, outdir=outdir)

    na_entry = next(lig for lig in result["ligands"] if lig["code"] == "NA")
    assert na_entry["bond_orders_restored"] is False
    assert na_entry["path"] == os.path.join(outdir, "structure_ligand_NA_A301.sdf")
    assert os.path.exists(na_entry["path"])

    mol = next(Chem.SDMolSupplier(na_entry["path"], sanitize=False))
    assert mol.GetNumAtoms() == 1
    assert mol.GetConformer().GetAtomPosition(0).x == pytest.approx(30.0)
