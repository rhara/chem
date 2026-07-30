import math
import re

import py3Dmol
import pytest
from rdkit import Chem

from chem.protein import SOLVENT_AND_IONS
from chem.view3d.render import (
    _build_view,
    _caption_lines,
    _chain_ids,
    _instance_block,
    _ligand_instances,
    _ligand_molblock,
    _ligand_resnames,
    _resolution,
    _template_mol,
    render_protein,
)


def _atom_line(serial, name, resname, chain, resseq, x, y, z, record="ATOM"):
    element = name.strip()[0]
    return (
        f"{record:<6}{serial:>5} {name:<4} {resname:>3} {chain:1}{resseq:>4}   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element:>2}"
    )


def _write_pdb(path, lines):
    path.write_text("\n".join(lines) + "\nEND\n")


def test_ligand_resnames_excludes_given_codes():
    pdb_text = "\n".join(
        [
            _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0, record="HETATM"),
            _atom_line(3, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
            _atom_line(4, "NA", "NA", "A", 301, 30.0, 30.0, 30.0, record="HETATM"),
        ]
    )
    assert _ligand_resnames(pdb_text, SOLVENT_AND_IONS) == ["LIG"]


def test_ligand_resnames_accepts_extra_excludes():
    pdb_text = "\n".join(
        [
            _atom_line(1, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0, record="HETATM"),
            _atom_line(2, "C1", "NAG", "A", 400, 40.0, 40.0, 40.0, record="HETATM"),
        ]
    )
    assert _ligand_resnames(pdb_text, SOLVENT_AND_IONS | {"NAG"}) == ["LIG"]


def test_ligand_resnames_empty_when_no_hetatm():
    pdb_text = _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)
    assert _ligand_resnames(pdb_text, SOLVENT_AND_IONS) == []


def _benzene_ring_lines(resname="LIG", chain="A", resseq=1):
    """Six HETATM lines laid out as a regular hexagon (1.4 Å sides, the
    aromatic C-C bond length) -- enough for RDKit to distance-perceive a
    six-membered ring."""
    lines = []
    for i in range(6):
        angle = math.radians(60 * i)
        x, y = 1.4 * math.cos(angle), 1.4 * math.sin(angle)
        lines.append(
            _atom_line(i + 1, f"C{i + 1}", resname, chain, resseq, x, y, 0.0, record="HETATM")
        )
    return lines


def _benzene_template():
    return Chem.RemoveHs(Chem.MolFromSmiles("c1ccccc1"))


def test_ligand_instances_dedupes_atoms_within_one_instance():
    pdb_text = "\n".join(
        [
            _atom_line(1, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0, record="HETATM"),
            _atom_line(2, "C2", "LIG", "A", 200, 11.0, 10.0, 10.0, record="HETATM"),
        ]
    )
    assert _ligand_instances(pdb_text, ["LIG"]) == [("LIG", "A", "200", "")]


def test_ligand_instances_distinguishes_by_chain_and_resnum_in_order():
    pdb_text = "\n".join(
        [
            _atom_line(1, "C1", "LIG", "B", 200, 0.0, 0.0, 0.0, record="HETATM"),
            _atom_line(2, "C1", "LIG", "A", 200, 1.0, 0.0, 0.0, record="HETATM"),
            _atom_line(3, "C1", "LIG", "B", 201, 2.0, 0.0, 0.0, record="HETATM"),
        ]
    )
    assert _ligand_instances(pdb_text, ["LIG"]) == [
        ("LIG", "B", "200", ""),
        ("LIG", "A", "200", ""),
        ("LIG", "B", "201", ""),
    ]


def test_ligand_instances_ignores_other_resnames():
    pdb_text = "\n".join(
        [
            _atom_line(1, "C1", "LIG", "A", 200, 0.0, 0.0, 0.0, record="HETATM"),
            _atom_line(2, "O1", "HOH", "A", 300, 1.0, 0.0, 0.0, record="HETATM"),
        ]
    )
    assert _ligand_instances(pdb_text, ["LIG"]) == [("LIG", "A", "200", "")]


def test_instance_block_selects_only_matching_atoms():
    pdb_text = "\n".join(
        [
            _atom_line(1, "C1", "LIG", "A", 200, 0.0, 0.0, 0.0, record="HETATM"),
            _atom_line(2, "C1", "LIG", "B", 200, 1.0, 0.0, 0.0, record="HETATM"),
            _atom_line(3, "O1", "HOH", "A", 300, 2.0, 0.0, 0.0, record="HETATM"),
        ]
    )
    block = _instance_block(pdb_text, "LIG", "A", "200", "")
    assert block.count("\n") == 2  # one atom line + trailing END
    assert "C1   LIG A 200" in block
    assert "HOH" not in block
    assert block.endswith("END\n")


def test_template_mol_parses_canonical_smiles(monkeypatch):
    calls = []

    class _FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "pdbx_chem_comp_descriptor": [
                    {"type": "SMILES", "descriptor": "c1ccccc1O"},
                    {"type": "SMILES_CANONICAL", "descriptor": "c1ccc(cc1)O"},
                ]
            }

    def fake_get(url, timeout):
        calls.append(url)
        return _FakeResponse()

    monkeypatch.setattr("chem.view3d.render.requests.get", fake_get)
    _template_mol.cache_clear()
    try:
        mol = _template_mol("PHN")
        assert mol is not None
        assert mol.GetNumAtoms() == 7  # phenol, Hs stripped
        assert calls == ["https://data.rcsb.org/rest/v1/core/chemcomp/PHN"]

        _template_mol("PHN")  # cached -- no second request
        assert calls == ["https://data.rcsb.org/rest/v1/core/chemcomp/PHN"]
    finally:
        _template_mol.cache_clear()


def test_template_mol_none_when_lookup_fails(monkeypatch):
    def fake_get(url, timeout):
        raise ConnectionError("no network")

    monkeypatch.setattr("chem.view3d.render.requests.get", fake_get)
    _template_mol.cache_clear()
    try:
        assert _template_mol("XXX") is None
    finally:
        _template_mol.cache_clear()


def test_ligand_molblock_kekulizes_aromatic_ring(monkeypatch):
    pdb_text = "\n".join(_benzene_ring_lines())
    monkeypatch.setattr("chem.view3d.render._template_mol", lambda resname: _benzene_template())

    molblock = _ligand_molblock(pdb_text, "LIG", "A", "1", "")

    assert molblock is not None
    bond_lines = molblock.splitlines()[10:16]  # counts line, then 6 atom lines, then 6 bond lines
    bond_orders = {int(line.split()[2]) for line in bond_lines}
    assert bond_orders == {1, 2}  # alternating single/double, not all-single


def test_ligand_molblock_none_when_template_lookup_fails(monkeypatch):
    pdb_text = "\n".join(_benzene_ring_lines())
    monkeypatch.setattr("chem.view3d.render._template_mol", lambda resname: None)

    assert _ligand_molblock(pdb_text, "LIG", "A", "1", "") is None


def test_ligand_molblock_none_when_atom_count_mismatches(monkeypatch):
    pdb_text = "\n".join(_benzene_ring_lines())
    # A template with a different atom count than the 6-carbon ring fixture.
    monkeypatch.setattr(
        "chem.view3d.render._template_mol", lambda resname: Chem.MolFromSmiles("CCO")
    )

    assert _ligand_molblock(pdb_text, "LIG", "A", "1", "") is None


def test_chain_ids_collects_distinct_atom_chains():
    pdb_text = "\n".join(
        [
            _atom_line(1, "CA", "ALA", "H", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "CA", "GLY", "L", 1, 1.0, 0.0, 0.0),
            _atom_line(3, "C1", "LIG", "H", 200, 10.0, 10.0, 10.0, record="HETATM"),
        ]
    )
    assert _chain_ids(pdb_text) == ["H", "L"]


def test_chain_ids_empty_when_no_atom_records():
    pdb_text = _atom_line(1, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0, record="HETATM")
    assert _chain_ids(pdb_text) == []


def test_resolution_parses_remark_2():
    pdb_text = "REMARK   2 RESOLUTION.    1.16 ANGSTROMS.\n" + _atom_line(
        1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0
    )
    assert _resolution(pdb_text) == "1.16 Å"


def test_resolution_na_when_remark_missing():
    pdb_text = _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)
    assert _resolution(pdb_text) == "N/A"


def test_caption_lines_includes_all_fields():
    pdb_text = "REMARK   2 RESOLUTION.    2.00 ANGSTROMS.\n" + _atom_line(
        1, "CA", "ALA", "H", 1, 0.0, 0.0, 0.0
    )
    lines = _caption_lines("data/1ABC.pdb", pdb_text, ["LIG"])
    assert lines == ["PDB ID: 1ABC", "Chain: H", "Ligand: LIG", "Resolution: 2.00 Å"]


def test_caption_lines_defaults_when_no_chain_or_ligand():
    pdb_text = "junk with no ATOM records"
    lines = _caption_lines("data/1ABC.pdb", pdb_text, [])
    assert lines == ["PDB ID: 1ABC", "Chain: N/A", "Ligand: none", "Resolution: N/A"]


def test_build_view_falls_back_to_plain_stick_when_no_template(tmp_path, monkeypatch):
    # No network / unknown HET code -- _ligand_molblock returns None for every
    # instance, so _build_view must fall back to a plain elemental-colored stick.
    monkeypatch.setattr("chem.view3d.render._ligand_molblock", lambda *args: None)

    pdb_text = "\n".join(
        [
            _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0, record="HETATM"),
            _atom_line(3, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
        ]
    )
    ligand_resnames = _ligand_resnames(pdb_text, SOLVENT_AND_IONS)

    view = _build_view(pdb_text, ligand_resnames, 600, 500, "spectrum", (50, 90))
    assert isinstance(view, py3Dmol.view)
    assert '"resn": "LIG"' in view.startjs
    assert '"resi": 200' in view.startjs
    assert '"colorscheme": "magentaCarbon"' in view.startjs
    assert "addModel" not in view.startjs.split("addModel", 1)[1]  # only the protein's own addModel
    assert "HOH" not in view.startjs.split("addStyle")[1]


def test_build_view_adds_a_bond_order_model_when_template_matches(monkeypatch):
    # _ligand_molblock succeeding means _build_view must addModel() the
    # returned MolBlock (bond orders / aromatic rings) instead of falling
    # back to a plain distance-only stick selector.
    monkeypatch.setattr(
        "chem.view3d.render._ligand_molblock",
        lambda pdb_text, resname, chain, resnum, icode: "FAKE MOLBLOCK",
    )

    pdb_text = "\n".join(
        [
            _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0, record="HETATM"),
        ]
    )
    ligand_resnames = _ligand_resnames(pdb_text, SOLVENT_AND_IONS)

    view = _build_view(pdb_text, ligand_resnames, 600, 500, "spectrum", (50, 90))
    assert 'addModel("FAKE MOLBLOCK","mol")' in view.startjs
    assert '"colorscheme": "magentaCarbon"' in view.startjs
    assert '"model": -1' in view.startjs
    assert '"resn": "LIG"' not in view.startjs  # no distance-only fallback selector


def test_build_view_no_ligand_style_when_only_solvent():
    pdb_text = "\n".join(
        [
            _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
        ]
    )
    ligand_resnames = _ligand_resnames(pdb_text, SOLVENT_AND_IONS)

    view = _build_view(pdb_text, ligand_resnames, 600, 500, "spectrum", (50, 90))
    assert "stick" not in view.startjs


def test_build_view_spectrum_coloring():
    pdb_text = _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)
    view = _build_view(pdb_text, [], 600, 500, "spectrum", (50, 90))
    assert '"color": "spectrum"' in view.startjs
    assert '"colorscheme": "roygb"' in view.startjs
    assert '"prop": "b"' not in view.startjs


def test_build_view_bfactor_coloring_uses_given_range():
    pdb_text = _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)
    view = _build_view(pdb_text, [], 600, 500, "bfactor", (30, 95))
    assert '"prop": "b"' in view.startjs
    assert '"gradient": "roygb"' in view.startjs
    assert '"min": 30' in view.startjs
    assert '"max": 95' in view.startjs


def test_render_protein_rejects_bad_coloring(tmp_path):
    path = tmp_path / "1ABC.pdb"
    _write_pdb(path, [_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)])
    with pytest.raises(ValueError):
        render_protein(str(path), coloring="rainbow")


def test_render_protein_frames_and_returns_none(tmp_path, monkeypatch):
    lines = [
        "REMARK   2 RESOLUTION.    1.50 ANGSTROMS.",
        _atom_line(1, "CA", "ALA", "H", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "C1", "LIG", "H", 200, 10.0, 10.0, 10.0, record="HETATM"),
    ]
    path = tmp_path / "1ABC.pdb"
    _write_pdb(path, lines)

    events = []
    monkeypatch.setattr(
        "chem.view3d.render.display", lambda obj: events.append(("display", obj.data))
    )
    monkeypatch.setattr(
        py3Dmol.view, "insert", lambda self, containerid: events.append(("insert", containerid))
    )
    # No real network calls to the RCSB Chemical Component Dictionary here --
    # this test is about the frame/caption plumbing, not bond-order lookup.
    monkeypatch.setattr("chem.view3d.render._template_mol", lambda resname: None)

    result = render_protein(str(path), width=600, height=500)

    assert result is None
    # The bordered placeholder + caption must be displayed before the view is
    # insert()ed into it, since insert()'s JS looks up the placeholder by id.
    assert events[0][0] == "display"
    html = events[0][1]
    assert events[1][0] == "insert"

    assert "border:1px solid #ccc" in html
    assert "width:600px" in html
    assert "height:500px" in html
    assert "PDB ID: 1ABC" in html
    assert "Chain: H" in html
    assert "Ligand: LIG" in html
    assert "Resolution: 1.50 Å" in html

    # The placeholder id in the displayed HTML must match what insert() received.
    match = re.search(r'id="(chem-view3d-[0-9a-f]+)"', html)
    assert match is not None
    assert events[1][1] == match.group(1)
