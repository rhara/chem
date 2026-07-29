import py3Dmol

from chem.protein import SOLVENT_AND_IONS
from chem.view3d.render import (
    _build_view,
    _caption,
    _chain_ids,
    _ligand_resnames,
    _resolution,
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


def test_caption_includes_all_fields():
    pdb_text = "REMARK   2 RESOLUTION.    2.00 ANGSTROMS.\n" + _atom_line(
        1, "CA", "ALA", "H", 1, 0.0, 0.0, 0.0
    )
    caption = _caption("data/1ABC.pdb", pdb_text, ["LIG"])
    assert "PDB ID: 1ABC" in caption
    assert "Chain: H" in caption
    assert "Ligand: LIG" in caption
    assert "Resolution: 2.00 Å" in caption


def test_caption_defaults_when_no_chain_or_ligand():
    pdb_text = "junk with no ATOM records"
    caption = _caption("data/1ABC.pdb", pdb_text, [])
    assert "Chain: N/A" in caption
    assert "Ligand: none" in caption
    assert "Resolution: N/A" in caption


def test_build_view_has_ligand_style(tmp_path):
    pdb_text = "\n".join(
        [
            _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0, record="HETATM"),
            _atom_line(3, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
        ]
    )
    ligand_resnames = _ligand_resnames(pdb_text, SOLVENT_AND_IONS)

    view = _build_view(pdb_text, ligand_resnames, 600, 500)
    assert isinstance(view, py3Dmol.view)
    assert '"resn": ["LIG"]' in view.startjs
    assert '"color": "magenta"' in view.startjs
    assert "HOH" not in view.startjs.split("addStyle")[1]


def test_build_view_no_ligand_style_when_only_solvent():
    pdb_text = "\n".join(
        [
            _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
            _atom_line(2, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
        ]
    )
    ligand_resnames = _ligand_resnames(pdb_text, SOLVENT_AND_IONS)

    view = _build_view(pdb_text, ligand_resnames, 600, 500)
    assert "stick" not in view.startjs


def test_render_protein_shows_view_then_caption_and_returns_none(tmp_path, monkeypatch):
    lines = [
        "REMARK   2 RESOLUTION.    1.50 ANGSTROMS.",
        _atom_line(1, "CA", "ALA", "H", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "C1", "LIG", "H", 200, 10.0, 10.0, 10.0, record="HETATM"),
    ]
    path = tmp_path / "1ABC.pdb"
    _write_pdb(path, lines)

    events = []
    monkeypatch.setattr(py3Dmol.view, "show", lambda self: events.append("view"))
    monkeypatch.setattr(
        "chem.view3d.render.display", lambda obj: events.append(("caption", obj.data))
    )

    result = render_protein(str(path))

    assert result is None
    assert events[0] == "view"
    assert events[1][0] == "caption"
    assert "PDB ID: 1ABC" in events[1][1]
    assert "Chain: H" in events[1][1]
    assert "Ligand: LIG" in events[1][1]
    assert "Resolution: 1.50 Å" in events[1][1]
