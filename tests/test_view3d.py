import py3Dmol

from chem.protein import SOLVENT_AND_IONS
from chem.view3d.render import _ligand_resnames, render_protein


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


def test_render_protein_returns_view_with_ligand_style(tmp_path):
    lines = [
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0, record="HETATM"),
        _atom_line(3, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
    ]
    path = tmp_path / "structure.pdb"
    _write_pdb(path, lines)

    view = render_protein(str(path))
    assert isinstance(view, py3Dmol.view)
    assert '"resn": ["LIG"]' in view.startjs
    assert "HOH" not in view.startjs.split("addStyle")[1]


def test_render_protein_no_ligand_style_when_only_solvent(tmp_path):
    lines = [
        _atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
    ]
    path = tmp_path / "structure.pdb"
    _write_pdb(path, lines)

    view = render_protein(str(path))
    assert "stick" not in view.startjs
