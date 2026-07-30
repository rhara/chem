import pytest

import chem.protein.pocket as pk


def _atom_line(serial, name, resname, chain, resseq, x, y, z, record="ATOM", icode=" "):
    element = name.strip()[0]
    return (
        f"{record:<6}{serial:>5} {name:<4} {resname:>3} {chain:1}{resseq:>4}{icode:1}   "
        f"{x:>8.3f}{y:>8.3f}{z:>8.3f}{1.0:>6.2f}{0.0:>6.2f}          {element:>2}"
    )


def _write_pdb(path, lines):
    path.write_text("\n".join(lines) + "\nEND\n")


@pytest.fixture
def structure_with_ligand(tmp_path):
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0),
        _atom_line(3, "C1", "LIG", "A", 200, 10.0, 10.0, 10.0, record="HETATM"),
        _atom_line(4, "C2", "LIG", "A", 200, 11.0, 10.0, 10.0, record="HETATM"),
        _atom_line(5, "C3", "LIG", "A", 200, 12.0, 10.0, 10.0, record="HETATM"),
        _atom_line(6, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
        _atom_line(7, "NA", "NA", "A", 301, 30.0, 30.0, 30.0, record="HETATM"),
    ]
    path = tmp_path / "structure.pdb"
    _write_pdb(path, lines)
    return str(path)


def test_auto_detect_ligand_excludes_solvent_and_ions(structure_with_ligand):
    atoms, label = pk._auto_detect_ligand(structure_with_ligand)
    assert len(atoms) == 3
    assert label.startswith("LIG")


def test_auto_detect_ligand_raises_when_only_solvent(tmp_path):
    lines = [
        _atom_line(1, "O1", "HOH", "A", 300, 20.0, 20.0, 20.0, record="HETATM"),
    ]
    path = tmp_path / "no_ligand.pdb"
    _write_pdb(path, lines)
    with pytest.raises(ValueError):
        pk._auto_detect_ligand(str(path))


def test_explicit_code_ligand(structure_with_ligand):
    atoms, label = pk._explicit_code_ligand(structure_with_ligand, "lig")
    assert len(atoms) == 3


def test_explicit_code_ligand_not_found(structure_with_ligand):
    with pytest.raises(ValueError):
        pk._explicit_code_ligand(structure_with_ligand, "ZZZ")


def test_resolve_ligand_atoms_dispatches_on_shape(structure_with_ligand):
    atoms, _ = pk._resolve_ligand_atoms(structure_with_ligand, None)
    assert len(atoms) == 3

    atoms, _ = pk._resolve_ligand_atoms(structure_with_ligand, "LIG")
    assert len(atoms) == 3

    with pytest.raises(ValueError):
        pk._resolve_ligand_atoms(structure_with_ligand, "not-a-code-or-file")


def test_pick_best_pocket_chooses_nearest(tmp_path):
    pockets_dir = tmp_path / "pockets"
    pockets_dir.mkdir()
    _write_pdb(
        pockets_dir / "pocket1_atm.pdb",
        [_atom_line(1, "CA", "ALA", "A", 1, 100.0, 100.0, 100.0)],
    )
    _write_pdb(
        pockets_dir / "pocket2_atm.pdb",
        [_atom_line(1, "CA", "ALA", "A", 2, 0.0, 0.0, 0.0)],
    )
    import numpy as np

    ligand_atoms = np.array([[0.5, 0.0, 0.0]])
    pocket_id, atm_path = pk._pick_best_pocket(str(pockets_dir), ligand_atoms)
    assert pocket_id == 2
    assert atm_path.endswith("pocket2_atm.pdb")


def test_parse_pocket_residues_dedups(tmp_path):
    lines = [
        _atom_line(1, "N", "ALA", "A", 1, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "ALA", "A", 1, 1.5, 0.0, 0.0),
        _atom_line(3, "N", "GLY", "A", 2, 3.0, 0.0, 0.0),
    ]
    path = tmp_path / "pocket1_atm.pdb"
    _write_pdb(path, lines)
    residues = pk._parse_pocket_residues(str(path))
    assert residues == [
        {"chain": "A", "resnum": 1, "icode": "", "resname": "ALA"},
        {"chain": "A", "resnum": 2, "icode": "", "resname": "GLY"},
    ]


def test_parse_pocket_residues_keeps_insertion_code_variants_distinct(tmp_path):
    # Chymotrypsin-numbered serine proteases (trypsin, thrombin, factor Xa, ...)
    # commonly have insertion-code residues like 60A/60B sharing a resnum with
    # plain residue 60 -- these must not collapse into one entry.
    lines = [
        _atom_line(1, "CA", "LEU", "H", 60, 0.0, 0.0, 0.0),
        _atom_line(2, "CA", "TRP", "H", 60, 1.5, 0.0, 0.0, icode="A"),
        _atom_line(3, "CA", "TYR", "H", 60, 3.0, 0.0, 0.0, icode="B"),
    ]
    path = tmp_path / "pocket1_atm.pdb"
    _write_pdb(path, lines)
    residues = pk._parse_pocket_residues(str(path))
    assert residues == [
        {"chain": "H", "resnum": 60, "icode": "", "resname": "LEU"},
        {"chain": "H", "resnum": 60, "icode": "A", "resname": "TRP"},
        {"chain": "H", "resnum": 60, "icode": "B", "resname": "TYR"},
    ]


def test_parse_info_txt(tmp_path):
    content = """Pocket 1 :
\tScore : \t0.553
\tDruggability Score : \t0.027
\tVolume : \t237.521

Pocket 2 :
\tScore : \t0.254
\tDruggability Score : \t0.003
\tVolume : \t407.329
"""
    path = tmp_path / "info.txt"
    path.write_text(content)
    pockets = pk._parse_info_txt(str(path))
    assert pockets[1]["Score"] == 0.553
    assert pockets[1]["Druggability Score"] == 0.027
    assert pockets[2]["Volume"] == 407.329


def test_pocket_atm_paths_maps_id_to_path(tmp_path):
    pockets_dir = tmp_path / "pockets"
    pockets_dir.mkdir()
    _write_pdb(pockets_dir / "pocket1_atm.pdb", [_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)])
    _write_pdb(pockets_dir / "pocket12_atm.pdb", [_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)])
    paths = pk._pocket_atm_paths(str(pockets_dir))
    assert set(paths) == {1, 12}
    assert paths[1].endswith("pocket1_atm.pdb")
    assert paths[12].endswith("pocket12_atm.pdb")


def test_pocket_result_shape_no_vert_file(tmp_path):
    path = tmp_path / "pocket1_atm.pdb"
    _write_pdb(path, [_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)])
    info = {"Score": 0.5, "Druggability Score": 0.6, "Volume": 200.0}
    result = pk._pocket_result(1, str(path), info)
    assert result == {
        "pocket_id": 1,
        "score": 0.5,
        "druggability_score": 0.6,
        "volume": 200.0,
        "residues": [{"chain": "A", "resnum": 1, "icode": "", "resname": "ALA"}],
        "spheres": [],
        "info": info,
    }


def _write_pqr(path, spheres):
    lines = [
        f"ATOM  {i + 1:>5}    C STP     1    {s['x']:>8.3f}{s['y']:>8.3f}{s['z']:>8.3f}    0.00    {s['radius']:.2f}"
        for i, s in enumerate(spheres)
    ]
    path.write_text("\n".join(lines) + "\nEND\n")


def test_parse_pocket_spheres(tmp_path):
    path = tmp_path / "pocket1_vert.pqr"
    _write_pqr(
        path,
        [
            {"x": 1.0, "y": 2.0, "z": 3.0, "radius": 3.5},
            {"x": -1.5, "y": 0.0, "z": 4.25, "radius": 3.6},
        ],
    )
    spheres = pk._parse_pocket_spheres(str(path))
    assert spheres == [
        {"x": 1.0, "y": 2.0, "z": 3.0, "radius": 3.5},
        {"x": -1.5, "y": 0.0, "z": 4.25, "radius": 3.6},
    ]


def test_pocket_result_includes_spheres_from_sibling_vert_file(tmp_path):
    atm_path = tmp_path / "pocket1_atm.pdb"
    _write_pdb(atm_path, [_atom_line(1, "CA", "ALA", "A", 1, 0.0, 0.0, 0.0)])
    _write_pqr(tmp_path / "pocket1_vert.pqr", [{"x": 1.0, "y": 2.0, "z": 3.0, "radius": 3.5}])
    result = pk._pocket_result(1, str(atm_path), {"Score": 0.1})
    assert result["spheres"] == [{"x": 1.0, "y": 2.0, "z": 3.0, "radius": 3.5}]


@pytest.fixture
def fpocket_output(tmp_path):
    """A fake fpocket output directory (pockets/pocketN_atm.pdb + N_out/stem_info.txt)
    with three pockets of differing (and one missing) druggability score."""
    out_dir = tmp_path / "structure_out"
    pockets_dir = out_dir / "pockets"
    pockets_dir.mkdir(parents=True)
    for pocket_id in (1, 2, 3):
        _write_pdb(
            pockets_dir / f"pocket{pocket_id}_atm.pdb",
            [_atom_line(1, "CA", "ALA", "A", pocket_id, 0.0, 0.0, 0.0)],
        )
    (out_dir / "structure_info.txt").write_text(
        "Pocket 1 :\n\tScore : \t0.1\n\tDruggability Score : \t0.2\n\tVolume : \t100.0\n\n"
        "Pocket 2 :\n\tScore : \t0.3\n\tDruggability Score : \t0.9\n\tVolume : \t300.0\n\n"
        # No Druggability Score line -- fpocket sometimes omits it.
        "Pocket 3 :\n\tScore : \t0.05\n\tVolume : \t50.0\n"
    )
    return str(out_dir)


def test_list_pockets_returns_every_pocket_sorted_by_druggability(tmp_path, fpocket_output, monkeypatch):
    monkeypatch.setattr(pk, "_run_fpocket", lambda structure, workdir: (fpocket_output, "structure"))
    dummy = tmp_path / "structure.pdb"
    dummy.write_text("")

    results = pk.list_pockets(str(dummy))

    assert [r["pocket_id"] for r in results] == [2, 1, 3]
    assert results[0]["druggability_score"] == 0.9
    assert results[-1]["druggability_score"] is None


def test_run_fpocket_missing_binary_raises_helpful_error(tmp_path, monkeypatch):
    monkeypatch.setattr(pk.shutil, "which", lambda *_: None)

    def fake_run(*args, **kwargs):
        raise FileNotFoundError("fpocket")

    monkeypatch.setattr(pk.subprocess, "run", fake_run)
    dummy = tmp_path / "structure.pdb"
    dummy.write_text("END\n")
    with pytest.raises(RuntimeError, match="fpocket executable not found"):
        pk._run_fpocket(str(dummy), str(tmp_path))
