import glob
import os
import re
import shutil
import subprocess
import sys
import tempfile

import numpy as np
from Bio.PDB import MMCIFParser, PDBParser
from rdkit import Chem

from ..verbosity import is_quiet, logged

# Common crystallization solvent/ion/additive HET codes to exclude when
# auto-detecting "the" ligand from a structure's HETATM records. Not exhaustive.
SOLVENT_AND_IONS = frozenset(
    {
        "HOH", "WAT", "DOD",
        "NA", "CL", "K", "MG", "CA", "ZN", "MN", "FE", "FE2", "CO", "NI", "CU", "CD",
        "LI", "RB", "CS", "BR", "IOD", "NH4",
        "SO4", "PO4", "GOL", "EDO", "PEG", "PG4", "1PE", "P6G", "MPD", "FMT", "ACT",
        "DMS", "TRS", "BME", "EPE", "HEPES", "IPA", "UNX", "UNL",
    }
)

_POCKET_HEADER_RE = re.compile(r"^Pocket\s+(\d+)\s*:\s*$")


def _load_structure(path):
    ext = os.path.splitext(path)[1].lower()
    parser = MMCIFParser(QUIET=True) if ext in (".cif", ".mmcif") else PDBParser(QUIET=True)
    return parser.get_structure(os.path.splitext(os.path.basename(path))[0], path)


def _hetero_residues(structure_path):
    """Every non-water heteroatom residue (ligands, ions, crystallization additives)
    in a structure file.
    """
    structure = _load_structure(structure_path)
    residues = []
    for chain in next(structure.get_models()):
        for res in chain:
            hetfield = res.id[0]
            if hetfield in (" ", "W"):
                continue
            residues.append(res)
    return residues


def _residue_atom_coords(residue):
    return np.array([a.coord for a in residue.get_atoms()])


def _residue_label(residue):
    return f"{residue.get_resname()} {residue.parent.id}{residue.id[1]}"


def _auto_detect_ligand(structure_path):
    candidates = [r for r in _hetero_residues(structure_path) if r.get_resname().strip() not in SOLVENT_AND_IONS]
    if not candidates:
        raise ValueError(
            f"no non-solvent/ion HETATM group found in {structure_path}; pass ligand explicitly"
        )
    best = max(candidates, key=lambda r: len(list(r.get_atoms())))
    return _residue_atom_coords(best), _residue_label(best)


def _explicit_code_ligand(structure_path, code):
    code = code.upper()
    matches = [r for r in _hetero_residues(structure_path) if r.get_resname().strip() == code]
    if not matches:
        raise ValueError(f"no HETATM group with code '{code}' found in {structure_path}")
    best = matches[0]
    return _residue_atom_coords(best), _residue_label(best)


def _load_external_ligand_coords(path):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".pdb":
        mol = Chem.MolFromPDBFile(path, sanitize=False, removeHs=False)
    elif ext in (".sdf", ".mol"):
        mol = Chem.MolFromMolFile(path, sanitize=False, removeHs=False)
    elif ext == ".mol2":
        mol = Chem.MolFromMol2File(path, sanitize=False, removeHs=False)
    else:
        raise ValueError(f"unsupported ligand file extension: {ext}")
    if mol is None or mol.GetNumConformers() == 0:
        raise ValueError(f"could not read 3D coordinates from ligand file {path}")
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return coords, os.path.basename(path)


def _resolve_ligand_atoms(structure_path, ligand):
    if ligand is None:
        return _auto_detect_ligand(structure_path)
    if os.path.isfile(ligand):
        return _load_external_ligand_coords(ligand)
    if re.match(r"^[A-Za-z0-9]{1,3}$", ligand):
        return _explicit_code_ligand(structure_path, ligand)
    raise ValueError(f"could not interpret ligand={ligand!r} as a file path or a 1-3 character PDB HET code")


def _run_fpocket(structure, workdir):
    local_path = os.path.join(workdir, os.path.basename(structure))
    if os.path.abspath(local_path) != os.path.abspath(structure):
        shutil.copy(structure, local_path)
    try:
        subprocess.run(
            ["fpocket", "-f", os.path.basename(local_path)],
            cwd=workdir,
            check=True,
            capture_output=True,
            text=True,
        )
    except FileNotFoundError as e:
        raise RuntimeError(
            "fpocket executable not found; install it via conda-forge, "
            "e.g. `mamba install -n chem -c conda-forge fpocket`"
        ) from e
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"fpocket failed on {structure}:\n{e.stderr}") from e

    stem = os.path.splitext(os.path.basename(local_path))[0]
    return os.path.join(workdir, f"{stem}_out"), stem


def _pick_best_pocket(pockets_dir, ligand_atoms):
    """The fpocket pocket whose lining atoms have the smallest minimum distance to
    ligand_atoms. Returns (pocket_id, path to its pocketN_atm.pdb).
    """
    parser = PDBParser(QUIET=True)
    best_id, best_path, best_dist = None, None, None
    for atm_path in sorted(glob.glob(os.path.join(pockets_dir, "pocket*_atm.pdb"))):
        pocket_id = int(re.search(r"pocket(\d+)_atm\.pdb$", atm_path).group(1))
        structure = parser.get_structure("pocket", atm_path)
        atoms = np.array([a.coord for a in structure.get_atoms()])
        dist = np.linalg.norm(atoms[:, None, :] - ligand_atoms[None, :, :], axis=-1).min()
        if best_dist is None or dist < best_dist:
            best_id, best_path, best_dist = pocket_id, atm_path, dist
    if best_id is None:
        raise ValueError("fpocket did not find any pockets")
    return best_id, best_path


def _parse_pocket_residues(atm_path):
    """Unique residues in a pocketN_atm.pdb file, in first-seen order.

    Includes the PDB insertion code: chymotrypsin-numbered serine proteases (e.g.
    trypsin, thrombin, factor Xa) commonly have insertion-code residues like 60A-60I
    sharing a resnum with plain residue 60, which would otherwise collide.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("pocket", atm_path)
    seen = set()
    residues = []
    for chain in next(structure.get_models()):
        for res in chain:
            icode = res.id[2].strip()
            key = (chain.id, res.id[1], icode, res.get_resname())
            if key in seen:
                continue
            seen.add(key)
            residues.append(
                {"chain": chain.id, "resnum": res.id[1], "icode": icode, "resname": res.get_resname()}
            )
    return residues


def _parse_info_txt(path):
    """Parse fpocket's {stem}_info.txt into {pocket_id: {field_name: value}}."""
    pockets = {}
    current = None
    with open(path) as f:
        for line in f:
            header = _POCKET_HEADER_RE.match(line.strip())
            if header:
                current = int(header.group(1))
                pockets[current] = {}
                continue
            if current is None or ":" not in line:
                continue
            key, _, value = line.partition(":")
            key, value = key.strip(), value.strip()
            if not key or not value:
                continue
            try:
                pockets[current][key] = float(value)
            except ValueError:
                pockets[current][key] = value
    return pockets


@logged
def find_pocket(structure, ligand=None, outdir=None):
    """Identify the fpocket pocket nearest a ligand, and its lining residues.

    Runs fpocket on `structure` (a PDB file -- fpocket requires legacy PDB format;
    the output of chem.protein.align is a natural fit), resolves the ligand's 3D
    coordinates, and picks the fpocket pocket whose lining atoms have the smallest
    minimum distance to those coordinates.

    structure: path to a PDB file.
    ligand: how to locate the ligand.
        - None (default): auto-detect the largest non-solvent/ion HETATM group in
          `structure` (e.g. a co-crystallized inhibitor from an RCSB download).
        - a 1-3 character PDB HET code (e.g. "STI"): use that HETATM group in
          `structure`, disambiguating when several ligand-like groups are present.
        - a path to an external ligand file (.pdb/.sdf/.mol/.mol2, e.g. a docking
          pose) with 3D coordinates already in the same frame as `structure`, for
          structures that don't contain the ligand themselves.
        Auto-detect is a size heuristic, not a drug-likeness check: a glycosylation
        sugar (e.g. NAG) or a modified residue on a bound peptide (e.g. TYS) can
        have more atoms than a small-molecule inhibitor and win instead. Pass an
        explicit HET code when a structure has such larger non-solvent HETATM
        groups alongside the intended ligand.
    outdir: optional directory to keep fpocket's full raw output (pockets/,
        *_info.txt, ...) in; if None, a temporary directory is used and discarded.

    Returns a dict:
        pocket_id: fpocket's pocket number
        score / druggability_score / volume: convenience fields pulled from
            fpocket's info file for the selected pocket
        residues: list of {"chain", "resnum", "icode", "resname"} lining the pocket
            (icode is "" when the residue has no PDB insertion code)
        info: the full raw fpocket score dict for the selected pocket
    """
    cleanup = None
    if outdir is not None:
        os.makedirs(outdir, exist_ok=True)
        workdir = outdir
    else:
        cleanup = tempfile.TemporaryDirectory()
        workdir = cleanup.name

    try:
        out_dir, stem = _run_fpocket(structure, workdir)
        ligand_atoms, ligand_label = _resolve_ligand_atoms(structure, ligand)

        pocket_id, atm_path = _pick_best_pocket(os.path.join(out_dir, "pockets"), ligand_atoms)
        residues = _parse_pocket_residues(atm_path)
        info = _parse_info_txt(os.path.join(out_dir, f"{stem}_info.txt"))[pocket_id]

        if not is_quiet():
            print(
                f"selected pocket {pocket_id} (score={info.get('Score')}) near ligand "
                f"{ligand_label}: {len(residues)} residues",
                file=sys.stderr,
            )

        return {
            "pocket_id": pocket_id,
            "score": info.get("Score"),
            "druggability_score": info.get("Druggability Score"),
            "volume": info.get("Volume"),
            "residues": residues,
            "info": info,
        }
    finally:
        if cleanup is not None:
            cleanup.cleanup()
