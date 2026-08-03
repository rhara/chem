import os
import tempfile

import numpy as np
import requests
from Bio.PDB import PDBIO, Select
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, QED
from rdkit.Geometry import Point3D

from ..protein.pocket import SOLVENT_AND_IONS, _hetero_residues
from ..verbosity import logged

RCSB_CHEMCOMP_API = "https://data.rcsb.org/rest/v1/core/chemcomp"

# Which chem_comp descriptor to use as the bond-order template, in priority
# order -- OpenEye's canonical form is typically the cleanest, CACTVS is a
# solid fallback. Stereochemistry in the template doesn't matter here (only
# the bond graph is used), so plain "SMILES" entries are fine too.
_DESCRIPTOR_PRIORITY = [
    ("OpenEye OEToolkits", "SMILES_CANONICAL"),
    ("OpenEye OEToolkits", "SMILES"),
    ("CACTVS", "SMILES_CANONICAL"),
    ("CACTVS", "SMILES"),
]


class _ResidueSelect(Select):
    def __init__(self, residue):
        self._residue = residue

    def accept_residue(self, residue):
        return residue is self._residue


def _instance_info(residue):
    _, resnum, icode = residue.id
    return {
        "code": residue.get_resname().strip(),
        "chain": residue.get_parent().id,
        "resnum": resnum,
        "icode": icode.strip(),
    }


def list_ligand_instances(structure_path, exclude=SOLVENT_AND_IONS):
    """Every non-excluded HETATM residue instance in a structure file, one
    entry per physical occurrence -- e.g. two copies of the same ligand code
    bound to different chains produce two entries, not one. Each entry is a
    dict: `{"code", "chain", "resnum", "icode"}`, in file order.
    """
    return [
        _instance_info(r)
        for r in _hetero_residues(structure_path)
        if r.get_resname().strip() not in exclude
    ]


def list_ligand_codes(structure_path, exclude=SOLVENT_AND_IONS):
    """Every distinct non-excluded HETATM residue code (3-letter PDB chemical
    component id) in a structure file, e.g. ["S54"]. Multiple copies of the
    same code (e.g. one per chain) collapse to a single entry here -- use
    `list_ligand_instances` to enumerate every physical occurrence instead.
    """
    codes = {i["code"] for i in list_ligand_instances(structure_path, exclude=exclude)}
    return sorted(codes)


def _pick_ligand_residue(structure_path, code, chain=None, resnum=None, icode=None):
    matches = [r for r in _hetero_residues(structure_path) if r.get_resname().strip() == code.upper()]
    if chain is not None:
        matches = [r for r in matches if r.get_parent().id == chain]
    if resnum is not None:
        matches = [r for r in matches if r.id[1] == resnum]
    if icode is not None:
        matches = [r for r in matches if r.id[2].strip() == icode]
    if not matches:
        raise ValueError(
            f"no HETATM group with code '{code}'"
            f"{f', chain={chain!r}' if chain is not None else ''}"
            f"{f', resnum={resnum!r}' if resnum is not None else ''}"
            f"{f', icode={icode!r}' if icode is not None else ''}"
            f" found in {structure_path}"
        )
    # Several copies of the same ligand (e.g. one per chain) are chemically
    # identical when unresolved, so with no chain/resnum/icode given to pin
    # down a specific one, pick the most complete instance.
    return max(matches, key=lambda r: len(list(r.get_atoms())))


def _residue_to_raw_mol(residue, code):
    # Reuse the residue's own structure/model/chain tree (rather than re-parsing
    # structure_path) since PDBIO.save's Select filter matches by object identity.
    structure = residue.get_parent().get_parent().get_parent()
    writer = PDBIO()
    writer.set_structure(structure)
    fd, tmp_path = tempfile.mkstemp(suffix=".pdb")
    os.close(fd)
    try:
        writer.save(tmp_path, _ResidueSelect(residue))
        # removeHs=True is a no-op unless sanitize=True (which the bond-order-less
        # raw mol can't survive yet), so strip explicit Hs -- when present in the
        # structure at all -- separately.
        mol = Chem.MolFromPDBFile(tmp_path, sanitize=False, removeHs=False)
    finally:
        os.remove(tmp_path)
    if mol is None or mol.GetNumAtoms() == 0:
        raise ValueError(f"RDKit could not parse ligand residue {code}")
    return Chem.RemoveHs(mol, sanitize=False)


def _fetch_template_smiles_candidates(code):
    """Candidate template SMILES for a chemical component, best first."""
    resp = requests.get(f"{RCSB_CHEMCOMP_API}/{code.upper()}", timeout=30)
    resp.raise_for_status()
    descriptors = resp.json().get("pdbx_chem_comp_descriptor", [])
    by_key = {(d["program"], d["type"]): d["descriptor"] for d in descriptors}
    ordered = [by_key[key] for key in _DESCRIPTOR_PRIORITY if key in by_key]
    ordered += [d["descriptor"] for d in descriptors if d["type"] in ("SMILES", "SMILES_CANONICAL")]
    # de-duplicate while preserving order
    seen = set()
    candidates = [s for s in ordered if not (s in seen or seen.add(s))]
    if not candidates:
        raise ValueError(f"no SMILES descriptor found for chemical component '{code}'")
    return candidates


def _assign_bond_orders(mol, code):
    candidates = _fetch_template_smiles_candidates(code)
    tried = 0
    for smiles in candidates:
        template = Chem.MolFromSmiles(smiles)
        if template is None:
            continue
        # Some CCD SMILES spell out a stereo-defining H (e.g. "[H]/N=C(...)")
        # as an explicit atom; strip it back down to the heavy-atom graph that
        # matches what was extracted from the PDB file (removeHs=True).
        template = Chem.RemoveHs(template)
        if template.GetNumAtoms() != mol.GetNumAtoms():
            continue
        tried += 1
        try:
            fixed = AllChem.AssignBondOrdersFromTemplate(template, mol)
            Chem.SanitizeMol(fixed)
            return fixed
        except ValueError:
            continue
    raise ValueError(
        f"could not match a bond-order template to ligand '{code}' "
        f"({len(candidates)} candidate SMILES fetched, {tried} atom-count-compatible, "
        f"all failed)"
    )


@logged
def load_ligand(structure_path, ligand, chain=None, resnum=None, icode=None):
    """Extract a ligand from a structure file as a proper RDKit molecule.

    `ligand` is its 3-letter PDB chemical component code (see
    `list_ligand_codes`/`list_ligand_instances`). When a code occurs more than
    once (e.g. one copy per chain), pass `chain`/`resnum`/`icode` -- as found
    in a `list_ligand_instances` entry -- to pick that exact instance;
    otherwise the most complete matching instance is used.

    The residue's atoms and 3D coordinates come straight from the structure
    file, but PDB format has no bond-order information, so RDKit's initial
    guess is all single bonds with no aromaticity -- this is corrected by
    fetching the PDB Chemical Component Dictionary's ideal SMILES for
    `ligand` and using it as a bond-order template
    (`rdkit.Chem.AllChem.AssignBondOrdersFromTemplate`). Raises if the
    template doesn't match the extracted atoms (e.g. a covalently modified or
    incomplete residue).
    """
    residue = _pick_ligand_residue(structure_path, ligand, chain=chain, resnum=resnum, icode=icode)
    raw_mol = _residue_to_raw_mol(residue, ligand)
    return _assign_bond_orders(raw_mol, ligand)


@logged
def apply_transform(sdf_path, rotation, translation, outpath):
    """Apply a rotation+translation -- as returned by
    chem.protein.compute_transform -- to an SDF file's 3D coordinates and
    write the result to `outpath`.

    Complements chem.protein.apply_transform (for protein PDB chains): a
    transform computed from one chain of a complex (e.g. via
    chem.protein.compute_transform on the kinase chain that chem.protein.split
    wrote out) can be reapplied here to that same entry's ligand SDF(s), so
    the ligand ends up in the same coordinate frame as the aligned protein --
    reconstituting the bound complex without needing a sequence to align the
    ligand on.

    sdf_path: SDF file to transform (e.g. a chem.protein.split ligand path).
        Read with sanitize=False since some SDFs written by split() (those
        with bond_orders_restored=False) can't survive sanitization.
    rotation, translation: as returned by chem.protein.compute_transform's
        "rotation"/"translation" -- applied in Biopython's convention,
        `new_coord = old_coord @ rotation + translation`, matching
        chem.protein.apply_transform so the same pair of matrices moves both
        the protein and the ligand consistently.
    outpath: destination SDF file path; parent directory created if missing.

    Returns outpath.
    """
    mol = next(Chem.SDMolSupplier(sdf_path, sanitize=False))
    rot = np.asarray(rotation)
    tran = np.asarray(translation)

    conf = mol.GetConformer()
    coords = conf.GetPositions()
    new_coords = coords @ rot + tran
    for i, (x, y, z) in enumerate(new_coords):
        conf.SetAtomPosition(i, Point3D(x, y, z))

    outdir = os.path.dirname(outpath)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    writer = Chem.SDWriter(outpath)
    writer.write(mol)
    writer.close()
    return outpath


def qed(mol):
    """Quantitative Estimate of Drug-likeness (Bickerton et al., 2012), 0-1."""
    return QED.qed(mol)


def molecular_weight(mol):
    """Average molecular weight (g/mol)."""
    return Descriptors.MolWt(mol)
