import os
import sys

from Bio.PDB import PDBIO, Select
from rdkit import Chem

from ..verbosity import is_quiet, logged
from .pocket import WATER, _hetero_residues, _load_structure


class _ProteinSelect(Select):
    """Keeps polymer residues and water, drops every other HETATM residue --
    those are written out separately as SDF ligands by `split` instead.
    Optionally restricted to a single chain.
    """

    def __init__(self, exclude_residues, chain_id=None):
        self._exclude = exclude_residues
        self._chain_id = chain_id

    def accept_chain(self, chain):
        return self._chain_id is None or chain.id == self._chain_id

    def accept_residue(self, residue):
        return residue not in self._exclude


@logged
def split(structure_path, split_chains=False, outdir="split"):
    """Split a structure file into a ligand-free protein PDB and one SDF file
    per non-water HETATM ligand instance -- e.g. to prep a receptor/ligand
    pair for docking.

    structure_path: path to a PDB/CIF structure file.
    split_chains: if True, write one protein PDB per chain, its filename
        including the chain id, instead of a single PDB with every chain
        together. Default False.
    outdir: destination directory; created if missing.

    Returns a dict:
        "protein": the protein PDB path (split_chains=False), or a
            {chain_id: path} dict, one entry per chain (split_chains=True).
            Water is kept (crystallographic waters are routinely useful
            downstream); every other HETATM residue -- real ligands, ions,
            crystallization additives, glycosylation sugars, alike -- is
            stripped, since those are exactly what end up in "ligands" below
            instead.
        "ligands": list of {"path", "code", "chain", "resnum", "icode"}, one
            entry per non-water HETATM residue *instance* (see
            chem.ligand.list_ligand_instances) -- e.g. two copies of the same
            ligand code bound to different chains produce two entries, each
            its own SDF file. Each molecule's 3D coordinates come straight
            from `structure_path`, with bond orders/aromaticity restored
            against the PDB Chemical Component Dictionary (see
            chem.ligand.load_ligand). An instance load_ligand can't resolve a
            template for (e.g. a covalently-linked peptidomimetic ligand, or
            incomplete crystallographic density) is skipped with a warning
            (unless quiet) rather than aborting the whole split.
    """
    # Deferred: chem.ligand.extract imports from chem.protein.pocket at
    # module load time, so importing it at this module's top level would
    # deadlock whichever of the two subpackages is imported first.
    from ..ligand.extract import list_ligand_instances, load_ligand

    os.makedirs(outdir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(structure_path))[0]

    structure = _load_structure(structure_path)
    model = next(structure.get_models())
    exclude_residues = set(_hetero_residues(structure_path))

    io = PDBIO()
    io.set_structure(structure)
    if split_chains:
        protein_result = {}
        for chain in model:
            path = os.path.join(outdir, f"{stem}_protein_{chain.id}.pdb")
            io.save(path, _ProteinSelect(exclude_residues, chain_id=chain.id))
            protein_result[chain.id] = path
    else:
        protein_result = os.path.join(outdir, f"{stem}_protein.pdb")
        io.save(protein_result, _ProteinSelect(exclude_residues))

    ligand_results = []
    for instance in list_ligand_instances(structure_path, exclude=WATER):
        try:
            mol = load_ligand(
                structure_path,
                instance["code"],
                chain=instance["chain"],
                resnum=instance["resnum"],
                icode=instance["icode"],
            )
        except ValueError as e:
            if not is_quiet():
                print(f"skipping ligand instance {instance}: {e}", file=sys.stderr)
            continue

        sdf_path = os.path.join(
            outdir,
            f"{stem}_ligand_{instance['code']}_{instance['chain']}{instance['resnum']}{instance['icode']}.sdf",
        )
        writer = Chem.SDWriter(sdf_path)
        writer.write(mol)
        writer.close()
        ligand_results.append({"path": sdf_path, **instance})

    if not is_quiet():
        print(f"wrote protein PDB and {len(ligand_results)} ligand SDF(s) to {outdir}", file=sys.stderr)

    return {"protein": protein_result, "ligands": ligand_results}
