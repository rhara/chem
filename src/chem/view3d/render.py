import os
import re

import py3Dmol
from IPython.display import HTML, display

from ..protein import SOLVENT_AND_IONS

_RESOLUTION_RE = re.compile(r"^REMARK\s+2\s+RESOLUTION\.\s+([\d.]+)\s+ANGSTROMS\.", re.MULTILINE)


def _ligand_resnames(pdb_text, exclude):
    """HET codes of every HETATM group in pdb_text, minus exclude."""
    return sorted(
        {line[17:20].strip() for line in pdb_text.splitlines() if line.startswith("HETATM")}
        - set(exclude)
    )


def _chain_ids(pdb_text):
    """Distinct chain ids (PDB column 22) across every ATOM record, sorted."""
    chains = {
        line[21]
        for line in pdb_text.splitlines()
        if line.startswith("ATOM") and len(line) > 21
    }
    return sorted(c for c in chains if c.strip())


def _resolution(pdb_text):
    """Experimental resolution ("N.NN Å") from a legacy-PDB REMARK 2 record, or
    "N/A" if absent -- NMR structures, AlphaFold predictions, and files written
    by chem.protein.align (which doesn't preserve header/REMARK records) have
    no such record.
    """
    match = _RESOLUTION_RE.search(pdb_text)
    return f"{match.group(1)} Å" if match else "N/A"


def _caption(path, pdb_text, ligand_resnames):
    pdb_id = os.path.splitext(os.path.basename(path))[0]
    chains = ", ".join(_chain_ids(pdb_text)) or "N/A"
    ligands = ", ".join(ligand_resnames) or "none"
    return (
        f"PDB ID: {pdb_id} &nbsp;|&nbsp; Chain: {chains} &nbsp;|&nbsp; "
        f"Ligand: {ligands} &nbsp;|&nbsp; Resolution: {_resolution(pdb_text)}"
    )


def render_protein(path, exclude=SOLVENT_AND_IONS, width=600, height=500):
    """Render a PDB structure file as an interactive py3Dmol view.

    Displays a caption (PDB id, chain ids, ligand HET codes, and experimental
    resolution if present) above a rainbow (N -> C) cartoon backbone, plus any
    HETATM ligand group not in `exclude` as pink sticks (cartoon alone only
    draws the polymer backbone).

    path: path to a PDB file.
    exclude: HET codes to leave off the ligand sticks. Defaults to
        chem.protein.SOLVENT_AND_IONS (water/ions/crystallization additives);
        pass a superset (e.g. `SOLVENT_AND_IONS | {"NAG", "TYS"}`) to also
        exclude structure-specific non-ligand HETATM groups such as
        glycosylation sugars or modified residues.
    width / height: viewer size in pixels.

    Returns the py3Dmol view. Use this call as a notebook cell's last
    expression to display it, or call `.show()` on the result explicitly.
    """
    with open(path) as f:
        pdb_text = f.read()

    ligand_resnames = _ligand_resnames(pdb_text, exclude)

    display(HTML(f"<b>{_caption(path, pdb_text, ligand_resnames)}</b>"))

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, "pdb")
    # "spectrum" alone defaults to 3Dmol.js's sinebow gradient, whose N-terminal
    # end drifts into purple/magenta; colorscheme="roygb" keeps it to blue (N) ->
    # cyan -> green -> yellow -> orange -> red (C), matching the usual convention.
    view.setStyle({"cartoon": {"color": "spectrum", "colorscheme": "roygb"}})
    if ligand_resnames:
        # Solid pink, not a "*Carbon" colorscheme preset (3Dmol.js has no pinkCarbon,
        # and yellow -- one of the 8 it does have -- collides with the cartoon's
        # roygb spectrum).
        view.addStyle({"resn": ligand_resnames}, {"stick": {"color": "pink"}})
    view.zoomTo()
    return view
