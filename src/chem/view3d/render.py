import os
import re

import py3Dmol
from IPython.display import HTML, display

from ..protein import SOLVENT_AND_IONS

_RESOLUTION_RE = re.compile(r"^REMARK\s+2\s+RESOLUTION\.\s+([\d.]+)\s+ANGSTROMS\.", re.MULTILINE)

COLORINGS = ("spectrum", "bfactor")


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


def _build_view(pdb_text, ligand_resnames, width, height, coloring, bfactor_range):
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, "pdb")
    if coloring == "bfactor":
        bmin, bmax = bfactor_range
        # e.g. AlphaFold stores per-residue pLDDT confidence in the B-factor column.
        view.setStyle({"cartoon": {"colorscheme": {"prop": "b", "gradient": "roygb", "min": bmin, "max": bmax}}})
    else:
        # "spectrum" alone defaults to 3Dmol.js's sinebow gradient, whose N-terminal
        # end drifts into purple/magenta; colorscheme="roygb" keeps it to blue (N) ->
        # cyan -> green -> yellow -> orange -> red (C), matching the usual convention.
        view.setStyle({"cartoon": {"color": "spectrum", "colorscheme": "roygb"}})
    if ligand_resnames:
        # Solid magenta -- bold and high-contrast against the cartoon's roygb
        # spectrum (unlike yellow, one of 3Dmol.js's 8 "*Carbon" presets, which
        # blends into it).
        view.addStyle({"resn": ligand_resnames}, {"stick": {"color": "magenta"}})
    view.zoomTo()
    return view


def render_protein(
    path,
    exclude=SOLVENT_AND_IONS,
    width=600,
    height=500,
    coloring="spectrum",
    bfactor_range=(50, 90),
):
    """Display a PDB structure file as an interactive py3Dmol view, followed by
    a caption.

    Shows a cartoon backbone -- colored per `coloring` -- plus any HETATM
    ligand group not in `exclude` as magenta sticks (cartoon alone only draws
    the polymer backbone), then -- below the view -- a caption with the PDB
    id, chain ids, ligand HET codes, and experimental resolution if present.

    path: path to a PDB file.
    exclude: HET codes to leave off the ligand sticks. Defaults to
        chem.protein.SOLVENT_AND_IONS (water/ions/crystallization additives);
        pass a superset (e.g. `SOLVENT_AND_IONS | {"NAG", "TYS"}`) to also
        exclude structure-specific non-ligand HETATM groups such as
        glycosylation sugars or modified residues.
    width / height: viewer size in pixels.
    coloring: "spectrum" (default) -- rainbow N -> C by residue position; or
        "bfactor" -- rainbow by the file's per-atom B-factor column, e.g.
        AlphaFold's per-residue pLDDT confidence.
    bfactor_range: (min, max) the "bfactor" gradient is scaled over; ignored
        for "spectrum". Defaults to AlphaFold's pLDDT confidence convention
        (50-90); pass the structure's own B-factor range for crystallographic
        temperature factors.

    Displays the view and caption directly as a side effect (no return
    value) -- just call it, no need to chain `.show()` or use it as a cell's
    last expression.
    """
    if coloring not in COLORINGS:
        raise ValueError(f"coloring must be one of {COLORINGS}")

    with open(path) as f:
        pdb_text = f.read()

    ligand_resnames = _ligand_resnames(pdb_text, exclude)

    view = _build_view(pdb_text, ligand_resnames, width, height, coloring, bfactor_range)
    view.show()
    display(HTML(f"<b>{_caption(path, pdb_text, ligand_resnames)}</b>"))
