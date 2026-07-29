import py3Dmol

from ..protein import SOLVENT_AND_IONS


def _ligand_resnames(pdb_text, exclude):
    """HET codes of every HETATM group in pdb_text, minus exclude."""
    return sorted(
        {line[17:20].strip() for line in pdb_text.splitlines() if line.startswith("HETATM")}
        - set(exclude)
    )


def render_protein(path, exclude=SOLVENT_AND_IONS, width=600, height=500):
    """Render a PDB structure file as an interactive py3Dmol view.

    Draws a rainbow (N -> C) cartoon backbone, plus any HETATM ligand group
    not in `exclude` as colored sticks (cartoon alone only draws the polymer
    backbone).

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

    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_text, "pdb")
    # "spectrum" alone defaults to 3Dmol.js's sinebow gradient, whose N-terminal
    # end drifts into purple/magenta; colorscheme="roygb" keeps it to blue (N) ->
    # cyan -> green -> yellow -> orange -> red (C), matching the usual convention.
    view.setStyle({"cartoon": {"color": "spectrum", "colorscheme": "roygb"}})
    if ligand_resnames:
        view.addStyle({"resn": ligand_resnames}, {"stick": {"colorscheme": "yellowCarbon"}})
    view.zoomTo()
    return view
