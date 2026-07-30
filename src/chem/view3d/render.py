import os
import re
import uuid
from functools import lru_cache

import py3Dmol
import requests
from IPython.display import HTML, display
from rdkit import Chem
from rdkit.Chem import AllChem

from ..protein import SOLVENT_AND_IONS

_RESOLUTION_RE = re.compile(r"^REMARK\s+2\s+RESOLUTION\.\s+([\d.]+)\s+ANGSTROMS\.", re.MULTILINE)

RCSB_CHEMCOMP_API = "https://data.rcsb.org/rest/v1/core/chemcomp"

# py3Dmol pins 3Dmol.js 2.5.4 by default, which has no special rendering for
# MDL bond type 4 ("aromatic") -- it silently draws it as a plain single
# bond. 2.5.5 added `aromaticStyle` (single bond + inscribed ring torus, or
# a dashed second line), so aromatic rings need this newer build to actually
# look aromatic.
_3DMOL_JS_URL = "https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js"

COLORINGS = ("spectrum", "bfactor")
STYLES = ("cartoon", "surface")


def _ligand_resnames(pdb_text, exclude):
    """HET codes of every HETATM group in pdb_text, minus exclude."""
    return sorted(
        {line[17:20].strip() for line in pdb_text.splitlines() if line.startswith("HETATM")}
        - set(exclude)
    )


def _ligand_instances(pdb_text, ligand_resnames):
    """(resname, chain, resnum, icode) for every distinct HETATM group in
    pdb_text whose resname is in ligand_resnames, in first-appearance order."""
    ligand_resnames = set(ligand_resnames)
    instances = []
    for line in pdb_text.splitlines():
        if not line.startswith("HETATM"):
            continue
        resname = line[17:20].strip()
        if resname not in ligand_resnames:
            continue
        key = (resname, line[21], line[22:26].strip(), line[26].strip())
        if key not in instances:
            instances.append(key)
    return instances


def _instance_block(pdb_text, resname, chain, resnum, icode):
    """The HETATM lines of one ligand instance, as a standalone PDB block."""
    lines = [
        line
        for line in pdb_text.splitlines()
        if line.startswith("HETATM")
        and line[17:20].strip() == resname
        and line[21] == chain
        and line[22:26].strip() == resnum
        and line[26].strip() == icode
    ]
    return "\n".join(lines) + "\nEND\n"


@lru_cache(maxsize=None)
def _template_mol(resname):
    """A Hs-stripped RDKit mol built from the canonical SMILES the RCSB
    Chemical Component Dictionary has on file for a 3-letter ligand HET code,
    or None if the lookup fails (unknown/custom code, or no network access).
    Cached since the same code often recurs across many ligand instances.
    """
    try:
        resp = requests.get(f"{RCSB_CHEMCOMP_API}/{resname}", timeout=15)
        resp.raise_for_status()
        descriptors = resp.json()["pdbx_chem_comp_descriptor"]
        smiles = next(
            d["descriptor"] for d in descriptors if d["type"] == "SMILES_CANONICAL"
        )
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[#1]"))
    except Exception:
        return None


def _ligand_molblock(pdb_text, resname, chain, resnum, icode):
    """A MolBlock for one ligand instance with bond orders assigned by
    matching its distance-perceived connectivity against the RCSB Chemical
    Component Dictionary's reference structure for its HET code, or None if
    that lookup or the match fails -- ligand atom positions from the PDB file
    are otherwise bond-order-blind (no CONECT bond-order records), so
    distance alone can't tell a ring's double bonds from its single ones.

    Aromatic rings are written unkekulized (MDL bond type 4, "aromatic")
    rather than as one arbitrarily-chosen single/double resonance form --
    3Dmol.js (2.5.5+) renders bond type 4 as a single bond plus an inscribed
    ring torus (`view3d._build_view`'s `aromaticStyle: "circle"`), the usual
    textbook depiction of a delocalized ring, instead of a kekulized
    depiction of one particular resonance structure.
    """
    template = _template_mol(resname)
    if template is None:
        return None
    block = _instance_block(pdb_text, resname, chain, resnum, icode)
    mol = Chem.MolFromPDBBlock(block, sanitize=False, removeHs=False)
    if mol is not None:
        mol = Chem.RemoveHs(mol, sanitize=False)
    if mol is None or mol.GetNumAtoms() != template.GetNumAtoms():
        return None
    try:
        assigned = AllChem.AssignBondOrdersFromTemplate(template, mol)
    except Exception:
        return None
    return Chem.MolToMolBlock(assigned, kekulize=False)


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


def _caption_lines(path, pdb_text, ligand_resnames):
    """[PDB id, chain ids, ligand HET codes, resolution] as plain-text lines,
    one field per line."""
    pdb_id = os.path.splitext(os.path.basename(path))[0]
    chains = ", ".join(_chain_ids(pdb_text)) or "N/A"
    ligands = ", ".join(ligand_resnames) or "none"
    return [
        f"PDB ID: {pdb_id}",
        f"Chain: {chains}",
        f"Ligand: {ligands}",
        f"Resolution: {_resolution(pdb_text)}",
    ]


def _build_view(pdb_text, ligand_resnames, width, height, coloring, bfactor_range, style="cartoon"):
    view = py3Dmol.view(width=width, height=height, js=_3DMOL_JS_URL)
    view.addModel(pdb_text, "pdb")
    if coloring == "bfactor":
        bmin, bmax = bfactor_range
        # e.g. AlphaFold stores per-residue pLDDT confidence in the B-factor column.
        color_style = {"colorscheme": {"prop": "b", "gradient": "roygb", "min": bmin, "max": bmax}}
    else:
        # "spectrum" alone defaults to 3Dmol.js's sinebow gradient, whose N-terminal
        # end drifts into purple/magenta; colorscheme="roygb" keeps it to blue (N) ->
        # cyan -> green -> yellow -> orange -> red (C), matching the usual convention.
        color_style = {"color": "spectrum", "colorscheme": "roygb"}
    if style == "surface":
        # A VDW (van der Waals) isosurface -- effectively the union of a Gaussian-ish
        # blob at every atom's position, computed by 3Dmol.js via marching cubes over
        # a grid built from atomic radii -- rather than individual per-atom spheres.
        # Restricted to non-HETATM atoms so a bound ligand's own stick rendering
        # (below) isn't swallowed by an opaque protein blob; opacity < 1 keeps it
        # from fully hiding a ligand buried in a pocket either way.
        view.addSurface("VDW", {**color_style, "opacity": 0.85}, {"hetflag": False})
    else:
        view.setStyle({"cartoon": color_style})
    # "magentaCarbon" -- one of 3Dmol.js's 8 "*Carbon" presets -- colors carbon
    # magenta (bold and high-contrast against the cartoon's roygb spectrum,
    # unlike yellowCarbon, which blends into it) while leaving every other
    # element (O, N, S, halogens, ...) at its normal CPK color, so heteroatoms
    # read the same way they do in any other molecular viewer. "circle" draws
    # aromatic (bond order 4) rings as a plain bond plus an inscribed torus --
    # the usual textbook symbol for a delocalized ring -- instead of guessing
    # at a single/double kekulized resonance form.
    ligand_style = {"stick": {"colorscheme": "magentaCarbon", "aromaticStyle": "circle"}}
    for resname, chain, resnum, icode in _ligand_instances(pdb_text, ligand_resnames):
        molblock = _ligand_molblock(pdb_text, resname, chain, resnum, icode)
        if molblock is not None:
            # A separate small-molecule model with real bond orders, so
            # aromatic rings render as such instead of a plain single-bonded
            # stick skeleton.
            view.addModel(molblock, "mol")
            view.setStyle({"model": -1}, ligand_style)
        else:
            # Bond-order lookup failed (unknown HET code, no network, ...) --
            # fall back to the old distance-only stick rendering for just
            # this instance.
            resi = f"{resnum}{icode}" if icode else int(resnum)
            view.addStyle({"resn": resname, "chain": chain, "resi": resi}, ligand_style)
    view.zoomTo()
    return view


def render_protein(
    path,
    exclude=SOLVENT_AND_IONS,
    width=600,
    height=500,
    coloring="spectrum",
    bfactor_range=(50, 90),
    style="cartoon",
):
    """Display a PDB structure file as an interactive py3Dmol view, followed by
    a caption.

    Shows the protein backbone -- as a cartoon or a solid volume per `style`,
    colored per `coloring` -- plus any HETATM ligand group not in `exclude` as
    sticks (neither backbone style draws ligands on its own), with a caption
    to its right listing the PDB id, chain ids, ligand HET codes, and
    experimental resolution if present. Ligand carbons are colored magenta
    for contrast against the backbone; every other element (O, N, S,
    halogens, ...) keeps its normal CPK color. Each ligand's bond orders are
    looked up from the RCSB Chemical Component Dictionary by HET code (one
    web request per distinct code, the PDB file's own HETATM records carry no
    bond-order information) and matched onto its 3D coordinates; aromatic
    rings render unkekulized, as a plain bond plus an inscribed ring torus
    (the usual textbook symbol for a delocalized ring) rather than one
    arbitrarily-chosen single/double resonance form. A ligand whose code
    isn't found there, whose atom count doesn't match the Dictionary's
    reference structure (e.g. a partially resolved ligand missing some
    atoms), or with no network access, falls back to plain single-bonded
    sticks.

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
    style: "cartoon" (default) -- ribbon backbone; or "surface" -- a solid
        van der Waals volume (the union of a smooth blob at every backbone
        atom, computed by 3Dmol.js via marching cubes over a grid of atomic
        radii, rather than a ribbon or individual per-atom spheres),
        translucent (opacity 0.85) so a ligand bound underneath still shows
        through. Restricted to non-HETATM atoms either way, so the ligand
        keeps its own separate stick rendering rather than being enveloped
        by the volume.

    The view sits in a light-gray-bordered frame -- exactly the area where
    3Dmol.js's mouse controls (rotate/zoom/pan) take over, so the border
    orients the user to where dragging does something different -- with the
    caption beside it on the right, one field per line.

    Displays the view and caption directly as a side effect (no return
    value) -- just call it, no need to chain `.show()` or use it as a cell's
    last expression.
    """
    if coloring not in COLORINGS:
        raise ValueError(f"coloring must be one of {COLORINGS}")
    if style not in STYLES:
        raise ValueError(f"style must be one of {STYLES}")

    with open(path) as f:
        pdb_text = f.read()

    ligand_resnames = _ligand_resnames(pdb_text, exclude)
    view = _build_view(pdb_text, ligand_resnames, width, height, coloring, bfactor_range, style)

    # A same-sized placeholder div, bordered up front so the frame is visible
    # immediately (3Dmol.js loads its viewer asynchronously from a CDN); once
    # loaded, view.insert() moves the actual 3Dmol viewer div into it.
    frame_id = f"chem-view3d-{uuid.uuid4().hex}"
    caption_html = "".join(
        f"<div><b>{line}</b></div>" for line in _caption_lines(path, pdb_text, ligand_resnames)
    )
    display(
        HTML(
            '<div style="display:flex; align-items:flex-start; gap:16px;">'
            f'<div id="{frame_id}" style="display:inline-block; border:1px solid #ccc; '
            f'width:{width}px; height:{height}px;"></div>'
            f"<div>{caption_html}</div>"
            "</div>"
        )
    )
    view.insert(frame_id)
