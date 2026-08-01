import requests

from .. import __version__
from ..ids import resolve_uniprot_accession_any

UNIPROT_API = "https://rest.uniprot.org/uniprotkb"

# Pharos target development levels: https://pharos.nih.gov/about
_PHAROS_LEVELS = {
    "Tclin": "target of an approved drug",
    "Tchem": "has potent small-molecule ligands",
    "Tbio": "biology characterized, no known drug/chemical probe",
    "Tdark": "understudied",
}


def _extract_properties(entry):
    protein_desc = entry["proteinDescription"]["recommendedName"]
    ec_number = protein_desc.get("ecNumbers", [{}])[0].get("value")

    gene = entry.get("genes", [{}])[0]
    gene_name = gene.get("geneName", {}).get("value")
    gene_synonyms = ", ".join(s["value"] for s in gene.get("synonyms", []))

    family = next(
        (t["value"] for c in entry["comments"] if c["commentType"] == "SIMILARITY" for t in c["texts"]),
        None,
    )
    function = next(
        (t["value"] for c in entry["comments"] if c["commentType"] == "FUNCTION" for t in c["texts"]),
        None,
    )
    subcellular_location = ", ".join(
        loc["location"]["value"]
        for c in entry["comments"]
        if c["commentType"] == "SUBCELLULAR LOCATION"
        for loc in c["subcellularLocations"]
    )

    # Matched loosely on "kinase" (not just the literal "Protein kinase" label UniProt
    # uses for CDK-family entries) so this stays meaningful for other kinase families;
    # non-kinase targets simply get None here.
    kinase_domain = next(
        (
            f["location"]
            for f in entry["features"]
            if f["type"] == "Domain" and "kinase" in f.get("description", "").lower()
        ),
        None,
    )
    kinase_domain_range = (
        f"{kinase_domain['start']['value']}-{kinase_domain['end']['value']}" if kinase_domain else None
    )
    active_site_residue = next(
        (f["location"]["start"]["value"] for f in entry["features"] if f["type"] == "Active site"), None
    )

    xrefs = {x["database"]: x for x in entry["uniProtKBCrossReferences"]}
    n_pdb_xrefs = sum(1 for x in entry["uniProtKBCrossReferences"] if x["database"] == "PDB")
    chembl_target_id = xrefs.get("ChEMBL", {}).get("id")

    pharos_raw = next(
        (p["value"] for p in xrefs.get("Pharos", {}).get("properties", []) if p["key"] == "DevelopmentLevel"),
        None,
    )
    pharos_development_level = f"{pharos_raw} ({_PHAROS_LEVELS[pharos_raw]})" if pharos_raw else None

    return {
        "entry_name": entry["uniProtkbId"],
        "accession": entry["primaryAccession"],
        "protein_name": protein_desc["fullName"]["value"],
        "gene_name": f"{gene_name} (synonyms: {gene_synonyms})" if gene_synonyms else gene_name,
        "organism": entry["organism"]["scientificName"],
        "sequence_length": entry["sequence"]["length"],
        "ec_number": ec_number,
        "family": family,
        "function": function,
        "subcellular_location": subcellular_location,
        "kinase_domain_range": kinase_domain_range,
        "active_site_residue": active_site_residue,
        "n_pdb_xrefs": n_pdb_xrefs,
        "has_alphafold_model": "AlphaFoldDB" in xrefs,
        "chembl_target_id": chembl_target_id,
        "has_bindingdb_entry": "BindingDB" in xrefs,
        "pharos_development_level": pharos_development_level,
        "protein_existence": entry["proteinExistence"],
        "annotation_score": entry["annotationScore"],
    }


def summary(id_):
    """Fetch a UniProt entry (by accession, e.g. "Q8IZL9"; entry name/mnemonic,
    e.g. "CDK20_HUMAN"; or ChEMBL target id, e.g. "CHEMBL3559690") and return a
    dict of properties useful for drug-discovery triage: identifiers (entry_name,
    accession), function/family/organism/localization, kinase domain range and
    catalytic active site residue (when annotated), and cross-references signaling
    data availability -- PDB structure count, whether an AlphaFold model or
    BindingDB entry exists, the ChEMBL target id, and the Pharos target
    development level (Tclin/Tchem/Tbio/Tdark, i.e. how druggable it already is).
    """
    accession = resolve_uniprot_accession_any(id_)
    resp = requests.get(f"{UNIPROT_API}/{accession}.json", timeout=30)
    resp.raise_for_status()
    return _extract_properties(resp.json())


def get_fasta(id_, email="user@example.com"):
    """Fetch a UniProt entry's sequence as a FASTA string. `id_` accepts a UniProt
    accession (e.g. "Q8IZL9"), entry name/mnemonic (e.g. "CDK20_HUMAN"), or ChEMBL
    target id (e.g. "CHEMBL3559690") -- same resolution as `summary`.

    `email` isn't required by UniProt's REST API, but is sent as a contact address
    in the request's `User-Agent` header per UniProt's own API usage guidelines
    (https://www.uniprot.org/help/api); pass your own if making many calls.
    """
    accession = resolve_uniprot_accession_any(id_)
    headers = {"User-Agent": f"chem/{__version__} ({email})"}
    resp = requests.get(f"{UNIPROT_API}/{accession}.fasta", headers=headers, timeout=30)
    resp.raise_for_status()
    return resp.text
