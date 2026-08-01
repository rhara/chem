"""Shared identifier resolution across chem subpackages: interconverting ChEMBL
target ids, UniProt accessions, and UniProt entry names (mnemonics) via the ChEMBL
and UniProt REST APIs.
"""

import re

import requests

CHEMBL_API = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT_API = "https://rest.uniprot.org/uniprotkb"

CHEMBL_ID_RE = re.compile(r"^CHEMBL\d+$", re.IGNORECASE)
# Official UniProtKB accession shape; anything else is treated as an entry name (mnemonic).
UNIPROT_ACCESSION_RE = re.compile(
    r"^([A-NR-Z][0-9][A-Z0-9]{3}[0-9]|[OPQ][0-9][A-Z0-9]{3}[0-9])$", re.IGNORECASE
)


def resolve_uniprot_accession(id_):
    """Resolve a UniProt accession (e.g. "P00734") or entry name/mnemonic
    (e.g. "THRB_HUMAN") to a canonical UniProt accession."""
    # UniProt's "accession" query field rejects anything not shaped like an accession
    # (400 Bad Request), so entry names (e.g. "THRB_HUMAN") must go through "id" instead.
    field = "accession" if UNIPROT_ACCESSION_RE.match(id_) else "id"
    resp = requests.get(
        f"{UNIPROT_API}/search",
        params={
            "query": f"{field}:{id_}",
            "fields": "accession,id",
            "format": "json",
            "size": 10,
        },
        timeout=30,
    )
    resp.raise_for_status()
    results = resp.json().get("results", [])
    if not results:
        raise ValueError(f"could not resolve '{id_}' to a UniProt accession")
    if field == "id":
        # UniProt's "id" query is a relevance-ranked text search, not an exact-match
        # lookup -- e.g. "id:CDK1_HUMAN" ranks O14519 (CDKA1_HUMAN, an unrelated
        # protein) first, ahead of the actual exact match P06493 (CDK1_HUMAN) at
        # position 2. Prefer a case-insensitive exact match on uniProtkbId among the
        # fetched results over blindly trusting rank; only fall back to the top hit
        # if none of them is an exact match (e.g. a deprecated/renamed mnemonic).
        exact = next((r for r in results if r["uniProtkbId"].lower() == id_.lower()), None)
        if exact is not None:
            return exact["primaryAccession"]
    return results[0]["primaryAccession"]


def resolve_target_chembl_id(id_):
    """Resolve a ChEMBL target id, UniProt accession, or UniProt entry name to a
    ChEMBL target id (e.g. "CHEMBL204")."""
    if CHEMBL_ID_RE.match(id_):
        return id_.upper()

    accession = resolve_uniprot_accession(id_)
    resp = requests.get(
        f"{CHEMBL_API}/target.json",
        params={"target_components__accession": accession, "format": "json"},
        timeout=30,
    )
    resp.raise_for_status()
    targets = resp.json().get("targets", [])
    if not targets:
        raise ValueError(
            f"no ChEMBL target found for UniProt accession {accession} (from id '{id_}')"
        )
    # Prefer a single-protein target over complexes/families sharing the same component.
    targets.sort(key=lambda t: t.get("target_type") != "SINGLE PROTEIN")
    return targets[0]["target_chembl_id"]


def resolve_uniprot_accession_any(id_):
    """Resolve a ChEMBL target id, UniProt accession, or UniProt entry name to a
    UniProt accession (e.g. "P00734")."""
    if not CHEMBL_ID_RE.match(id_):
        return resolve_uniprot_accession(id_)

    resp = requests.get(f"{CHEMBL_API}/target/{id_.upper()}.json", timeout=30)
    resp.raise_for_status()
    for component in resp.json().get("target_components", []):
        if component.get("component_type") == "PROTEIN" and component.get("accession"):
            return component["accession"]
    raise ValueError(f"no UniProt accession found for ChEMBL target '{id_}'")
