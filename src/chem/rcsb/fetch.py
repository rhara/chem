import os
import sys

import requests
from tqdm import tqdm

from ..ids import resolve_uniprot_accession_any
from ..verbosity import is_quiet, logged

RCSB_SEARCH_API = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_GRAPHQL_API = "https://data.rcsb.org/graphql"
RCSB_FILES_API = "https://files.rcsb.org/download"

FILETYPES = ("cif", "pdb", "both")

_SEARCH_PAGE_SIZE = 1000
_RESOLUTION_BATCH_SIZE = 200

_RESOLUTION_QUERY = """
query($ids: [String!]!) {
    entries(entry_ids: $ids) {
        rcsb_id
        rcsb_entry_info { resolution_combined }
    }
}
"""


def _search_entry_ids(accession):
    """Return every PDB entry id whose polymer entities are annotated with the given
    UniProt accession, via the RCSB Search API."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers"
                        ".reference_sequence_identifiers.database_accession",
                        "operator": "exact_match",
                        "value": accession,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "rcsb_polymer_entity_container_identifiers"
                        ".reference_sequence_identifiers.database_name",
                        "operator": "exact_match",
                        "value": "UniProt",
                    },
                },
            ],
        },
        "return_type": "entry",
        "request_options": {"paginate": {"start": 0, "rows": _SEARCH_PAGE_SIZE}},
    }

    entry_ids = []
    start = 0
    while True:
        query["request_options"]["paginate"]["start"] = start
        resp = requests.post(RCSB_SEARCH_API, json=query, timeout=60)
        if resp.status_code == 204:
            break  # no hits at all
        resp.raise_for_status()
        payload = resp.json()
        hits = payload.get("result_set", [])
        if not hits:
            break
        entry_ids.extend(hit["identifier"] for hit in hits)
        start += len(hits)
        if start >= payload.get("total_count", start):
            break
    return entry_ids


def _fetch_resolutions(entry_ids):
    """Return {entry_id: resolution_in_angstrom}, batched via the RCSB GraphQL API.
    Entries without a resolution (e.g. NMR structures) map to None.
    """
    resolutions = {}
    batches = [
        entry_ids[i : i + _RESOLUTION_BATCH_SIZE]
        for i in range(0, len(entry_ids), _RESOLUTION_BATCH_SIZE)
    ]
    for batch in tqdm(
        batches, desc="fetching resolution metadata", unit="batch", disable=is_quiet()
    ):
        resp = requests.post(
            RCSB_GRAPHQL_API,
            json={"query": _RESOLUTION_QUERY, "variables": {"ids": batch}},
            timeout=60,
        )
        resp.raise_for_status()
        for entry in resp.json()["data"]["entries"]:
            combined = (entry.get("rcsb_entry_info") or {}).get("resolution_combined")
            resolutions[entry["rcsb_id"]] = min(combined) if combined else None
    return resolutions


def _select_entries(resolutions, resolution_thres):
    """Filter an {entry_id: resolution_or_None} mapping down to the entry ids to
    download.

    resolution_thres=None keeps everything (including entries without a resolution,
    e.g. NMR structures). Otherwise only entries with a resolution <= resolution_thres
    are kept, and resolution-less entries are always dropped.
    """
    if resolution_thres is None:
        return list(resolutions)
    return [
        entry_id
        for entry_id, res in resolutions.items()
        if res is not None and res <= resolution_thres
    ]


def _download_one(entry_id, outdir, filetype):
    """Download the requested file(s) for one PDB entry into outdir, skipping any
    file that already exists on disk.

    Returns True if at least one file was written or already present (a format can
    be legitimately unavailable, e.g. no legacy .pdb file for some large cryo-EM
    structures).
    """
    extensions = ("cif", "pdb") if filetype == "both" else (filetype,)
    have_any = False
    for ext in extensions:
        path = os.path.join(outdir, f"{entry_id}.{ext}")
        if os.path.exists(path):
            have_any = True
            continue
        resp = requests.get(f"{RCSB_FILES_API}/{entry_id}.{ext}", timeout=60)
        if resp.status_code == 404:
            if not is_quiet():
                print(f"no .{ext} file available for {entry_id}, skipping", file=sys.stderr)
            continue
        resp.raise_for_status()
        with open(path, "wb") as f:
            f.write(resp.content)
        have_any = True
    return have_any


@logged
def download_structures(id, resolution_thres=None, outdir="data", filetype="cif"):
    """Download RCSB PDB structure files for a target into a directory.

    Resolves id to a UniProt accession, finds every PDB entry whose polymer entities
    are annotated with that accession (via the RCSB Search API), optionally filters
    by resolution, and downloads each qualifying entry's structure file(s). A file
    already present in outdir (same entry id and extension) is left as-is and not
    re-downloaded, so repeated calls only fetch what's missing.

    id: ChEMBL target id (e.g. "CHEMBL204"), UniProt accession (e.g. "P00734"),
        or UniProt entry name (e.g. "THRB_HUMAN").
    resolution_thres: optional maximum resolution in Angstrom, inclusive. When set,
        entries without a resolution (e.g. NMR structures) are excluded; when None
        (default), all entries are kept regardless of resolution.
    outdir: destination directory; created if missing.
    filetype: "cif" (default), "pdb", or "both" -- which structure file format(s) to
        download for each entry.

    Returns the number of entries for which at least one file is present (whether
    newly downloaded or already on disk from a previous call).
    """
    if filetype not in FILETYPES:
        raise ValueError(f"filetype must be one of {FILETYPES}")

    accession = resolve_uniprot_accession_any(id)
    entry_ids = _search_entry_ids(accession)
    if not entry_ids:
        raise ValueError(
            f"no PDB entries found for UniProt accession {accession} (from id '{id}')"
        )

    resolutions = _fetch_resolutions(entry_ids)
    selected = _select_entries(resolutions, resolution_thres)

    os.makedirs(outdir, exist_ok=True)
    n_written = 0
    for entry_id in tqdm(
        selected, desc=f"downloading structures to {outdir}", unit="entry", disable=is_quiet()
    ):
        if _download_one(entry_id, outdir, filetype):
            n_written += 1

    if not is_quiet():
        print(f"wrote {n_written} structures to {outdir}", file=sys.stderr)
    return n_written
