import os
import sys

import requests
from tqdm import tqdm

from ..ids import resolve_uniprot_accession_any
from ..verbosity import is_quiet, logged

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"

FILETYPES = ("cif", "pdb", "both")


def _fetch_predictions(accession):
    """Return the list of AlphaFold DB prediction entries for a UniProt accession
    (usually one, but very large proteins may be split into multiple fragments, and
    some targets without an official prediction have community-submitted
    alternatives instead). Each entry carries its own cifUrl/pdbUrl.
    """
    resp = requests.get(f"{ALPHAFOLD_API}/{accession}", timeout=30)
    if resp.status_code == 404:
        return []
    resp.raise_for_status()
    return resp.json()


def _select_entries(entries, plddt_thres):
    """Filter AlphaFold prediction entries by an optional minimum average pLDDT
    (globalMetricValue), inclusive. plddt_thres=None keeps everything.
    """
    if plddt_thres is None:
        return entries
    return [e for e in entries if e.get("globalMetricValue", 0) >= plddt_thres]


def _download_one(entry, outdir, filetype):
    """Download the requested file(s) for one AlphaFold prediction entry into outdir,
    skipping any file that already exists on disk.

    Returns True if at least one file is present (whether just downloaded or already
    on disk).
    """
    entry_id = entry["entryId"]
    urls = {"cif": entry.get("cifUrl"), "pdb": entry.get("pdbUrl")}
    extensions = ("cif", "pdb") if filetype == "both" else (filetype,)
    have_any = False
    for ext in extensions:
        path = os.path.join(outdir, f"{entry_id}.{ext}")
        if os.path.exists(path):
            have_any = True
            continue
        url = urls.get(ext)
        if not url:
            if not is_quiet():
                print(f"no .{ext} file available for {entry_id}, skipping", file=sys.stderr)
            continue
        resp = requests.get(url, timeout=60)
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
def download_structures(id, plddt_thres=None, outdir="data", filetype="cif"):
    """Download AlphaFold DB predicted structure files for a target into a directory.

    Resolves id to a UniProt accession and downloads every AlphaFold DB prediction
    entry for it. A file already present in outdir (same entry id and extension) is
    left as-is and not re-downloaded, so repeated calls only fetch what's missing.

    id: ChEMBL target id (e.g. "CHEMBL204"), UniProt accession (e.g. "P00734"),
        or UniProt entry name (e.g. "THRB_HUMAN").
    plddt_thres: optional minimum average pLDDT confidence (0-100), inclusive. When
        None (default), every prediction entry is kept regardless of confidence.
    outdir: destination directory; created if missing.
    filetype: "cif" (default), "pdb", or "both" -- which structure file format(s) to
        download for each entry.

    Returns the number of entries for which at least one file is present (whether
    newly downloaded or already on disk from a previous call).
    """
    if filetype not in FILETYPES:
        raise ValueError(f"filetype must be one of {FILETYPES}")

    accession = resolve_uniprot_accession_any(id)
    entries = _fetch_predictions(accession)
    if not entries:
        raise ValueError(
            f"no AlphaFold DB predictions found for UniProt accession {accession} (from id '{id}')"
        )

    selected = _select_entries(entries, plddt_thres)

    os.makedirs(outdir, exist_ok=True)
    n_written = 0
    for entry in tqdm(
        selected, desc=f"downloading structures to {outdir}", unit="entry", disable=is_quiet()
    ):
        if _download_one(entry, outdir, filetype):
            n_written += 1

    if not is_quiet():
        print(f"wrote {n_written} structures to {outdir}", file=sys.stderr)
    return n_written
