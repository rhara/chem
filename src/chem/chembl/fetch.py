import csv
import os
import statistics
import sys
from collections import defaultdict

import requests
from chembl_structure_pipeline import standardizer as csp
from rdkit import Chem, RDLogger
from rdkit.Chem.Descriptors import MolWt
from tqdm import tqdm

from ..ids import CHEMBL_API, resolve_target_chembl_id
from ..verbosity import is_quiet, logged

# chembl_structure_pipeline's standardize_mol/get_parent_mol run RDKit's C++
# Normalizer/Uncharger, which log "Running Normalizer" / "Running Uncharger" etc.
# to rdApp.info; that's internal noise, not something callers need to see.
RDLogger.DisableLog("rdApp.info")
RDLogger.DisableLog("rdApp.debug")


def _warm_up_standardizer():
    # The Normalizer's one-time "Initializing Normalizer" message is printed straight
    # to the process's stderr fd, bypassing RDLogger entirely, so it can only be
    # caught by redirecting fd 2 for this throwaway first call.
    devnull_fd = os.open(os.devnull, os.O_WRONLY)
    saved_stderr_fd = os.dup(2)
    try:
        os.dup2(devnull_fd, 2)
        # Must actually need neutralizing, or the Normalizer's lazy singleton never
        # gets constructed and the one-time init message is deferred to later calls.
        mol = Chem.MolFromSmiles("[NH3+]CCC(=O)[O-]")
        mol = csp.standardize_mol(mol)
        csp.get_parent_mol(mol)
    finally:
        os.dup2(saved_stderr_fd, 2)
        os.close(saved_stderr_fd)
        os.close(devnull_fd)


_warm_up_standardizer()

# One row per activity record (normalize_smiles=False).
ACTIVITY_FIELDS = [
    "molecule_chembl_id",
    "assay_chembl_id",
    "target_chembl_id",
    "document_chembl_id",
    "pchembl_value",
    "smiles",
    "mw",
]

# One row per unique normalized compound, duplicates aggregated (normalize_smiles=True).
AGGREGATED_FIELDS = [
    "parent_chembl_id",
    "target_chembl_id",
    "smiles",
    "mw",
    "n",
    "pchembl_mean",
    "pchembl_median",
    "pchembl_std",
]


def _normalize_smiles(smiles):
    """Standardize and desalt via the ChEMBL Structure Pipeline, returning the parent
    compound's canonical smiles (https://github.com/chembl/ChEMBL_Structure_Pipeline).
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = csp.standardize_mol(mol)
    mol, _excluded = csp.get_parent_mol(mol)
    return Chem.MolToSmiles(mol)


def _fetch_activity_pages(target_chembl_id):
    url = f"{CHEMBL_API}/activity.json"
    params = {
        "target_chembl_id": target_chembl_id,
        "pchembl_value__isnull": "false",
        "format": "json",
        "limit": 1000,
    }
    pbar = None
    try:
        while url:
            resp = requests.get(url, params=params, timeout=60)
            resp.raise_for_status()
            payload = resp.json()
            page = payload["activities"]
            meta = payload["page_meta"]
            if pbar is None:
                pbar = tqdm(
                    total=meta.get("total_count"),
                    desc=f"fetching activities for {target_chembl_id}",
                    unit="rec",
                    disable=is_quiet(),
                )
            pbar.update(len(page))
            yield from page
            next_path = meta.get("next")
            url = f"https://www.ebi.ac.uk{next_path}" if next_path else None
            params = None  # the next URL already carries the full query string
    finally:
        if pbar is not None:
            pbar.close()


def _aggregate_by_smiles(records):
    """Group records (dicts with molecule_chembl_id/target_chembl_id/smiles/mw/pchembl_value)
    by smiles and reduce each group to one row with pChEMBL count/mean/median/std.
    """
    groups = defaultdict(list)
    for r in records:
        groups[r["smiles"]].append(r)

    rows = []
    for smiles, items in groups.items():
        values = [it["pchembl_value"] for it in items]
        rows.append(
            [
                items[0]["molecule_chembl_id"],
                items[0]["target_chembl_id"],
                smiles,
                items[0]["mw"],
                len(values),
                round(statistics.fmean(values), 3),
                round(statistics.median(values), 3),
                round(statistics.pstdev(values), 3),
            ]
        )
    return rows


@logged
def download_activities(id, mw=None, normalize_smiles=False, output="activities.tsv"):
    """Download ChEMBL bioactivity records for a target into a tsv/csv file.

    Only records with a pChEMBL value are kept (activities without one are dropped).

    If normalize_smiles is True, compounds are desalted (see below) and then grouped by
    the resulting smiles: duplicate compounds are collapsed into a single row reporting
    the pChEMBL count, mean, median, and (population) standard deviation, plus a
    parent_chembl_id column holding one representative molecule_chembl_id for the group.
    If normalize_smiles is False, one row is written per activity record instead.

    id: ChEMBL target id (e.g. "CHEMBL204"), UniProt accession (e.g. "P00734"),
        or UniProt entry name (e.g. "THRB_HUMAN").
    mw: optional [lower, upper] molecular weight range, inclusive on both ends.
    normalize_smiles: if True, standardize and desalt each compound via the ChEMBL
        Structure Pipeline (get its parent compound) before computing mw, and aggregate
        duplicate compounds (see above).
    output: destination file path; delimiter is chosen from the extension
        (".tsv" -> tab, ".csv" -> comma, default tab).

    Returns the number of rows written.
    """
    if mw is not None:
        if len(mw) != 2:
            raise ValueError("mw must be [lower, upper]")
        lo, hi = mw

    target_chembl_id = resolve_target_chembl_id(id)
    delimiter = "," if str(output).lower().endswith(".csv") else "\t"

    records = []
    for act in _fetch_activity_pages(target_chembl_id):
        pchembl_value = act.get("pchembl_value")
        if pchembl_value is None:
            continue
        smiles = act.get("canonical_smiles")
        if not smiles:
            continue
        if normalize_smiles:
            smiles = _normalize_smiles(smiles)
            if smiles is None:
                continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            continue
        weight = MolWt(mol)
        if mw is not None and not (lo <= weight <= hi):
            continue
        records.append(
            {
                "molecule_chembl_id": act.get("molecule_chembl_id"),
                "assay_chembl_id": act.get("assay_chembl_id"),
                "target_chembl_id": act.get("target_chembl_id"),
                "document_chembl_id": act.get("document_chembl_id"),
                "pchembl_value": float(pchembl_value),
                "smiles": smiles,
                "mw": round(weight, 2),
            }
        )

    with open(output, "w", newline="") as f:
        writer = csv.writer(f, delimiter=delimiter)
        if normalize_smiles:
            writer.writerow(AGGREGATED_FIELDS)
            rows = _aggregate_by_smiles(records)
        else:
            writer.writerow(ACTIVITY_FIELDS)
            rows = [
                [
                    r["molecule_chembl_id"],
                    r["assay_chembl_id"],
                    r["target_chembl_id"],
                    r["document_chembl_id"],
                    r["pchembl_value"],
                    r["smiles"],
                    r["mw"],
                ]
                for r in records
            ]
        writer.writerows(rows)

    n_written = len(rows)
    if not is_quiet():
        print(f"wrote {n_written} rows to {output}", file=sys.stderr)
    return n_written
