import sys
import time

import requests
from tqdm import tqdm

from ..verbosity import is_quiet

EBI_BLAST_API = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"

_TERMINAL_STATES = ("FINISHED", "FAILURE", "ERROR", "NOT_FOUND")


def _submit(sequence, database, email, matrix, expect, max_hits, title):
    resp = requests.post(
        f"{EBI_BLAST_API}/run",
        data={
            "email": email,
            "program": "blastp",
            "stype": "protein",
            "sequence": sequence,
            "database": database,
            "matrix": matrix,
            "exp": str(expect),
            "alignments": str(max_hits),
            "scores": str(max_hits),
            "title": title,
        },
        timeout=30,
    )
    resp.raise_for_status()
    return resp.text.strip()


def _status(job_id):
    resp = requests.get(f"{EBI_BLAST_API}/status/{job_id}", timeout=30)
    resp.raise_for_status()
    return resp.text.strip()


def _wait(job_id, poll_interval, timeout):
    """Poll job_id until it reaches a terminal state. There's no way to know a
    BLAST job's runtime up front, so progress (elapsed polls + latest status) is
    reported via tqdm rather than blocking silently.
    """
    start = time.time()
    with tqdm(desc=f"waiting for {job_id}", unit="poll", disable=is_quiet()) as pbar:
        while True:
            state = _status(job_id)
            pbar.set_postfix_str(state)
            pbar.update(1)
            if state in _TERMINAL_STATES:
                return state
            if time.time() - start > timeout:
                raise TimeoutError(f"{job_id} did not reach a terminal state within {timeout}s (last: {state})")
            time.sleep(poll_interval)


def _parse_hits(result_json):
    hits = []
    for h in result_json.get("hits", []):
        hsp = h["hit_hsps"][0]
        hits.append(
            {
                "accession": h["hit_acc"],
                "description": h["hit_def"],
                "identity_pct": round(hsp["hsp_identity"], 1),
                "align_len": hsp["hsp_align_len"],
                "evalue": hsp["hsp_expect"],
            }
        )
    return hits


def blastp(
    sequence,
    database,
    email,
    matrix="BLOSUM62",
    expect=1e-10,
    max_hits=50,
    poll_interval=10,
    timeout=600,
    title="chem.blast",
):
    """Run a protein BLAST (blastp) search via EBI's Job Dispatcher REST API
    (https://www.ebi.ac.uk/jdispatcher/docs/webservices/) and return the ranked
    hit list.

    Submits the search, polls until it finishes, then fetches and parses the
    JSON result.

    sequence: a protein sequence, plain or FASTA-formatted (EBI accepts either).
    database: which EBI-hosted database to search, e.g. "pdb" (every PDB chain
        -- hits are chain-level, e.g. "1ABC_A", so a homolog with an
        experimentally solved structure appears once per chain) or
        "uniprotkb_swissprot" (reviewed UniProt entries, one hit per protein).
        See the webservices docs above for the full list EBI hosts.
    email: contact email required by EBI's Job Dispatcher API (their
        abuse-prevention/contact policy -- not stored or used for anything
        else, and deliberately not passed to `chem.verbosity.logged`-style
        call logging, see below).
    matrix: substitution matrix, default "BLOSUM62".
    expect: E-value threshold (upper bound, inclusive), default 1e-10.
    max_hits: maximum number of hits to request, default 50 (also EBI's own
        per-search cap).
    poll_interval: seconds between status polls while the job runs, default 10.
    timeout: seconds to wait for the job to reach a terminal state before
        raising TimeoutError, default 600.
    title: job title EBI records for the search, default "chem.blast".

    Returns a list of dicts, one per hit, in EBI's own ranking order (best
    first):
        accession: the hit's database accession (e.g. a PDB chain id "1ABC_A"
            for database="pdb", or a UniProt accession for
            database="uniprotkb_swissprot").
        description: the hit's full description line.
        identity_pct: percent sequence identity over the aligned region.
        align_len: aligned region length (residues).
        evalue: the alignment's E-value.

    Raises RuntimeError if the search job itself fails (state FAILURE/ERROR/
    NOT_FOUND), and TimeoutError if it doesn't reach a terminal state within
    `timeout` seconds.

    Unlike every other public `chem` function, this one is *not* decorated
    with `chem.verbosity.logged`: that decorator logs every bound argument
    (including `email`) to stderr on each call, which would echo the caller's
    contact email into notebook output/logs for no benefit. `CHEM_QUIETNESS`
    still silences the tqdm progress bar and the one-line hit-count summary
    printed on success.
    """
    job_id = _submit(sequence, database, email, matrix, expect, max_hits, title)
    state = _wait(job_id, poll_interval, timeout)
    if state != "FINISHED":
        raise RuntimeError(f"BLAST search {job_id} against {database!r} ended in state {state}")

    resp = requests.get(f"{EBI_BLAST_API}/result/{job_id}/json", timeout=30)
    resp.raise_for_status()
    hits = _parse_hits(resp.json())

    if not is_quiet():
        print(f"{job_id}: {len(hits)} hit(s) against {database!r}", file=sys.stderr)

    return hits
