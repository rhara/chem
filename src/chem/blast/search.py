import sys
import time

import requests
from tqdm import tqdm

from ..verbosity import is_quiet

EBI_BLAST_API = "https://www.ebi.ac.uk/Tools/services/rest/ncbiblast"

_TERMINAL_STATES = ("FINISHED", "FAILURE", "ERROR", "NOT_FOUND")

# EBI's ncbiblast "exp" parameter isn't a free-form float: only these exact
# strings are accepted (fetched from /parameterdetails/exp). plain str(expect)
# breaks for several of them -- e.g. str(1e-3) == "0.001", which EBI rejects
# with a 400 -- so every accepted value is looked up by its numeric value here
# instead of reformatted on the fly.
_EXPECT_VALUES = {
    1e-200: "1e-200",
    1e-100: "1e-100",
    1e-50: "1e-50",
    1e-10: "1e-10",
    1e-5: "1e-5",
    1e-4: "1e-4",
    1e-3: "1e-3",
    1e-2: "1e-2",
    1e-1: "1e-1",
    1.0: "1.0",
    10: "10",
    100: "100",
    1000: "1000",
}

# EBI's "alignments"/"scores" parameters (both driven by max_hits here) are
# likewise a fixed enum, not an arbitrary integer.
_MAX_HITS_VALUES = (0, 5, 10, 20, 50, 100, 150, 200, 250, 500, 750, 1000)


def _format_expect(expect):
    try:
        return _EXPECT_VALUES[expect]
    except (KeyError, TypeError):
        raise ValueError(
            f"expect={expect!r} is not one of EBI's accepted E-value thresholds: "
            f"{sorted(_EXPECT_VALUES, reverse=True)}"
        ) from None


def _format_max_hits(max_hits):
    if max_hits not in _MAX_HITS_VALUES:
        raise ValueError(f"max_hits={max_hits!r} is not one of EBI's accepted values: {_MAX_HITS_VALUES}")
    return str(int(max_hits))


def _submit(sequence, database, email, matrix, expect, max_hits, title):
    hits_str = _format_max_hits(max_hits)
    resp = requests.post(
        f"{EBI_BLAST_API}/run",
        data={
            "email": email,
            "program": "blastp",
            "stype": "protein",
            "sequence": sequence,
            "database": database,
            "matrix": matrix,
            "exp": _format_expect(expect),
            "alignments": hits_str,
            "scores": hits_str,
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
    email="user@example.com",
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
        call logging, see below). Defaults to the placeholder
        "user@example.com"; pass your own if you're making many calls.
    matrix: substitution matrix, default "BLOSUM62".
    expect: E-value threshold (upper bound, inclusive), default 1e-10 --
        notably stricter than EBI's own tool default (10). Not a free-form
        float: EBI's API only accepts one of a fixed set of values (raises
        ValueError otherwise) -- 1e-200, 1e-100, 1e-50, 1e-10, 1e-5, 1e-4,
        1e-3, 1e-2, 1e-1, 1.0, 10, 100, 1000. Pick a larger one (e.g. 1e-3,
        1.0) to surface more distant/lower-identity homologs.
    max_hits: maximum number of hits to request, default 50. Also a fixed
        EBI enum (ValueError otherwise) -- 0, 5, 10, 20, 50, 100, 150, 200,
        250, 500, 750, 1000. EBI's own default is also 50, but a search
        often has more hits than that within the `expect` cutoff -- raise
        this to see further down the ranked list rather than just the
        closest matches.
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
