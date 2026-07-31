import pytest

import chem.blast.search as bs


def test_parse_hits_extracts_top_hsp_fields():
    result_json = {
        "hits": [
            {
                "hit_acc": "1ABC_A",
                "hit_def": "PDB:1ABC_A mol:protein length:300  Some Kinase",
                "hit_hsps": [{"hsp_identity": 44.756, "hsp_align_len": 288, "hsp_expect": 7e-82}],
            }
        ]
    }
    hits = bs._parse_hits(result_json)
    assert hits == [
        {
            "accession": "1ABC_A",
            "description": "PDB:1ABC_A mol:protein length:300  Some Kinase",
            "identity_pct": 44.8,
            "align_len": 288,
            "evalue": 7e-82,
        }
    ]


def test_parse_hits_empty_when_no_hits():
    assert bs._parse_hits({"hits": []}) == []


def test_submit_posts_expected_params(monkeypatch):
    captured = {}

    class FakeResponse:
        text = "  ncbiblast-job-123  "

        def raise_for_status(self):
            pass

    def fake_post(url, data, timeout):
        captured["url"] = url
        captured["data"] = data
        captured["timeout"] = timeout
        return FakeResponse()

    monkeypatch.setattr(bs.requests, "post", fake_post)

    job_id = bs._submit("MSEQ", "pdb", "me@example.com", "BLOSUM62", 1e-10, 50, "my-title")

    assert job_id == "ncbiblast-job-123"
    assert captured["url"] == f"{bs.EBI_BLAST_API}/run"
    assert captured["data"] == {
        "email": "me@example.com",
        "program": "blastp",
        "stype": "protein",
        "sequence": "MSEQ",
        "database": "pdb",
        "matrix": "BLOSUM62",
        "exp": "1e-10",
        "alignments": "50",
        "scores": "50",
        "title": "my-title",
    }


def test_status_strips_response_text(monkeypatch):
    class FakeResponse:
        text = "  RUNNING  "

        def raise_for_status(self):
            pass

    monkeypatch.setattr(bs.requests, "get", lambda url, timeout: FakeResponse())
    assert bs._status("job-1") == "RUNNING"


def test_wait_polls_until_terminal_state(monkeypatch):
    states = iter(["RUNNING", "RUNNING", "FINISHED"])
    monkeypatch.setattr(bs, "_status", lambda job_id: next(states))
    monkeypatch.setattr(bs.time, "sleep", lambda seconds: None)

    result = bs._wait("job-1", poll_interval=10, timeout=600)
    assert result == "FINISHED"


def test_wait_raises_timeout_error(monkeypatch):
    monkeypatch.setattr(bs, "_status", lambda job_id: "RUNNING")
    monkeypatch.setattr(bs.time, "sleep", lambda seconds: None)

    # Fake a clock that jumps straight past the timeout on the second read.
    clock = iter([0, 1000])
    monkeypatch.setattr(bs.time, "time", lambda: next(clock))

    with pytest.raises(TimeoutError):
        bs._wait("job-1", poll_interval=10, timeout=600)


def test_blastp_returns_parsed_hits_on_success(monkeypatch):
    monkeypatch.setattr(bs, "_submit", lambda *args: "job-1")
    monkeypatch.setattr(bs, "_wait", lambda job_id, poll_interval, timeout: "FINISHED")

    class FakeResponse:
        def raise_for_status(self):
            pass

        def json(self):
            return {
                "hits": [
                    {
                        "hit_acc": "P12345",
                        "hit_def": "SP:P12345 SOME_HUMAN",
                        "hit_hsps": [{"hsp_identity": 99.0, "hsp_align_len": 100, "hsp_expect": 0.0}],
                    }
                ]
            }

    monkeypatch.setattr(bs.requests, "get", lambda url, timeout: FakeResponse())

    hits = bs.blastp("MSEQ", "uniprotkb_swissprot", "me@example.com")
    assert hits == [
        {"accession": "P12345", "description": "SP:P12345 SOME_HUMAN", "identity_pct": 99.0, "align_len": 100, "evalue": 0.0}
    ]


def test_blastp_raises_runtime_error_on_job_failure(monkeypatch):
    monkeypatch.setattr(bs, "_submit", lambda *args: "job-1")
    monkeypatch.setattr(bs, "_wait", lambda job_id, poll_interval, timeout: "FAILURE")

    with pytest.raises(RuntimeError):
        bs.blastp("MSEQ", "pdb", "me@example.com")
