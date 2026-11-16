"""fetch_bvbrc_annotations.py (both PICRUSt2-comparison copies), offline.

Re-score D10 and the Run 2 network rules: the PGFam fetch covers the references of BOTH selection
files by default (the depth-20 legacy arms were projected without their families until 2026-11),
requests carry prefetch_gene_families' User-Agent, and two requests are at least --min-interval
seconds apart.
"""
import importlib.util
import io
import json
import os
import re

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
COPIES = ["region_validation/picrust2_comparison/scripts/fetch_bvbrc_annotations.py",
          "region_validation_d100/picrust2_comparison/scripts/fetch_bvbrc_annotations.py"]


def _load(path):
    name = "fba_" + path.split("/")[0]
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROOT, path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_copies_are_identical():
    a, b = (open(os.path.join(ROOT, p), "rb").read() for p in COPIES)
    assert a == b


@pytest.mark.parametrize("path", COPIES)
def test_user_agent_is_prefetch_gene_families(path):
    mod = _load(path)
    src = open(os.path.join(ROOT, "prefetch_gene_families.py")).read()
    ua = re.search(r'^USER_AGENT = "([^"]+)"', src, re.M).group(1)
    assert mod.UA == ua


class _Resp(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False


@pytest.mark.parametrize("path", COPIES)
def test_requests_are_rate_limited_and_carry_the_user_agent(monkeypatch, path):
    mod = _load(path)
    clock = {"t": 100.0}
    sleeps, seen = [], []

    def fake_urlopen(req, timeout=None):
        seen.append((clock["t"], req.get_header("User-agent")))
        clock["t"] += 0.2                      # the request itself takes 0.2 s
        return _Resp(json.dumps([{"family_id": "PGF_1", "family_product": "x"}]).encode())

    def fake_sleep(s):
        sleeps.append(s)
        clock["t"] += s

    monkeypatch.setattr(mod.urllib.request, "urlopen", fake_urlopen)
    monkeypatch.setattr(mod.time, "sleep", fake_sleep)
    monkeypatch.setattr(mod.time, "monotonic", lambda: clock["t"])
    monkeypatch.setattr(mod, "MIN_INTERVAL", 1.5)
    mod._LAST[0] = 0.0
    for _ in range(3):
        mod.rql_post("protein_family_ref", "in(family_id,(PGF_1))")
    starts = [t for t, _ in seen]
    assert all(b - a >= 1.5 - 1e-9 for a, b in zip(starts, starts[1:]))
    assert sleeps == pytest.approx([1.3, 1.3])
    assert {ua for _, ua in seen} == {mod.UA}


def _write(path, obj):
    path.write_text(json.dumps(obj))
    return str(path)


@pytest.mark.parametrize("path", COPIES)
def test_selected_reps_unions_every_file(tmp_path, path):
    mod = _load(path)
    a = _write(tmp_path / "a.json", {"g1__V4": {"include": {"reps": ["r1"]},
                                                "exclude_0.987": {"reps": ["r2"]},
                                                "best_nonself_identity": 0.9},
                                     "_meta": {"reps": ["nope"]}})
    b = _write(tmp_path / "b.json", {"g1__V4": {"include": {"reps": ["r3"]},
                                                "exclude_0.987": {"reps": []}}})
    assert mod.selected_reps([a, b]) == {"r1", "r2", "r3"}


@pytest.mark.parametrize("path", COPIES)
def test_pgfams_default_selection_is_both_files(tmp_path, monkeypatch, path):
    mod = _load(path)
    truth = _write(tmp_path / "truth.json", {"g1__V4": {"src_genome_id": "g1"}})
    monkeypatch.setattr(mod.C.PATHS, "truth_json", truth)
    got = {}

    class Stop(Exception):
        pass

    def capture(paths):
        got["paths"] = paths
        raise Stop

    monkeypatch.setattr(mod, "selected_reps", capture)
    monkeypatch.setattr(mod.sys, "argv", ["fetch_bvbrc_annotations.py", "pgfams"])
    with pytest.raises(Stop):
        mod.main()
    names = [os.path.basename(p) for p in got["paths"]]
    assert names == ["prfba_selection.json", "prfba_selection_legacy.json"]
    monkeypatch.setattr(mod.sys, "argv", ["x", "pgfams", "--selection", "/a.json,/b.json"])
    with pytest.raises(Stop):
        mod.main()
    assert got["paths"] == ["/a.json", "/b.json"]
