"""Offline tests for manuscript_numbers/derived_numbers.py (rescore plan item W6): the
pre/post input resolver and the pure recipes, on tiny fixtures. Run:
    ~/Documents/py_venv/bin/python -m pytest -q tests/test_derived_numbers.py
"""
import importlib.util
import os

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_spec = importlib.util.spec_from_file_location(
    "derived_numbers", os.path.join(ROOT, "manuscript_numbers", "derived_numbers.py"))
D = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(D)


def _write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        fh.write(text)


@pytest.fixture
def tree(tmp_path):
    repo, w = tmp_path / "repo", tmp_path / "w"
    # a.csv: moved aside to pre_repair, tracked at the tag, and live (rewritten)
    _write(str(w / "pre_repair" / "d" / "a.csv"), "pre_repair")
    _write(str(w / "start" / "d" / "a.csv"), "start")
    _write(str(repo / "d" / "a.csv"), "live-new")
    # b.csv: tracked at the tag, not moved aside
    _write(str(w / "start" / "d" / "b.csv"), "start-b")
    _write(str(repo / "d" / "b.csv"), "live-b")
    # c.csv: git-ignored, never moved aside
    _write(str(repo / "d" / "c.csv"), "live-c")
    # globbed batches: ignored, live only
    _write(str(repo / "g" / "b000" / "x.csv"), "0")
    _write(str(repo / "g" / "b001" / "x.csv"), "1")
    tracked = {"d/a.csv", "d/b.csv"}
    return str(repo), str(w), tracked


def test_resolver_pre_order(tree):
    repo, w, tracked = tree
    r = D.Resolver("pre", repo, w, tracked=tracked)
    assert open(r.path("d/a.csv")).read() == "pre_repair"
    assert open(r.path("d/b.csv")).read() == "start-b"
    assert open(r.path("d/c.csv")).read() == "live-c"
    assert r.inputs["d/a.csv"]["source"] == "pre_repair"
    assert r.inputs["d/b.csv"]["source"] == "start"
    assert r.inputs["d/c.csv"]["source"] == "live"
    assert "warning" in r.inputs["d/c.csv"] and "warning" not in r.inputs["d/b.csv"]
    assert len(r.inputs["d/a.csv"]["sha256"]) == 64
    got = r.glob("g/b*/x.csv")
    assert [open(p).read() for p in got] == ["0", "1"]
    assert r.inputs["g/b000/x.csv"]["source"] == "live"


def test_resolver_post_reads_live(tree):
    repo, w, tracked = tree
    r = D.Resolver("post", repo, w, tracked=tracked)
    assert open(r.path("d/a.csv")).read() == "live-new"
    assert open(r.path("d/b.csv")).read() == "live-b"
    assert not r.warnings


def test_resolver_missing_input_fails(tree):
    repo, w, tracked = tree
    with pytest.raises(SystemExit):
        D.Resolver("pre", repo, w, tracked=tracked).path("d/missing.csv")


def test_glob_prefers_pre_repair(tree):
    repo, w, tracked = tree
    _write(os.path.join(w, "pre_repair", "g", "b000", "x.csv"), "old0")
    r = D.Resolver("pre", repo, w, tracked=tracked)
    assert [open(p).read() for p in r.glob("g/b*/x.csv")] == ["old0"]


def test_selection_size_stats():
    recs = [("m1", "0", "1"), ("m1", "0", "1"),          # duplicate amplicon row of one md5
            ("m2", "0", "2"), ("m3", "0", "3"), ("m4", "0", "0"),   # m4 abstained
            ("m1", "1", "1")]
    s = D.selection_size_stats(recs)
    assert s["answered"] == 4 and s["single"] == 2
    assert s["fold0_answered"] == 3 and s["fold0_single"] == 1
    assert s["fold0_two"] == 1 and s["fold0_three_plus"] == 1
    assert s["mean"] == pytest.approx((1 + 2 + 3 + 1) / 4)


def test_tie_set_coherence():
    lin = {"a": {"phylum": "P1", "genus": "G1"}, "b": {"phylum": "P1", "genus": "G1"},
           "c": {"phylum": "P2", "genus": "G2"}, "d": {"phylum": "", "genus": None}}
    sets = [("a", "b"), ("a", "b"), ("a", "c"), ("a", "d")]
    t = D.tie_set_coherence(sets, lin)
    assert t["n_records"] == 4 and t["n_unique"] == 3
    assert t["span_records"] == pytest.approx(1 / 4)     # only {a,c}; an empty phylum is not distinct
    assert t["span_unique"] == pytest.approx(1 / 3)
    assert t["genus_unique"] == pytest.approx(1 / 3)     # only {a,b}


def test_pair_mag_amp_compare_arms_rules():
    mag = {"g1": ("ok", True, 0.8), "g2": ("ok", False, None),   # g2: pred_empty -> not answered
           "g3": ("abstained", False, None), "g4": ("truth_empty", False, None),
           "g5": ("ok", True, 0.6)}
    amp = {("g1", "R"): ("ok", 0.7), ("g2", "R"): ("ok", 0.5), ("g3", "R"): ("abstained", None),
           ("g4", "R"): ("ok", 0.9), ("g5", "R"): ("pred_empty", None)}
    s = D.pair_mag_amp(mag, amp, "R")
    assert s["n"] == 1 and s["delta"] == pytest.approx(0.1)
    assert s["n_pop"] == 4                                  # g4 truth_empty dropped
    assert s["pop_mag"] == pytest.approx((0.8 + 0 + 0 + 0.6) / 4)
    assert s["pop_amp"] == pytest.approx((0.7 + 0.5 + 0 + 0) / 4)
    assert s["cov_mag"] == pytest.approx(2 / 4) and s["cov_amp"] == pytest.approx(2 / 4)


def test_exclusion_and_af_sweep():
    hits = {"m1": [dict(ref="100.1", ani=100.0, af_query=99.0),     # self
                   dict(ref="100.7", ani=96.0, af_query=80.0),      # same taxon
                   dict(ref="200.1", ani=97.0, af_query=70.0),      # ANI >= 95
                   dict(ref="300.1", ani=88.0, af_query=10.0)],     # survives, AF 10
            "m2": [dict(ref="400.1", ani=82.0, af_query=40.0)]}     # survives, below ANI floor
    gid = {"m1": "100.1", "m2": "500.1"}
    assert [h["ref"] for h in D.exclude_hits(hits["m1"], "100.1")] == ["300.1"]
    sw = {s["label"]: s["k"] for s in D.af_sweep(["m1", "m2"], hits, gid)}
    labels = [v[0] for v in D.AF_VARIANTS]
    assert sw[labels[0]] == 0          # shipped: m1 fails AF, m2 fails the ANI floor
    assert sw[labels[1]] == 1          # floor 0: m2 qualifies
    assert sw[labels[2]] == 1          # AF 5, floor 83: m1 qualifies (ANI 88, AF 10)
    assert sw[labels[5]] == 2          # any surviving hit


def test_af_sweep_gate_relaxation():
    hits = {"m1": [dict(ref="300.1", ani=88.0, af_query=10.0)]}
    out = D.af_sweep(["m1"], hits, {"m1": "100.1"})
    k = {(s["floor"], s["gate"]): s["k"] for s in out}
    assert k[(83.0, 15.0)] == 0 and k[(83.0, 5.0)] == 1 and k[(83.0, 2.0)] == 1


def test_wald_ci_reproduces_the_printed_interval():
    d, lo, hi = D.wald_diff_ci(2718 / 3942, 3942, 6594 / 9580, 9580)
    assert D._fmt("sf3", d) == "+0.001"
    assert (D._fmt("sf3", lo), D._fmt("sf3", hi)) == ("-0.016", "+0.018")


def test_formatting_and_number():
    assert D._fmt("int", 50842) == "50,842"
    assert D._fmt("pct1", 0.86802) == "86.8"
    assert D._fmt("e1", 2.3469e-4) == "2.3e-04"
    n = D.number("k", 0.8017, "f3", "0.802", "74", None, "r", ["a", "a"])
    assert n["match"] and n["inputs"] == ["a"]
    assert not D.number("k", 0.8012, "f3", "0.802", "74", None, "r", [])["match"]


def test_modal_fraction():
    from collections import Counter
    assert D.modal_fraction(Counter({0.5: 10, 0.333333: 7}))[0] == 0.5
    assert D.modal_fraction(Counter()) == (None, 0)
