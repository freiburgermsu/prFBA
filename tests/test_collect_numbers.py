"""Offline tests for manuscript_numbers/collect.py (W7, the manuscript number differ): selectors,
reductions, number formatting against the printed string, the tex-line check, CSV recipes, the
derived-number lookup, and old (tag) vs new (live) resolution in a throw-away git repository.
Run:
    ~/Documents/py_venv/bin/python -m pytest -q tests/test_collect_numbers.py
"""
import gzip
import json
import os
import subprocess
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "manuscript_numbers"))
import collect as C  # noqa: E402

DOC = {"cells": [{"region": "V4", "regime": "exclude", "est": -0.2905},
                 {"region": "FL", "regime": "exclude", "est": -0.2418},
                 {"region": "FL", "regime": "include", "est": -0.0159}],
       "rows": [{"c": 0.01, "v": 1}, {"c": 0.01, "v": 2}],
       "by": {"a": {"x": 0.70}, "b": {"x": 0.66}, "c": {"x": None}},
       "dist": {"1": 10, "7": 3, "77": 1}}


# ------------------------------------------------------------------ selectors
def test_select_keys_where_each():
    assert C.select(DOC, ["cells", {"where": {"region": "V4", "regime": "exclude"}}, "est"]) == -0.2905
    assert C.select(DOC, ["cells", 1, "est"]) == -0.2418
    assert C.select(DOC, ["by", {"each": ["a", "b"]}, "x"]) == [0.70, 0.66]
    assert C.select(DOC, ["by", {"each": "*"}, "x"]) == [0.70, 0.66, None]
    with pytest.raises(C.SelectorError):
        C.select(DOC, ["cells", {"where": {"region": "FL"}}])      # two matches: must be unique
    with pytest.raises(C.SelectorError):
        C.select(DOC, ["by", "zz"])
    assert C.select(DOC, ["rows", {"where": {"c": 0.01}, "pick": "last"}, "v"]) == 2


def test_reduce():
    vals = C.select(DOC, ["cells", {"each": "*"}, "est"])
    assert C.reduce_values(vals, "min") == -0.2905
    assert C.reduce_values(vals, "minabs") == -0.0159
    assert C.reduce_values(vals, "maxabs") == -0.2905
    assert C.reduce_values(vals, "len") == 3
    assert C.reduce_values(DOC["dist"], "maxkey") == 77
    with pytest.raises(TypeError):
        C.reduce_values([0.1, None], "max")
    assert C.reduce_values([0.1, None], "max", skip_null=True) == 0.1


# ------------------------------------------------------------------ formatting
@pytest.mark.parametrize("v,fmt,out", [
    (0.80176, "f3", "0.802"), (0.10299, "+f3", "+0.103"), (-0.09151, "absf3", "0.092"),
    (0.8893, "pct1", "88.9"), (0.4449, "+pct0", "+44"), (66.383, "pp1", "66.4"),
    (67038, "int", "67,038"), (2.3468e-4, "e1", "2.3e-4"), (1.38e-79, "e0", "1e-79"),
    (-0.0, "f3", "0.000"),
])
def test_fmt_value(v, fmt, out):
    assert C.fmt_value(v, fmt) == out


@pytest.mark.parametrize("v,printed,fmt,ok,why", [
    (0.8017, "0.802", "f3", True, None),
    (67038, "67{,}038", "int", True, None),
    (0.1112, "11.1\\%", "pct1", True, None),
    (0.103, "$+0.103$", "+f3", True, None),
    (0.3075, "0.308", "f3", True, "boundary"),      # 4-dp value on the rounding boundary of the printed one
    (-0.0195, "-0.020", "f3", True, "boundary"),
    (0.4449, "45", "pct0", False, "rounding"),
    (0.3075, "0.3071", "f4", False, "differs"),
    (0.103, "-0.103", "f3", False, "differs"),      # the sign is checked
    (1e-200, "1e-161", "lt", True, None),
    (2e-161, "1e-161", "lt", False, "bound violated"),
    (2.35e-4, "2.3\\times10^{-4}", "e1", True, None),
    (True, "True", "str", True, None),
])
def test_matches_printed(v, printed, fmt, ok, why):
    assert C.matches_printed(v, printed, fmt) == (ok, why)


# ------------------------------------------------------------------ tex check
TEX = ["intro",
       "conditional EC F1 0.802 against 0.798; $+0.103$ and 67{,}038 amplicons",
       "p=$2.3\\times10^{-4}$, every comparison $p < 10^{-161}$, eight amplicons",
       "F10.555 is not a match for 0.555"]


def test_tex_line_set():
    assert C.tex_line_set("420-442; 1175 (57.6-59.2%)", pad=0) == set(range(420, 443)) | {1175}
    assert C.tex_line_set("L77 ('-0.123 ([-0.150]) ...'); 605", pad=1) == {76, 77, 78, 604, 605, 606}


def test_printed_in_tex():
    lines = C.tex_line_set("2-3", pad=0)
    assert C.printed_in_tex(TEX, "0.802", "f3", lines) is True
    assert C.printed_in_tex(TEX, "0.80", "f2", lines) is False          # no partial-number matches
    assert C.printed_in_tex(TEX, "+0.103", "+f3", lines) is True
    assert C.printed_in_tex(TEX, "67,038", "int", lines) is True
    assert C.printed_in_tex(TEX, "2.3e-4", "e1", lines) is True
    assert C.printed_in_tex(TEX, "1e-161", "lt", lines) is True
    assert C.printed_in_tex(TEX, "8", "int", lines, tex_text="eight amplicons") is True
    assert C.printed_in_tex(TEX, "0.798", "f3", {1}) is False          # outside the cited lines
    assert C.printed_in_tex(TEX, "0.555", "f3", {4}) is True           # the second, free-standing one
    assert C.printed_in_tex(TEX, "True", "str", lines) is None


# ------------------------------------------------------------------ old (tag) vs new (live), recipes, derived
def _git(repo, *a):
    subprocess.run(["git", "-C", repo, *a], check=True, capture_output=True,
                   env=dict(os.environ, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t", GIT_COMMITTER_NAME="t",
                            GIT_COMMITTER_EMAIL="t@t"))


@pytest.fixture
def repo(tmp_path):
    r = tmp_path / "repo"
    (r / "data").mkdir(parents=True)
    (r / "manuscript_numbers").mkdir()
    (r / ".gitignore").write_text("data/cells.csv.gz\n")
    (r / "data" / "summary.json").write_text(json.dumps({"all": {"f1": {"value": 0.8017}}, "n": 67038}))
    with gzip.open(r / "data" / "cells.csv.gz", "wt") as fh:
        fh.write("amplicon,method,outcome,n_reps,precision\n"
                 "a__V4,prfba,ok,1,0.9\na__V4,picrust2,ok,,0.5\nb__V4,prfba,abstained,0,\n"
                 "b__V4,picrust2,ok,,0.3\nc__FL,prfba,abstained,2,\nc__FL,picrust2,ok,,0.7\n")
    _git(r, "init", "-q")
    _git(r, "add", ".")
    _git(r, "commit", "-qm", "start")
    _git(r, "tag", "t0")
    # the "re-score": the tracked summary changes and the ignored cell table is moved aside to $W/pre_repair
    (r / "data" / "summary.json").write_text(json.dumps({"all": {"f1": {"value": 0.8101}}, "n": 67713}))
    w = tmp_path / "W"
    (w / "pre_repair" / "data").mkdir(parents=True)
    os.replace(r / "data" / "cells.csv.gz", w / "pre_repair" / "data" / "cells.csv.gz")
    with gzip.open(r / "data" / "cells.csv.gz", "wt") as fh:
        fh.write("amplicon,method,outcome,n_reps,precision\n"
                 "a__V4,prfba,ok,1,0.9\na__V4,picrust2,ok,,0.5\nb__V4,prfba,ok,1,0.8\n"
                 "b__V4,picrust2,ok,,0.3\nc__FL,prfba,abstained,2,\nc__FL,picrust2,ok,,0.7\n")
    (r / "manuscript_numbers" / "derived_numbers.pre.json").write_text(json.dumps(
        {"_meta": {}, "numbers": [{"key": "as10_n", "printed": "59,607", "tex_lines": "521", "value": 59607}]}))
    return str(r), str(w)


def _spec(items):
    return {"_meta": {}, "items": {k: dict({"stage": "S", "inv": "i", "key": k, "tex_lines": "1", "claim": "c"}, **v)
                                   for k, v in items.items()}}


def test_old_vs_new_end_to_end(repo, tmp_path):
    r, w = repo
    abst = {"method": "prfba", "outcome": "abstained"}
    spec = _spec({
        "f1": {"values": [{"label": "F1", "artifact": "data/summary.json", "selector": ["all", "f1", "value"],
                           "printed": "0.802", "fmt": "f3"}]},
        "n": {"values": [{"label": "n", "artifact": "data/summary.json", "selector": ["n"], "printed": "67,038", "fmt": "int"}]},
        "abst": {"values": [{"label": "abstained", "artifact": "data/cells.csv.gz", "recipe": "csv_count",
                             "args": {"filter": abst}, "printed": "2", "fmt": "int"},
                            {"label": "abstained with reps", "artifact": "data/cells.csv.gz", "recipe": "csv_count",
                             "args": {"filter": dict(abst, n_reps=[">", 0])}, "printed": "1", "fmt": "int"},
                            {"label": "regions", "artifact": "data/cells.csv.gz", "recipe": "csv_nunique",
                             "args": {"filter": abst, "column": "reg", "derive": {"reg": ["split", "amplicon", "__", 1]}},
                             "printed": "2", "fmt": "int"},
                            {"label": "PICRUSt2 precision where prFBA abstains", "artifact": "data/cells.csv.gz",
                             "recipe": "csv_mean", "args": {"filter": {"method": "picrust2"}, "column": "precision",
                                                            "within": {"filter": abst, "on": "amplicon"}},
                             "printed": "0.50", "fmt": "f2"}]},
        "share": {"values": [{"label": "1-coverage", "expr": "1-a/b", "printed": "50", "fmt": "pct0",
                              "vars": {"a": {"artifact": "data/cells.csv.gz", "recipe": "csv_count",
                                             "args": {"filter": {"method": "prfba", "outcome": "ok"}}},
                                       "b": {"artifact": "data/cells.csv.gz", "recipe": "csv_count",
                                             "args": {"filter": {"method": "prfba"}}}}}]},
        "derived": {"values": [{"label": "W6 by key", "artifact": "derived:as10_n", "printed": "59,607", "fmt": "int"},
                               {"label": "W6 by printed value", "artifact": "derived:other_name", "printed": "59,607",
                                "fmt": "int"}]},
        "manual": {"values": [], "manual": "no stored source"},
    })
    sp = tmp_path / "spec.json"
    sp.write_text(json.dumps(spec))
    out_j, out_m = tmp_path / "diff.json", tmp_path / "diff.md"
    assert C.main(["--spec", str(sp), "--tag", "t0", "--w", w, "--root", r, "--tex", "",
                   "--out-json", str(out_j), "--out-md", str(out_m), "--quiet"]) == 0
    d = json.loads(out_j.read_text())
    rows = {(o["key"], v["label"]): v for o in d["items"] for v in o["values"]}
    f1 = rows[("f1", "F1")]
    assert (f1["old"], f1["new"], f1["changed"], f1["changed_printed"]) == (0.8017, 0.8101, True, True)
    assert f1["old_matches_printed"] is True
    n = rows[("n", "n")]
    assert (n["old_fmt"], n["new_fmt"]) == ("67,038", "67,713")
    ab = rows[("abst", "abstained")]
    assert (ab["old"], ab["new"]) == (2, 1)           # old from $W/pre_repair, new from the live tree
    assert rows[("abst", "abstained with reps")]["old"] == 1
    assert rows[("abst", "regions")]["old"] == 2
    assert rows[("abst", "PICRUSt2 precision where prFBA abstains")]["old"] == pytest.approx(0.5)
    assert rows[("abst", "PICRUSt2 precision where prFBA abstains")]["new"] == pytest.approx(0.7)
    sh = rows[("share", "1-coverage")]
    assert (sh["old_fmt"], sh["new_fmt"]) == ("67", "33")
    assert sh["old_matches_printed"] is False and sh["mismatch"] == "differs"
    for lab in ("W6 by key", "W6 by printed value"):
        dv = rows[("derived", lab)]
        assert dv["old"] == 59607 and dv["new_status"].startswith("pending")   # no post file yet: never an error
    src = d["_meta"]["sources"]
    assert src["old:data/summary.json"]["source"] in ("tag-git-show", "tag-worktree")
    assert src["old:data/cells.csv.gz"]["source"] == "pre_repair"
    assert src["new:data/cells.csv.gz"]["source"] == "live"
    assert {o["key"]: o["extractor"] for o in d["items"]}["manual"] == "manual"
    assert d["summary"]["n_errors"] == 0
    md = out_m.read_text()
    assert "**CHANGED**" in md and "(manual)" in md


def test_ignored_without_pre_repair_copy_is_flagged(repo, tmp_path):
    r, w = repo
    os.remove(os.path.join(w, "pre_repair", "data", "cells.csv.gz"))
    spec = _spec({"x": {"values": [{"label": "rows", "artifact": "data/cells.csv.gz", "recipe": "csv_count",
                                    "args": {}, "printed": "6", "fmt": "int"}]}})
    sp = tmp_path / "spec.json"
    sp.write_text(json.dumps(spec))
    out_j = tmp_path / "d.json"
    C.main(["--spec", str(sp), "--tag", "t0", "--w", w, "--root", r, "--tex", "", "--out-json", str(out_j),
            "--out-md", str(tmp_path / "d.md"), "--quiet"])
    d = json.loads(out_j.read_text())
    assert d["_meta"]["ignored_read_live_as_old"] == ["data/cells.csv.gz"]


def test_committed_spec_covers_every_dependent_item():
    spec = json.load(open(os.path.join(ROOT, "manuscript_numbers", "claims_spec.json")))
    items = spec["items"]
    assert len(items) == 178
    for iid, it in items.items():
        assert it["values"] or it.get("manual"), f"{iid}: neither extractors nor a manual reason"
        for v in it["values"]:
            assert v.get("printed") is not None and v.get("fmt"), iid
            assert ("artifact" in v) or ("expr" in v), iid
