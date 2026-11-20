"""Tests of the re-score tools in this folder (ctl_run.py, build_snapshot.py, attrib_diff.py).

    ~/Documents/py_venv/bin/python -m pytest -q /home/freiburger/Documents/prFBA_rescore_2026-11/tools/test_tools.py
"""
import json
import os
import subprocess
import sys
import textwrap

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
REAL_CACHE = "/home/freiburger/Documents/prFBA/bvbrc_cache/genome_gene_families.json"


def _run(args, **kw):
    return subprocess.run([PY] + args, capture_output=True, text=True, **kw)


# ----------------------------------------------------------------------------- ctl_run
@pytest.fixture()
def tree(tmp_path):
    sd = tmp_path / "scripts"
    sd.mkdir()
    (sd / "_common.py").write_text(textwrap.dedent(f"""
        class PATHS:
            pgfam_cache = {REAL_CACHE!r}
            out_json = {str(tmp_path / 'committed.json')!r}
    """))
    (sd / "helper.py").write_text(textwrap.dedent(f"""
        CACHE = {REAL_CACHE!r}
        def first_bytes(path={REAL_CACHE!r}):
            return open(path).read(20)
    """))
    snap = tmp_path / "snap.json"
    snap.write_text('{"g1": ["PGF_1", "PGF_2"]}')
    return tmp_path, sd, snap


def test_ctl_run_redirects_cache_defaults_and_outputs(tree):
    tmp, sd, snap = tree
    out = tmp / "ctl_out.json"
    (sd / "job.py").write_text(textwrap.dedent("""
        import json, sys, multiprocessing as mp
        import _common as C
        import helper
        OUT = C.PATHS.out_json
        def work(x):
            return helper.first_bytes()[:5] + str(x)
        def main():
            with mp.get_context("fork").Pool(2) as p:
                r = p.map(work, [1, 2])
            json.dump({"cache": C.PATHS.pgfam_cache, "helper": helper.CACHE, "r": r,
                       "argv": sys.argv[1:]}, open(OUT, "w"))
        if __name__ == "__main__":
            main()
    """))
    rec = tmp / "rec.json"
    p = _run([os.path.join(HERE, "ctl_run.py"), "--common", str(sd), "--cache", str(snap),
              "--redirect", f"PATHS.out_json={out}", "--record", str(rec),
              "--", str(sd / "job.py"), "--x", "1"])
    assert p.returncode == 0, p.stderr
    d = json.load(open(out))
    assert d["cache"] == str(snap) and d["helper"] == str(snap)
    assert d["r"] == ['{"g1"1', '{"g1"2']          # forked workers read the snapshot
    assert d["argv"] == ["--x", "1"]
    assert not (tmp / "committed.json").exists()
    r = json.load(open(rec))
    assert r["rc"] == 0 and r["forbidden_opens"] == []
    assert any(x["where"].endswith("helper.first_bytes") for x in r["rebound"])


def test_ctl_run_forbids_the_live_cache(tree):
    tmp, sd, snap = tree
    (sd / "sneaky.py").write_text(
        "import os\nopen(os.path.join('/home/freiburger/Documents/prFBA', 'bvbrc_cache', "
        "'genome_gene_families.json')).read(5)\n")
    p = _run([os.path.join(HERE, "ctl_run.py"), "--common", str(sd), "--cache", str(snap),
              "--", str(sd / "sneaky.py")])
    assert p.returncode != 0
    assert "forbidden" in p.stderr


def test_ctl_run_forbids_snapshot_writes_unless_allowed(tree):
    tmp, sd, snap = tree
    (sd / "writer.py").write_text(
        "import _common as C\nopen(C.PATHS.pgfam_cache + '.tmp', 'w').write('x')\n")
    p = _run([os.path.join(HERE, "ctl_run.py"), "--common", str(sd), "--cache", str(snap),
              "--", str(sd / "writer.py")])
    assert p.returncode != 0 and "forbidden" in p.stderr
    p = _run([os.path.join(HERE, "ctl_run.py"), "--common", str(sd), "--cache", str(snap),
              "--allow-cache-writes", "--", str(sd / "writer.py")])
    assert p.returncode == 0, p.stderr
    # the atomic-write pattern: a differently named temp file os.replace()d over the snapshot
    other = tmp / "elsewhere.tmp"
    (sd / "replacer.py").write_text(
        "import os, _common as C\nopen(%r, 'w').write('{}')\nos.replace(%r, C.PATHS.pgfam_cache)\n"
        % (str(other), str(other)))
    before = snap.read_text()
    p = _run([os.path.join(HERE, "ctl_run.py"), "--common", str(sd), "--cache", str(snap),
              "--", str(sd / "replacer.py")])
    assert p.returncode != 0 and "os.rename" in p.stderr and snap.read_text() == before
    (sd / "remover.py").write_text("import os, _common as C\nos.remove(C.PATHS.pgfam_cache)\n")
    p = _run([os.path.join(HERE, "ctl_run.py"), "--common", str(sd), "--cache", str(snap),
              "--", str(sd / "remover.py")])
    assert p.returncode != 0 and snap.exists()


def test_ctl_run_entry_mode_main_redirect(tree):
    tmp, sd, snap = tree
    out = tmp / "entry_out.txt"
    (sd / "consts.py").write_text(textwrap.dedent(f"""
        import multiprocessing as mp
        OUT = {str(tmp / 'committed.txt')!r}
        def sq(x):
            return x * x
        def main():
            with mp.get_context("fork").Pool(2) as p:
                r = p.map(sq, [2, 3])
            open(OUT, "w").write(str(r))
        if __name__ == "__main__":
            main()
    """))
    p = _run([os.path.join(HERE, "ctl_run.py"), "--entry", "main",
              "--redirect", f"main:OUT={out}", "--", str(sd / "consts.py")])
    assert p.returncode == 0, p.stderr
    assert out.read_text() == "[4, 9]"
    assert not (tmp / "committed.txt").exists()
    # a JSON-valued redirect (dict of paths, as the E3 gate's REF_CELLS)
    (sd / "dictjob.py").write_text("import json\nREF = {'a': 'x'}\nOUT = 'nowhere'\n"
                                   "def main():\n    json.dump(REF, open(OUT, 'w'))\n")
    o2 = tmp / "dict_out.json"
    p = _run([os.path.join(HERE, "ctl_run.py"), "--entry", "main", "--redirect", f"main:OUT={o2}",
              "--redirect", 'main:REF={"a": "/w/start/cells.csv"}', "--", str(sd / "dictjob.py")])
    assert p.returncode == 0, p.stderr
    assert json.load(open(o2)) == {"a": "/w/start/cells.csv"}


def test_ctl_run_rebinds_pfunc_estimator_defaults(tmp_path):
    snap = tmp_path / "snap.json"
    snap.write_text("{}")
    cat = tmp_path / "cat.json"
    cat.write_text("{}")
    job = tmp_path / "job.py"
    job.write_text(textwrap.dedent("""
        import sys
        sys.path.insert(0, "/home/freiburger/Documents/prFBA")
        import pfunc_estimator as P
        print("CAT", P.load_clade_catalogue.__defaults__[0])
        print("CUR", P.catalogue_is_current.__defaults__[0])
    """))
    p = _run([os.path.join(HERE, "ctl_run.py"), "--cache", str(snap), "--catalogue", str(cat),
              "--", str(job)])
    assert p.returncode == 0, p.stderr
    assert f"CAT {cat}" in p.stdout and f"CUR {snap}" in p.stdout


# ----------------------------------------------------------------------------- build_snapshot
def test_build_snapshot_inverts_a_repair(tmp_path):
    # a repaired cache in the exact serialisation of the shared cache (json.dumps defaults)
    old = {"a": ["PGF_1", "PGF_2"], "b": [], "c": ["PGF_9"], "d": ["PGF_3"]}
    new = dict(old, b=["PGF_7", "PGF_8"], c=["PGF_9", "PGF_10"])
    oldf, newf = tmp_path / "old.json", tmp_path / "cache.json"
    oldf.write_text(json.dumps(old))
    newf.write_text(json.dumps(new))
    import hashlib
    sha_old = hashlib.sha256(oldf.read_bytes()).hexdigest()
    corr = {"files": {str(newf): {"before": {"size": oldf.stat().st_size, "mtime": 1794418222.7086723,
                                             "sha256": sha_old},
                                  "after": {"size": newf.stat().st_size, "sha256": "x"},
                                  "n_entries_replaced": 2}},
            "entries": [{"file": str(newf), "status": "applied", "reason": "empty_list_poison", "genome": "b"},
                        {"file": str(newf), "status": "applied", "reason": "truncated_by_patric_cap",
                         "genome": "c"}],
            "backup_old_values": {str(newf): {"b": [], "c": ["PGF_9"]}}}
    cf = tmp_path / "corr.json"
    cf.write_text(json.dumps(corr))
    out = tmp_path / "snap.json"
    p = _run([os.path.join(HERE, "build_snapshot.py"), "--cache", str(newf), "--corrections", str(cf),
              "--out", str(out), "--log-key", str(newf)])
    assert p.returncode == 0, p.stdout + p.stderr
    assert out.read_bytes() == oldf.read_bytes()
    assert os.stat(out).st_mtime == 1794418222.7086723
    side = json.load(open(str(out) + ".stamp.json"))
    assert side["sha256_equals_target"] is True and side["entry_level_check"]["ok"]
    # the d20 variant: also reset 'd' to []
    rk = tmp_path / "reset.json"
    rk.write_text('["d"]')
    out20 = tmp_path / "snap20.json"
    p = _run([os.path.join(HERE, "build_snapshot.py"), "--cache", str(newf), "--corrections", str(cf),
              "--out", str(out20), "--log-key", str(newf), "--d20", "--reset-keys", str(rk)])
    assert p.returncode == 0, p.stdout + p.stderr
    assert json.load(open(out20)) == dict(old, d=[])


# ----------------------------------------------------------------------------- attrib_diff
def _write_csv(path, header, rows):
    path.write_text(",".join(header) + "\n" + "\n".join(",".join(map(str, r)) for r in rows) + "\n")


def test_attrib_diff_mag(tmp_path):
    rep_ids = tmp_path / "rep.json"
    rep_ids.write_text(json.dumps({"repaired_2026_11_13": ["R1", "R2"], "poisoned_d20": ["P1"]}))
    h = ["mag", "genome_id", "completeness", "contamination", "regime", "outcome", "ani_class", "tier",
         "best_ani", "n_hits", "n_selected", "anchor", "recall", "precision", "f1", "jaccard",
         "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    base = [["m1", "R1", 1, 0, "exclude", "ok", ">=95", "s", 99, 3, 1, "A", .5, .5, .5, .3] + [1] * 6,
            ["m2", "G2", 1, 0, "exclude", "ok", ">=95", "s", 99, 3, 1, "R2", .5, .5, .5, .3] + [1] * 6,
            ["m3", "G3", 1, 0, "exclude", "ok", ">=95", "s", 99, 3, 2, "A", .5, .5, .5, .3] + [1] * 6,
            ["m4", "G4", 1, 0, "exclude", "ok", ">=95", "s", 99, 3, 1, "A", .5, .5, .5, .3] + [1] * 6,
            ["m5", "P1", 1, 0, "exclude", "ok", ">=95", "s", 99, 3, 1, "A", .5, .5, .5, .3] + [1] * 6]
    new = [list(r) for r in base]
    for r in new[:3]:
        r[14] = .9                                             # source / anchor / multi-pick changes
    old_p, new_p = tmp_path / "old.csv", tmp_path / "new.csv"
    _write_csv(old_p, h, base)
    _write_csv(new_p, h, new)
    out = tmp_path / "r.json"
    p = _run([os.path.join(HERE, "attrib_diff.py"), "--table", "mag", "--old", str(old_p),
              "--new", str(new_p), "--repaired", str(rep_ids), "--out", str(out)])
    assert p.returncode == 0, p.stdout + p.stderr
    r = json.load(open(out))
    assert r["by_reason_first"] == {"anchor_repaired": 1, "multi_pick_unstored": 1, "source_repaired": 1}
    # an unrepaired single-pick change fails; with --d20 the P1 change is attributed
    new[3][14] = .9
    new[4][14] = .9
    _write_csv(new_p, h, new)
    p = _run([os.path.join(HERE, "attrib_diff.py"), "--table", "mag", "--old", str(old_p),
              "--new", str(new_p), "--repaired", str(rep_ids), "--out", str(out)])
    assert p.returncode == 1
    r = json.load(open(out))
    assert r["n_unattributed"] == 2 and r["unattributed_examples"][0]["mag"] == "m4"
    p = _run([os.path.join(HERE, "attrib_diff.py"), "--table", "mag", "--old", str(old_p),
              "--new", str(new_p), "--repaired", str(rep_ids), "--d20", "--out", str(out)])
    assert json.load(open(out))["n_unattributed"] == 1
    # an identity column change is a failure on its own
    new = [list(r) for r in base]
    new[0][8] = 98
    _write_csv(new_p, h, new)
    p = _run([os.path.join(HERE, "attrib_diff.py"), "--table", "mag", "--old", str(old_p),
              "--new", str(new_p), "--repaired", str(rep_ids), "--out", str(out)])
    assert p.returncode == 1 and json.load(open(out))["n_identity_changed"] == 1


def test_attrib_diff_genecap_with_reps_pickle(tmp_path):
    import pickle
    rep_ids = tmp_path / "rep.json"
    rep_ids.write_text(json.dumps({"repaired_2026_11_13": ["R9"], "poisoned_d20": []}))
    h = ["ampliconKey", "src_genome_id", "src_taxon_id", "region", "domain", "tier", "self_mode",
         "n_selected", "n_reps", "best_identity", "ident_bin", "best_nonself_identity", "truth_size",
         "pred_size", "f1", "outcome"]
    base = [["G1__V4", "G1", 1, "V4", "Bacteria", "s", "exclude_0.987", 1, 1, 1, "x", .99, 10, 10, .5, "ok"],
            ["G2__V4", "G2", 2, "V4", "Bacteria", "s", "exclude_0.987", 1, 1, 1, "x", .99, 10, 10, .5, "ok"]]
    new = [list(r) for r in base]
    new[0][13], new[1][13] = 12, 12
    old_p, new_p = tmp_path / "o.csv", tmp_path / "n.csv"
    _write_csv(old_p, h, base)
    _write_csv(new_p, h, new)
    pk = tmp_path / "reps.pkl"
    pickle.dump({"cells": {"G1__V4": [{"self_mode": "exclude_0.987", "_reps": ("R9",)}],
                           "G2__V4": [{"self_mode": "exclude_0.987", "_reps": ("X",)}]}}, open(pk, "wb"))
    out = tmp_path / "r.json"
    p = _run([os.path.join(HERE, "attrib_diff.py"), "--table", "genecap", "--old", str(old_p),
              "--new", str(new_p), "--repaired", str(rep_ids), "--reps-pickle", str(pk), "--out", str(out)])
    r = json.load(open(out))
    assert p.returncode == 1 and r["by_reason_first"] == {"rep_repaired": 1, "unattributed": 1}
    p = _run([os.path.join(HERE, "attrib_diff.py"), "--table", "genecap", "--old", str(old_p),
              "--new", str(new_p), "--repaired", str(rep_ids), "--out", str(out)])
    assert json.load(open(out))["by_reason_first"] == {"unattributed_reps_unknown": 2}


def test_ctl_run_does_not_rebind_itself(tmp_path):
    """A script in ctl_run's own directory (e.g. dump_genecap_reps.py) puts that directory among
    the rebinding roots; ctl_run's own REAL_CACHE must stay the live path (it is stamped)."""
    import shutil
    shutil.copy(os.path.join(HERE, "ctl_run.py"), tmp_path / "ctl_run.py")
    snap = tmp_path / "snap.json"
    snap.write_text("{}")
    (tmp_path / "job.py").write_text("print('ran')\n")
    rec = tmp_path / "rec.json"
    p = _run([str(tmp_path / "ctl_run.py"), "--cache", str(snap), "--record", str(rec),
              "--", str(tmp_path / "job.py")])
    assert p.returncode == 0, p.stderr
    r = json.load(open(rec))
    assert r["real_cache_before"] == r["real_cache_after"] and r["rc"] == 0
    assert not any("ctl_run" in x["where"] or x["where"].startswith("__main__") for x in r["rebound"])
