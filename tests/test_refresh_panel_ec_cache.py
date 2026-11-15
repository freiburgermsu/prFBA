"""region_validation/shared_reference/scripts/refresh_panel_ec_cache.py (re-score W3): Path A
refresh of the answer-key EC cache, on a tiny synthetic data dir."""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "region_validation", "shared_reference", "scripts"))
import refresh_panel_ec_cache as RP  # noqa: E402


def _setup(tmp_path, gfs):
    d = tmp_path / "data"
    d.mkdir()
    fn = tmp_path / "gfs.json"
    fn.write_text(json.dumps(gfs))
    (d / "panel_eligible.json").write_text(json.dumps(
        {"meta": {}, "genomes": {"g1": {}, "g2": {}}, "dropped": {"g3": "empty_ec", "g4": "no_function_set"}}))
    (d / "panel_ec_universe.json").write_text(json.dumps({"meta": {}, "ec": ["1.1.1.1", "2.7.1.1"]}))
    return d, fn


OLD_GFS = {"g1": {"ec": ["EC:1.1.1.1", "EC:2.7.-.-"]}, "x": {"ec": ["EC:9.9.9.9"]},
           "g2": {"ec": ["EC:2.7.1.1"]}, "g3": {"ec": []}}


def test_control_reproduces_the_cache_and_refresh_reports(tmp_path):
    d, fn = _setup(tmp_path, OLD_GFS)
    # the committed cache, as plan_folds builds it
    import plan_folds as PF
    PF._cache_panel_ec(str(fn), {"g1", "g2", "g3", "g4"}, str(d / "_cache_panel_ec.json"), force=True)
    committed = (d / "_cache_panel_ec.json").read_bytes()
    out = tmp_path / "ctl.json"
    rc = RP.main(["--functions", str(fn), "--data-dir", str(d), "--out", str(out),
                  "--expect-equal", str(d / "_cache_panel_ec.json")])
    assert rc == 0
    assert out.read_bytes() == committed
    rep = json.load(open(str(out) + ".refresh.json"))
    assert rep["n_changed"] == 0 and rep["guarded_unchanged"] and rep["expect_equal"]["equal_bytes"]

    # the regenerated function file: g3 gains ECs (empty -> repaired), g4 now has a function set
    new_gfs = dict(OLD_GFS, g3={"ec": ["EC:1.1.1.1", "EC:3.1.1.1"]}, g4={"ec": ["EC:1.1.1.1"]})
    fn.write_text(json.dumps(new_gfs))
    rc = RP.main(["--functions", str(fn), "--data-dir", str(d)])
    assert rc == 0
    new = json.load(open(d / "_cache_panel_ec.json"))
    assert set(new) == {"g1", "g2", "g3"}                   # key set unchanged: g4 stays out
    assert new["g3"]["ec"] == ["1.1.1.1", "3.1.1.1"] and new["g1"]["ec"] == ["1.1.1.1"]
    rep = json.load(open(str(d / "_cache_panel_ec.json") + ".refresh.json"))
    assert rep["n_changed"] == 1 and rep["changed"]["g3"]["gained"] == 2
    assert rep["no_function_set_now_covered"] == ["g4"]
    assert rep["guarded_unchanged"]


def test_missing_genome_aborts(tmp_path):
    d, fn = _setup(tmp_path, OLD_GFS)
    import plan_folds as PF
    PF._cache_panel_ec(str(fn), {"g1", "g2", "g3"}, str(d / "_cache_panel_ec.json"), force=True)
    fn.write_text(json.dumps({k: v for k, v in OLD_GFS.items() if k != "g2"}))
    try:
        RP.main(["--functions", str(fn), "--data-dir", str(d), "--dry-run"])
    except SystemExit as e:
        assert "no function set" in str(e)
    else:
        raise AssertionError("expected an abort")
