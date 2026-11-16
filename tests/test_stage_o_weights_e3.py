"""region_validation/shared_reference/scripts/stage_o_weights_e3.py: Stage O's E3 weights file,
derived from Stage E's stage_o_weights.json (isotonic map removed, estimator relabelled)."""
import hashlib
import json
import os
import sys

import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = os.path.join(ROOT, "region_validation", "shared_reference", "scripts")
PFE = os.path.join(ROOT, "region_validation", "shared_reference", "data", "pfunc_estimators")
sys.path.insert(0, SCRIPTS)
import stage_o_weights_e3 as S  # noqa: E402

ISO = {"x": [0.0, 0.5, 1.0], "y": [0.01, 0.48, 0.943375]}


def _ci(est, lo, hi):
    return {"est": est, "lo": lo, "hi": hi}


def _fixture(tmp_path, winner="E4(E3)", label="E4", iso=ISO, sha="abc", K=100, uni=None):
    amps = {"a__V4": {"md5": "m1", "fold": 0, "neighbours": [["g1", 1.0], ["g2", 0.5]],
                      "clade": {"rank": "family", "members": ["g1", "g3"]}, "alpha": 2.0,
                      "n_reps": 1, "reps": ["g1"]},
            "b__V1-V3": {"md5": "m2", "fold": 3, "neighbours": [["g4", 1.0]], "clade": None,
                         "alpha": 0.0, "n_reps": 1, "reps": ["g4"]}}
    src = {"meta": {"estimator": label, "K": K, "tau": 0.02, "dedup": "genome",
                    "clade_rank_rule": "family -> order; ...", "alpha": 2.0, "isotonic": iso,
                    "script_sha256": sha, "generated": "2026-11-02T23:06:05-0600"},
           "amplicons": amps}
    pc = {"universe_brier": uni or _ci(-0.0002568894, -0.0007098106, 0.0001442025),
          "union_brier": _ci(-0.0005459983, -0.002385615, 0.00101706),
          "macro_F1_q_ge_0.5": _ci(-0.004248453, -0.006655402, -0.001694771)}
    scores = {"meta": {"script_sha256": "abc"},
              "winner": {"name": winner, "base_family": winner[3:5] if winner.startswith("E4(") else winner,
                         "carried_forward_variant": "E3|K=100|tau=0.02|genome|family|alpha=2",
                         "isotonic_all_folds": iso},
              "composites": {"d100": {"all": {"contrasts_pairwise": {"E4(E3) - E3": pc}}}}}
    p_src, p_sc = tmp_path / "stage_o_weights.json", tmp_path / "scores.json"
    p_src.write_text(json.dumps(src, indent=2))
    p_sc.write_text(json.dumps(scores))
    return p_src, p_sc, src


def test_derivation(tmp_path):
    p_src, p_sc, src = _fixture(tmp_path)
    out = tmp_path / "E3.json"
    rc = S.main(["--src", str(p_src), "--scores", str(p_sc), "--out", str(out),
                 "--generated-e3", "2026-11-03T08:50:15-06:00"])
    assert rc == 0
    raw = out.read_bytes()
    obj = json.loads(raw)
    assert raw == json.dumps(obj).encode()                     # compact, no trailing newline
    m = obj["meta"]
    assert list(m) == ["estimator", "K", "tau", "dedup", "clade_rank_rule", "alpha", "isotonic",
                       "script_sha256", "generated", "derived_from", "generated_e3"]
    assert m["estimator"] == "E3" and m["isotonic"] is None
    assert {k: m[k] for k in ("K", "tau", "dedup", "alpha", "script_sha256", "generated")} == \
        {k: src["meta"][k] for k in ("K", "tau", "dedup", "alpha", "script_sha256", "generated")}
    assert obj["amplicons"] == src["amplicons"] and list(obj["amplicons"]) == list(src["amplicons"])
    d = m["derived_from"]
    assert d["file"] == "stage_o_weights.json"
    assert d["sha256"] == hashlib.sha256(p_src.read_bytes()).hexdigest()
    assert d["change"] == (
        "isotonic map removed: Stage O carries E3, not E4(E3). E4 vs E3 universe Brier -0.0003 "
        "[-0.0007, 0.0001] and union Brier -0.0005 [-0.0024, 0.0010] cover 0, macro F1 -0.0042 "
        "[-0.0067, -0.0017]; the map is fitted on all five folds (in-sample for these 2 "
        "amplicons), on per-EC values (applied at role level it extrapolates), and caps p at "
        "0.943 so no reaction reaches the cost floor.")
    assert m["generated_e3"] == "2026-11-03T08:50:15-06:00"
    # --expect-equal: equal bytes -> 0; anything else -> 1
    assert S.main(["--src", str(p_src), "--scores", str(p_sc), "--out", str(tmp_path / "b.json"),
                   "--generated-e3", "2026-11-03T08:50:15-06:00", "--expect-equal", str(out)]) == 0
    assert S.main(["--src", str(p_src), "--scores", str(p_sc), "--out", str(tmp_path / "c.json"),
                   "--generated-e3", "2026-11-04T00:00:00-06:00", "--expect-equal", str(out)]) == 1


def test_change_text_follows_the_numbers(tmp_path):
    p_src, p_sc, _ = _fixture(tmp_path, uni=_ci(-0.0011, -0.0019, -0.0004))
    out = tmp_path / "E3.json"
    S.main(["--src", str(p_src), "--scores", str(p_sc), "--out", str(out)])
    ch = json.loads(out.read_bytes())["meta"]["derived_from"]["change"]
    assert "universe Brier -0.0011 [-0.0019, -0.0004] (excludes 0) and union Brier" in ch
    assert "(covers 0)" in ch


@pytest.mark.parametrize("kw,msg", [
    (dict(winner="E1", label="E1", iso=None), "not E3-based"),
    (dict(label="E3"), "meta.estimator"),
    (dict(sha="other"), "different pfunc_estimators.py"),
    (dict(K=50), "meta.K"),
])
def test_refusals(tmp_path, kw, msg):
    p_src, p_sc, _ = _fixture(tmp_path, **kw)
    with pytest.raises(SystemExit) as e:
        S.main(["--src", str(p_src), "--scores", str(p_sc), "--out", str(tmp_path / "x.json")])
    assert msg in str(e.value)


def test_isotonic_mismatch_refused(tmp_path):
    p_src, p_sc, _ = _fixture(tmp_path)
    sc = json.loads(p_sc.read_text())
    sc["winner"]["isotonic_all_folds"] = {"x": [0.0, 1.0], "y": [0.0, 0.9]}
    p_sc.write_text(json.dumps(sc))
    with pytest.raises(SystemExit, match="isotonic"):
        S.main(["--src", str(p_src), "--scores", str(p_sc), "--out", str(tmp_path / "x.json")])


@pytest.mark.skipif(not all(os.path.exists(os.path.join(PFE, f)) for f in
                            ("stage_o_weights.json", "scores.json", "stage_o_weights_E3.json")),
                    reason="Stage E outputs not present")
def test_committed_files_are_consistent(tmp_path):
    """The committed E3 file is exactly this script's output from the committed Stage E files."""
    ref = os.path.join(PFE, "stage_o_weights_E3.json")
    ts = json.load(open(ref))["meta"]["generated_e3"]
    rc = S.main(["--src", os.path.join(PFE, "stage_o_weights.json"),
                 "--scores", os.path.join(PFE, "scores.json"), "--out", str(tmp_path / "E3.json"),
                 "--generated-e3", ts, "--expect-equal", ref])
    assert rc == 0
