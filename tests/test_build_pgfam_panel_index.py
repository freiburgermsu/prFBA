"""mag_validation/build_pgfam_panel_index.py (re-score W1): rebuild of the MAG-arm panel PGFam
index from the shared cache, on a tiny synthetic cache and corrections log."""
import json
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "mag_validation"))
import build_pgfam_panel_index as B  # noqa: E402


def _fixture(tmp_path):
    # pre-repair values; "g2" was an empty-list poison entry and "g3" a truncated one
    old_vals = {"g3": ["PGF_5", "PGF_1"], "g1": ["PGF_2", "PGF_1", "PGF_2"], "g2": [],
                "g4": ["PGF_9", "PGF_1"]}
    new_vals = dict(old_vals, g2=["PGF_7", "PGF_1"], g3=["PGF_5", "PGF_1", "PGF_8"])
    cache = dict(new_vals, x1=["PGF_100"])                     # non-panel genome
    cpath = tmp_path / "cache.json"
    cpath.write_text(json.dumps(cache))
    key = os.path.realpath(cpath)
    corr = {"entries": [
        {"genome": "g2", "file": key, "status": "applied", "reason": "empty_list_poison",
         "old_size": 0, "new_size": 2},
        {"genome": "g3", "file": key, "status": "applied", "reason": "truncated_by_patric_cap",
         "old_size": 2, "new_size": 3}],
        "backup_old_values": {key: {"g2": [], "g3": ["PGF_5", "PGF_1"]}}}
    corr_path = tmp_path / "corr.json"
    corr_path.write_text(json.dumps(corr))
    # the old index, in the original builder's encoding (order g3, g1, g2, g4)
    order = ["g3", "g1", "g2", "g4"]
    n_vocab, genomes = B.encode(order, old_vals)
    old_path = tmp_path / "old_index.json"
    old_path.write_text(json.dumps({"n_vocab": n_vocab, "genomes": genomes}))
    return cpath, corr_path, old_path, order, new_vals


def test_encode_first_appearance_sorted_unique():
    n, g = B.encode(["b", "a"], {"a": ["PGF_3", "PGF_1"], "b": ["PGF_2", "PGF_1", "PGF_2"]})
    assert n == 3
    assert g == {"b": [0, 1], "a": [0, 2]}      # PGF_1 -> 0, PGF_2 -> 1 (from b), PGF_3 -> 2


def test_build_new_index_and_checks(tmp_path):
    cpath, corr_path, old_path, order, new_vals = _fixture(tmp_path)
    out = tmp_path / "new.json"
    side = B.build(str(old_path), str(out), str(cpath), str(corr_path))
    d = json.load(open(out))
    assert list(d["genomes"]) == order
    assert {g: len(v) for g, v in d["genomes"].items()} == {g: len(set(new_vals[g])) for g in order}
    assert side["len_differs_from_old"] == ["g2", "g3"]
    assert side["cache_before"]["sha256"] == side["cache_after"]["sha256"]
    ok, rep = B.check(str(old_path), str(out), str(cpath), str(corr_path), n_pairs=50)
    assert ok, rep
    assert rep["n_repaired_in_panel"] == 2 and rep["n_len_eq_old_among_repaired"] == 0


def test_restore_old_values_reproduces_the_old_index_bytes(tmp_path):
    cpath, corr_path, old_path, order, _ = _fixture(tmp_path)
    out = tmp_path / "ctl.json"
    side = B.build(str(old_path), str(out), str(cpath), str(corr_path), restore_old=True)
    assert side["n_restored"] == 2
    assert out.read_bytes() == old_path.read_bytes()
    ok, rep = B.check(str(old_path), str(out), str(cpath), str(corr_path), restore_old=True, n_pairs=50)
    assert ok, rep


def test_check_catches_a_wrong_index(tmp_path):
    cpath, corr_path, old_path, order, _ = _fixture(tmp_path)
    bad = json.load(open(old_path))
    bad["genomes"]["g1"] = bad["genomes"]["g1"][:1]              # an unrepaired genome changed
    bpath = tmp_path / "bad.json"
    bpath.write_text(json.dumps(bad))
    ok, rep = B.check(str(old_path), str(bpath), str(cpath), str(corr_path), n_pairs=50)
    assert not ok
    assert rep["len_eq_cache_set"] is False and rep["len_eq_old_outside_repaired"] is False
