#!/usr/bin/env python
"""One-off refresh of the shared-reference answer-key EC cache (`_cache_panel_ec.json`) under the
pre-registered panel (2026-11 re-score, Path A; plan item W3).

    python refresh_panel_ec_cache.py --functions <genome_function_sets.json> --data-dir <SR data dir>
                                     [--out PATH] [--expect-equal PATH] [--dry-run]

`plan_folds.py --phase panel` builds `_cache_panel_ec.json` (`{genome_id: {"ec": [...],
"n_ec_raw": n}}` for the panel genomes that have a function set) and re-decides eligibility at the
same time. Path A keeps the panel, the folds and eligibility (`panel_eligible.json`,
`panel_ec_universe.json`) as pre-registered and refreshes only the EC sets, from the regenerated
d20 `genome_function_sets.json`. This script does exactly that and nothing else:

* the genome set is the EXISTING cache's key set (asserted unchanged afterwards); the projection
  is `plan_folds._cache_panel_ec` itself (same norm_ec, same sort, same n_ec_raw), called with
  force=True into a temp file in the same directory, then renamed over the target;
* `panel_eligible.json` and `panel_ec_universe.json` in the data dir must keep their sha256 and
  inode (in `data_family/` they are hard links of `data/`'s), checked before and after;
* it reports which genomes changed (ECs gained / lost), whether they are eligible, any new EC
  outside the frozen EC universe, and panel genomes dropped as `no_function_set` that the new
  function file now covers (eligibility is not re-decided: they stay dropped under Path A).

Control: run with the OLD function file and `--expect-equal <committed cache>`; the content must
be equal (it is then also written to --out, never over the committed file unless --out says so).
"""
import argparse
import hashlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import plan_folds as PF  # noqa: E402

GUARDED = ("panel_eligible.json", "panel_ec_universe.json")


def sha256_file(path, chunk=16 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def guard_stamp(data_dir):
    out = {}
    for name in GUARDED:
        p = os.path.join(data_dir, name)
        st = os.stat(p)
        out[name] = {"sha256": sha256_file(p), "inode": st.st_ino, "size": st.st_size, "mtime": st.st_mtime}
    return out


def diff_caches(old, new, eligible, universe):
    changed = {}
    for g in old:
        a, b = set(old[g]["ec"]), set(new[g]["ec"])
        if a != b or old[g].get("n_ec_raw") != new[g].get("n_ec_raw"):
            changed[g] = {"eligible": g in eligible, "n_ec_old": len(a), "n_ec_new": len(b),
                          "gained": len(b - a), "lost": len(a - b),
                          "n_ec_raw_old": old[g].get("n_ec_raw"), "n_ec_raw_new": new[g].get("n_ec_raw"),
                          "new_outside_universe": sorted((b - a) - universe)}
    return changed


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--functions", required=True, help="genome_function_sets.json to project from")
    ap.add_argument("--data-dir", required=True, help="SR data dir holding _cache_panel_ec.json")
    ap.add_argument("--out", help="where to write (default: <data-dir>/_cache_panel_ec.json)")
    ap.add_argument("--expect-equal", help="control: the result must equal this cache in content")
    ap.add_argument("--report", help="write the JSON report here (default: <out>.refresh.json)")
    ap.add_argument("--dry-run", action="store_true", help="compute and report, write nothing")
    a = ap.parse_args(argv)

    t0 = time.time()
    target = os.path.join(a.data_dir, "_cache_panel_ec.json")
    out = a.out or target
    old = PF.load_json(target)
    keys = list(old)
    eligible_doc = PF.load_json(os.path.join(a.data_dir, "panel_eligible.json"))
    eligible = set(eligible_doc["genomes"])
    nofs = {g for g, r in eligible_doc.get("dropped", {}).items() if r == "no_function_set"}
    universe = set(PF.load_json(os.path.join(a.data_dir, "panel_ec_universe.json"))["ec"])
    g0 = guard_stamp(a.data_dir)
    fn_stamp = {"path": os.path.abspath(a.functions), "size": os.path.getsize(a.functions),
                "sha256": sha256_file(a.functions)}

    # project with plan_folds' own function, over the existing keys plus the no_function_set
    # drops (the latter only to report them); write to a temp file next to the output
    out_dir = os.path.dirname(os.path.abspath(out))
    os.makedirs(out_dir, exist_ok=True)
    tmp = os.path.join(out_dir, f".refresh_panel_ec.{os.getpid()}.json")
    try:
        got = PF._cache_panel_ec(a.functions, set(keys) | nofs, tmp, force=True)
    finally:
        if os.path.exists(tmp):
            os.remove(tmp)
    newly_covered = sorted(nofs & set(got))
    new = {g: v for g, v in got.items() if g in old}
    missing = sorted(set(keys) - set(new))
    if missing:
        raise SystemExit(f"[abort] {len(missing)} cached panel genomes have no function set in "
                         f"{a.functions}: {missing[:5]}")
    assert set(new) == set(keys)

    changed = diff_caches(old, new, eligible, universe)
    rep = {
        "what": "refresh of _cache_panel_ec.json for its existing key set (Path A)",
        "tool": os.path.abspath(__file__), "when": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "functions": fn_stamp, "data_dir": os.path.abspath(a.data_dir), "target": target, "out": out,
        "old_cache": {"sha256": sha256_file(target), "n": len(old)},
        "n_keys": len(new), "key_set_unchanged": True,
        "n_changed": len(changed), "n_changed_eligible": sum(1 for c in changed.values() if c["eligible"]),
        "n_lost_any": sum(1 for c in changed.values() if c["lost"]),
        "gained_range_eligible": ([min(c["gained"] for c in changed.values() if c["eligible"]),
                                   max(c["gained"] for c in changed.values() if c["eligible"])]
                                  if any(c["eligible"] for c in changed.values()) else None),
        "new_ec_outside_universe_eligible": sorted({e for c in changed.values() if c["eligible"]
                                                    for e in c["new_outside_universe"]}),
        "no_function_set_now_covered": newly_covered,
        "changed": changed,
    }
    if a.expect_equal:
        ref = PF.load_json(a.expect_equal)
        rep["expect_equal"] = {"path": os.path.abspath(a.expect_equal), "equal_content": new == ref,
                               "equal_key_order": list(new) == list(ref)}
    if not a.dry_run:
        PF.write_json(out, new)
        rep["written"] = {"path": out, "sha256": sha256_file(out)}
        if a.expect_equal:
            rep["expect_equal"]["equal_bytes"] = sha256_file(out) == sha256_file(a.expect_equal)
    g1 = guard_stamp(a.data_dir)
    rep["guarded_before"], rep["guarded_after"] = g0, g1
    rep["guarded_unchanged"] = g0 == g1
    rep["seconds"] = round(time.time() - t0, 1)
    rpath = a.report or ((out + ".refresh.json") if not a.dry_run else None)
    if rpath:
        with open(rpath, "w") as fh:
            json.dump(rep, fh, indent=1)
    ok = rep["guarded_unchanged"] and (not a.expect_equal or rep["expect_equal"]["equal_content"])
    print(f"[refresh_panel_ec] {len(new)} genomes; {len(changed)} changed "
          f"({rep['n_changed_eligible']} eligible; {rep['n_lost_any']} lost any EC); "
          f"{len(rep['new_ec_outside_universe_eligible'])} new ECs outside the universe; "
          f"{len(newly_covered)} no_function_set drops now covered; guarded files unchanged: "
          f"{rep['guarded_unchanged']}"
          + (f"; equal to {a.expect_equal}: {rep['expect_equal']}" if a.expect_equal else "")
          + f" ({rep['seconds']}s)", flush=True)
    if not rep["guarded_unchanged"]:
        print("[refresh_panel_ec] ERROR: panel_eligible.json / panel_ec_universe.json changed", file=sys.stderr)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
