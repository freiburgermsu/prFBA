#!/usr/bin/env python
"""Derive Stage O's weights file `stage_o_weights_E3.json` from Stage E's `stage_o_weights.json`.

    python stage_o_weights_e3.py [--src data/pfunc_estimators/stage_o_weights.json]
                                 [--scores data/pfunc_estimators/scores.json]
                                 [--out data/pfunc_estimators/stage_o_weights_E3.json]
                                 [--generated-e3 ISO8601] [--expect-equal PATH]

Stage E (`pfunc_estimators.py`) writes `stage_o_weights.json` for its winner. On 2026-11-03 the
winner was E4(E3): the E3 configuration plus an all-fold isotonic map in `meta.isotonic`.
Stage O carries E3 itself, not E4(E3) (PFUNC_ESTIMATORS_RESULTS.md, "E4 map in Stage O"): the
E4-vs-E3 contrasts do not separate the two, and the map is in-sample for the 360 Stage O
amplicons, is fitted on per-EC values (at role level it extrapolates) and caps p below 1. That
file was derived by hand at 08:50 on 2026-11-03. This script is its producer, so the derivation
can be re-run after Stage E is re-run (2026-11 re-score). The transformation is:

* `amplicons`: copied unchanged (same keys, same order, same entries);
* `meta.estimator`: the winner's base family ("E4" -> "E3");
* `meta.isotonic`: null;
* `meta.derived_from`: {file: the source's basename, sha256: of the source's bytes, change: one
  sentence whose numbers are read from `scores.json`: the E4(E3) - E3 pairwise contrasts in
  universe Brier, union-support Brier and macro F1 at q >= 0.5 (d100, stratum "all"; estimate
  [lo, hi], 4 dp), the amplicon count and the map's largest output value};
* `meta.generated_e3`: local time, ISO 8601 with seconds (`--generated-e3` overrides it, for the
  control);
* every other `meta` key is kept, in order. Serialisation: `json.dumps` with default separators,
  no indentation and no trailing newline (the form of the 2026-11-03 file).

It refuses unless the source and `scores.json` belong to the same Stage E run and describe an
E3-based winner: the winner's base family is E3, the source's estimator label is the one
Stage E gives that winner ("E4" for E4(E3), "E3" for E3), the source's `script_sha256` equals
`scores.json`'s, its K / tau / dedup / alpha equal the carried-forward variant, and its isotonic
map equals the winner's `isotonic_all_folds`. If Stage E ever selected E3 itself, there is no map
to remove and the file is written with an empty-change note.

Control: `--src`/`--scores` = the committed 2026-11-03 files, `--generated-e3` = the committed
file's `meta.generated_e3`, `--expect-equal` = the committed `stage_o_weights_E3.json`; the bytes
must be equal (exit 1 otherwise).

Interpreter: /home/freiburger/Documents/py_venv/bin/python (the only interpreter).
"""
import argparse
import datetime
import hashlib
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)
PFE = os.path.join(SR, "data", "pfunc_estimators")
METRICS = (("universe_brier", "universe Brier"), ("union_brier", "union Brier"),
           ("macro_F1_q_ge_0.5", "macro F1"))
COST_FLOOR_P = 0.99        # lin cost max(0.01, 1 - p) reaches its floor at p >= 0.99


def sha256_bytes(b):
    return hashlib.sha256(b).hexdigest()


def fmt_ci(c):
    return f"{c['est']:.4f} [{c['lo']:.4f}, {c['hi']:.4f}]"


def covers0(c):
    return c["lo"] <= 0.0 <= c["hi"]


def parse_variant(s):
    """'E3|K=100|tau=0.02|genome|family|alpha=2' -> dict (the fields the weights meta carries)."""
    parts = s.split("|")
    out = {"family": parts[0]}
    for p in parts[1:]:
        if p.startswith("K="):
            out["K"] = int(p[2:])
        elif p.startswith("tau="):
            out["tau"] = float(p[4:])
        elif p.startswith("alpha="):
            out["alpha"] = float(p[6:])
        elif p in ("genome", "species"):
            out["dedup"] = p
    return out


def check_pair(src, scores):
    """Refuse unless src is Stage E's weights file for an E3-based winner of the same run."""
    win = scores.get("winner") or {}
    name, base = win.get("name"), win.get("base_family")
    m = src.get("meta") or {}
    bad = []
    if base != "E3":
        bad.append(f"Stage E's winner is {name!r} (base family {base!r}), not E3-based")
    want_label = "E4" if str(name).startswith("E4") else base
    if m.get("estimator") != want_label:
        bad.append(f"source meta.estimator {m.get('estimator')!r} != {want_label!r} for winner {name!r}")
    if m.get("script_sha256") != (scores.get("meta") or {}).get("script_sha256"):
        bad.append("source and scores.json come from different pfunc_estimators.py versions")
    v = parse_variant(win.get("carried_forward_variant", ""))
    for k in ("K", "tau", "dedup", "alpha"):
        if k in v and m.get(k) != v[k]:
            bad.append(f"source meta.{k} {m.get(k)!r} != carried-forward variant's {v[k]!r}")
    iso_win = win.get("isotonic_all_folds")
    if m.get("isotonic") != iso_win:
        bad.append("source meta.isotonic differs from scores.json winner.isotonic_all_folds")
    if not isinstance(src.get("amplicons"), dict) or not src["amplicons"]:
        bad.append("source has no amplicons")
    if bad:
        raise SystemExit("[stage_o_weights_e3] refusing: " + "; ".join(bad))
    return name, base


def change_text(src, scores, name, base):
    iso = src["meta"]["isotonic"]
    if iso is None:
        return (f"no isotonic map in the source: Stage E selected {name} itself; entries and "
                f"meta copied unchanged")
    pc = scores["composites"]["d100"]["all"]["contrasts_pairwise"][f"{name} - {base}"]
    u, n, f = (pc[k] for k, _ in METRICS)
    if covers0(u) and covers0(n):
        brier = f"universe Brier {fmt_ci(u)} and union Brier {fmt_ci(n)} cover 0"
    else:
        brier = (f"universe Brier {fmt_ci(u)} ({'covers' if covers0(u) else 'excludes'} 0) and "
                 f"union Brier {fmt_ci(n)} ({'covers' if covers0(n) else 'excludes'} 0)")
    cap = max(iso["y"])
    tail = (f"caps p at {cap:.3f} so no reaction reaches the cost floor" if cap < COST_FLOOR_P
            else f"reaches p = {cap:.3f}")
    return (f"isotonic map removed: Stage O carries {base}, not {name}. E4 vs {base} {brier}, "
            f"macro F1 {fmt_ci(f)}; the map is fitted on all five folds (in-sample for these "
            f"{len(src['amplicons'])} amplicons), on per-EC values (applied at role level it "
            f"extrapolates), and {tail}.")


def derive(src_bytes, src_name, scores, generated_e3):
    src = json.loads(src_bytes)
    name, base = check_pair(src, scores)
    meta = dict(src["meta"])                       # key order kept; reassignments stay in place
    meta["estimator"] = base
    meta["isotonic"] = None
    meta["derived_from"] = {"file": src_name, "sha256": sha256_bytes(src_bytes),
                            "change": change_text(src, scores, name, base)}
    meta["generated_e3"] = generated_e3
    out = {"meta": meta, "amplicons": src["amplicons"]}
    for k in src:
        if k not in out:
            out[k] = src[k]
    return out


def serialise(obj):
    return json.dumps(obj).encode()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--src", default=os.path.join(PFE, "stage_o_weights.json"))
    ap.add_argument("--scores", default=os.path.join(PFE, "scores.json"))
    ap.add_argument("--out", default=os.path.join(PFE, "stage_o_weights_E3.json"))
    ap.add_argument("--generated-e3", default=None,
                    help="timestamp for meta.generated_e3 (default: now, local time)")
    ap.add_argument("--expect-equal", default=None,
                    help="compare the output with this file; exit 1 unless the bytes are equal")
    args = ap.parse_args(argv)
    with open(args.src, "rb") as fh:
        src_bytes = fh.read()
    with open(args.scores) as fh:
        scores = json.load(fh)
    ts = args.generated_e3 or datetime.datetime.now().astimezone().isoformat(timespec="seconds")
    obj = derive(src_bytes, os.path.basename(args.src), scores, ts)
    data = serialise(obj)
    tmp = args.out + ".tmp"
    with open(tmp, "wb") as fh:
        fh.write(data)
    os.replace(tmp, args.out)
    src = json.loads(src_bytes)
    rep = {"out": args.out, "sha256": sha256_bytes(data), "bytes": len(data),
           "src": args.src, "src_sha256": obj["meta"]["derived_from"]["sha256"],
           "estimator": obj["meta"]["estimator"], "n_amplicons": len(obj["amplicons"]),
           "amplicons_identical": obj["amplicons"] == src["amplicons"],
           "change": obj["meta"]["derived_from"]["change"], "generated_e3": ts}
    rc = 0
    if args.expect_equal:
        with open(args.expect_equal, "rb") as fh:
            ref = fh.read()
        rep["expect_equal"] = {"path": args.expect_equal, "sha256": sha256_bytes(ref),
                               "equal_bytes": ref == data,
                               "equal_content": json.loads(ref) == obj}
        rc = 0 if ref == data else 1
    print(json.dumps(rep, indent=1))
    return rc


if __name__ == "__main__":
    sys.exit(main())
