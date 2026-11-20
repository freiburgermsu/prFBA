#!/usr/bin/env python
"""W8 for the row-level tables Run 2 regenerates that attrib_diff.py has no profile for.

    attrib_rows.py --out report.json      (reads fixed paths; old = $W/pre_repair or git start tag)

Tables and the reasons that may explain a changed (or one-sided) row:
  selfrec_d100_legacy  region_validation_d100_legacy/data/self_recovery_failures.csv (key amplicon):
                       the amplicon's legacy selection record changed in R1.3 (every such record is
                       touched by a repaired genome: control/r13/sel_region_validation_d100_legacy.json)
  ceiling_d20          annotation_ceiling.json per_genome (key genome_id): genome repaired (d20 set)
  depth_probe_d20      depth_probe_results.json rows (key amplicon): source or a rep (old or new)
                       repaired, or a rep whose GFS entry changed through newly fetched products (D10)
  afgate               af_gate_counterfactual.json mags (key mag): source repaired, or a skani hit
                       reference repaired (the only genomes a re-selection can pick), MAG panel indices
  magarm_cells_d20     P2C data/mag_arm/cells.csv.gz (optional, when present)
  magarm_sel_d20       P2C data/mag_arm/prfba_mag_selections.csv.gz (optional, when present)
"""
import argparse
import csv
import gzip
import io
import json
import os
import subprocess
import sys
from collections import Counter

W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
R = "/home/freiburger/Documents/prFBA"
TAG = "rescore-2026-11-start"
P2C = "region_validation/picrust2_comparison/data"


def git_old(path):
    return subprocess.run(["git", "-C", R, "show", f"{TAG}:{path}"], capture_output=True,
                          text=True, check=True).stdout


def repaired(d20):
    r = json.load(open(f"{W}/expected/repaired_ids.json"))
    return set(r["repaired_2026_11_13"]) | (set(r["poisoned_d20"]) if d20 else set())


def compare(old, new, explain, name, max_ex=15):
    """old/new: {key: row}; explain(key, o, n) -> list of reasons."""
    keys = set(old) | set(new)
    rep = {"table": name, "n_old": len(old), "n_new": len(new),
           "n_only_old": len(set(old) - set(new)), "n_only_new": len(set(new) - set(old))}
    first, anyc, unexpl = Counter(), Counter(), []
    n_changed = 0
    for k in sorted(keys):
        o, n = old.get(k), new.get(k)
        if o == n:
            continue
        n_changed += 1
        rs = explain(k, o, n) or ["unattributed"]
        first[rs[0]] += 1
        for x in rs:
            anyc[x] += 1
        if rs[0] == "unattributed" and len(unexpl) < max_ex:
            unexpl.append({"key": k, "old": o, "new": n})
    rep.update(n_changed_or_one_sided=n_changed, by_reason_first=dict(first),
               by_reason_any=dict(anyc), n_unattributed=first.get("unattributed", 0),
               unattributed_examples=unexpl, ok=first.get("unattributed", 0) == 0)
    return rep


def t_selfrec():
    path = "region_validation_d100_legacy/data/self_recovery_failures.csv"
    rd = lambda txt: {r["amplicon"]: r for r in csv.DictReader(io.StringIO(txt))}  # noqa: E731
    old, new = rd(git_old(path)), rd(open(f"{R}/{path}").read())
    changed = set(json.load(open(f"{W}/control/r13/sel_region_validation_d100_legacy.keys.json"))["changed"])
    return compare(old, new, lambda k, o, n: ["selection_changed_R1.3_touched"] if k in changed else [],
                   "selfrec_d100_legacy")


def t_ceiling():
    rp = repaired(True)
    rows = lambda d: {g["genome_id"]: g for g in d["per_genome"]}  # noqa: E731
    old = rows(json.load(open(f"{W}/pre_repair/{P2C}/annotation_ceiling.json")))
    new = rows(json.load(open(f"{R}/{P2C}/annotation_ceiling.json")))
    return compare(old, new, lambda k, o, n: ["source_repaired"] if k in rp else [], "ceiling_d20")


def t_depth_probe():
    rp = repaired(True)
    np_ = set(json.load(open(f"{W}/control/r2/d20/gfs_changed_newproducts.json")))
    rows = lambda d: {r["amplicon"]: r for r in d["rows"]}  # noqa: E731
    old = rows(json.load(open(f"{W}/pre_repair/{P2C}/depth_probe/depth_probe_results.json")))
    new = rows(json.load(open(f"{R}/{P2C}/depth_probe/depth_probe_results.json")))

    def why(k, o, n):
        rs = []
        if k.split("__")[0] in rp:
            rs.append("source_repaired")
        reps = set((o or {}).get("reps", [])) | set((n or {}).get("reps", []))
        if reps & rp:
            rs.append("rep_repaired")
        if reps & np_:
            rs.append("rep_newproducts")
        return rs
    return compare(old, new, why, "depth_probe_d20")


def t_afgate():
    rp = repaired(False)
    old = json.loads(git_old("mag_validation/out/af_gate_counterfactual.json"))["mags"]
    new = json.load(open(f"{R}/mag_validation/out/af_gate_counterfactual.json"))["mags"]
    src, refs = {}, {}
    for i in range(20):
        for r in json.load(open(f"{R}/mag_validation/manifests/b{i:03d}.json"))["rows"]:
            src[r["mag"]] = r["genome_id"]
        with open(f"{R}/mag_validation/hits/b{i:03d}.tsv") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                q = os.path.basename(row["Query_file"])[:-4]
                if q in old or q in new:
                    refs.setdefault(q, set()).add(os.path.basename(row["Ref_file"])[:-4])

    def why(k, o, n):
        rs = []
        if src.get(k) in rp:
            rs.append("source_repaired")
        if refs.get(k, set()) & rp:
            rs.append("hit_ref_repaired")
        return rs
    return compare(old, new, why, "afgate")


def t_magarm():
    out = []
    old_p = f"{W}/pre_repair/{P2C}/mag_arm/prfba_mag_selections.csv.gz"
    new_p = f"{R}/{P2C}/mag_arm/prfba_mag_selections.csv.gz"
    if not (os.path.exists(old_p) and os.path.exists(new_p)):
        return out
    rp = repaired(True)
    rd = lambda p: {(r["mag"], r["regime"]): r for r in csv.DictReader(gzip.open(p, "rt"))}  # noqa: E731
    old, new = rd(old_p), rd(new_p)
    refs = {}
    import glob
    for hp in sorted(glob.glob(f"{R}/mag_validation/hits/*.tsv") + glob.glob(f"{R}/mag_validation/lowc_hits/*.tsv")):
        with open(hp) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                refs.setdefault(os.path.basename(row["Query_file"])[:-4], set()).add(
                    os.path.basename(row["Ref_file"])[:-4])

    def why_sel(k, o, n):
        rs = []
        if (n or o)["genome_id"] in rp:
            rs.append("source_repaired")
        picks = set(((o or {}).get("picks") or "").split(";")) | set(((n or {}).get("picks") or "").split(";"))
        if picks & rp:
            rs.append("pick_repaired")
        if refs.get(k[0], set()) & rp:
            rs.append("hit_ref_repaired")
        if (o or {}).get("matches_recorded") != (n or {}).get("matches_recorded") and \
                (o or {}).get("picks") == (n or {}).get("picks"):
            rs.append("recorded_table_rescored_R2.2")
        return rs
    out.append(compare(old, new, why_sel, "magarm_sel_d20"))
    old_c = f"{W}/pre_repair/{P2C}/mag_arm/cells.csv.gz"
    new_c = f"{R}/{P2C}/mag_arm/cells.csv.gz"
    if os.path.exists(old_c) and os.path.exists(new_c):
        hdr = next(csv.reader(gzip.open(new_c, "rt")))
        kcols = [c for c in ("mag", "genome_id", "completeness", "regime", "tier", "method") if c in hdr]

        def rdc(p):
            d = {}
            for r in csv.DictReader(gzip.open(p, "rt")):
                d[tuple(r.get(c) for c in kcols)] = r
            return d
        oc, nc = rdc(old_c), rdc(new_c)
        gcol = "genome_id" if "genome_id" in hdr else None
        np_ = set(json.load(open(f"{W}/control/r2/d20/gfs_changed_newproducts.json")))
        sel_changed = {k[0] for k in set(old) | set(new) if old.get(k) != new.get(k)}

        picks = {}
        for (mag, regime), row in new.items():
            picks[(mag, "prfba_" + regime)] = set(filter(None, (row.get("picks") or "").split(";")))

        def why_cell(k, o, n):
            r = n or o
            rs = []
            g = r.get(gcol) if gcol else None
            if g in rp:
                rs.append("source_repaired")
            pk = picks.get((r.get("mag"), r.get("method")), set())
            if pk & rp:
                rs.append("pick_repaired")
            if pk & np_:
                rs.append("pick_newproducts")
            if r.get("mag") in sel_changed:
                rs.append("selection_changed")
            if g in np_:
                rs.append("source_newproducts")
            if r.get("method") == "base_constant":
                # mag_score's calibration floor: the K most prevalent ECs over ALL panel truths
                # (K = median truth size); the 140 repaired panel genomes move it for every MAG
                rs.append("panel_prevalence_baseline")
            return rs
        rep = compare(oc, nc, why_cell, "magarm_cells_d20")
        rep["key_columns"] = kcols
        out.append(rep)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", nargs="*")
    a = ap.parse_args(argv)
    fns = {"selfrec": t_selfrec, "ceiling": t_ceiling, "depth_probe": t_depth_probe,
           "afgate": t_afgate, "magarm": t_magarm}
    reps = []
    for name, fn in fns.items():
        if a.only and name not in a.only:
            continue
        r = fn()
        reps.extend(r if isinstance(r, list) else [r])
    json.dump(reps, open(a.out, "w"), indent=1, default=str)
    for r in reps:
        print(json.dumps({k: r[k] for k in ("table", "n_old", "n_new", "n_only_old", "n_only_new",
                                            "n_changed_or_one_sided", "by_reason_first", "n_unattributed", "ok")}))
    return 0 if all(r["ok"] for r in reps) else 1


if __name__ == "__main__":
    sys.exit(main())
