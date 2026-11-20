#!/usr/bin/env python
"""Record-level comparison of two select_references selection.json files (re-score R1.2 C-L*,
R1.3/R1.4 legacy checks). Read-only, streaming (ijson), memory ~ a few GB.

    cmp_selection.py --old OLD --new NEW --out report.json
                     [--repaired $W/expected/repaired_ids.json] [--d20]
                     [--detail-max 5000] [--keys-out changed_keys.json]

Pass 1 streams each file once and keeps, per record key, sha256 of
orjson.dumps(record, OPT_SORT_KEYS) (serialisation-independent), the include-mode n_selected,
and whether any genome the record names (selected / provenance `genome_id`, `equivalent_genomes`)
is in the repaired set ("touched"). `_meta` is compared separately with `generated_from` masked;
knob and summary differences are listed.

Pass 2 (only when records differ) re-streams both files for up to --detail-max differing keys and
classifies each difference by its leaf-path signature (list indices generalised to []),
plus: selected genome list changed (ordered / as a set), n_selected changed, and the provenance
stages that changed (`no_gene_set` -> ... etc.).

Verdicts: `records_identical` (same key set, every digest equal); `changed_subset_of_touched`
(every differing record is touched in old or new).
"""
import argparse
import collections
import hashlib
import json
import os
import sys
import time

import ijson
import orjson

W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
MASK_META = ("generated_from",)


def genome_ids(obj, out):
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k == "genome_id" and v is not None:
                out.add(str(v))
            elif k == "equivalent_genomes" and isinstance(v, list):
                for g in v:
                    out.add(str(g.get("genome_id") if isinstance(g, dict) else g))
            else:
                genome_ids(v, out)
    elif isinstance(obj, list):
        for v in obj:
            genome_ids(v, out)


def pass1(path, repaired):
    digests, nsel, touched, meta = {}, {}, set(), None
    with open(path, "rb") as fh:
        for k, rec in ijson.kvitems(fh, "", use_float=True):
            if k == "_meta":
                meta = rec
                continue
            digests[k] = hashlib.sha256(orjson.dumps(rec, option=orjson.OPT_SORT_KEYS)).digest()
            if isinstance(rec, dict):
                nsel[k] = rec.get("n_selected")
                ids = set()
                genome_ids(rec, ids)
                if ids & repaired:
                    touched.add(k)
    return digests, nsel, touched, meta


def pass2(path, keys):
    out = {}
    with open(path, "rb") as fh:
        for k, rec in ijson.kvitems(fh, "", use_float=True):
            if k in keys:
                out[k] = rec
                if len(out) == len(keys):
                    break
    return out


def leaf_diffs(a, b, path, out):
    if isinstance(a, dict) and isinstance(b, dict):
        for k in set(a) | set(b):
            p = f"{path}.{k}"
            if k not in a:
                out.add(p + " (added)")
            elif k not in b:
                out.add(p + " (removed)")
            else:
                leaf_diffs(a[k], b[k], p, out)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.add(path + " (length)")
        for x, y in zip(a, b):
            leaf_diffs(x, y, path + "[]", out)
    elif a != b:
        out.add(path)


def sel_ids(rec):
    return [s.get("genome_id") for s in (rec.get("selected") or []) if isinstance(s, dict)]


def prov_map(rec):
    """Per-candidate (disposition, stage): the audit list is `decisions` (provenance is a dict of
    summary fields in these selections); older layouts kept it in a `provenance` list."""
    lst = rec.get("decisions")
    if not isinstance(lst, list):
        lst = rec.get("provenance") if isinstance(rec.get("provenance"), list) else []
    return {(p.get("genome_id"), p.get("rank")): (p.get("disposition"), p.get("stage"))
            for p in lst if isinstance(p, dict)}


def meta_diff(mo, mn):
    mo = {k: v for k, v in (mo or {}).items() if k not in MASK_META}
    mn = {k: v for k, v in (mn or {}).items() if k not in MASK_META}
    d = {}
    for k in sorted(set(mo) | set(mn)):
        if mo.get(k) != mn.get(k):
            if isinstance(mo.get(k), dict) and isinstance(mn.get(k), dict):
                sub = {kk: [mo[k].get(kk), mn[k].get(kk)] for kk in sorted(set(mo[k]) | set(mn[k]))
                       if mo[k].get(kk) != mn[k].get(kk)}
                d[k] = sub
            else:
                d[k] = [mo.get(k), mn.get(k)]
    return d


def mean(xs):
    xs = [x for x in xs if isinstance(x, (int, float))]
    return round(sum(xs) / len(xs), 6) if xs else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--repaired", default=os.path.join(W, "expected", "repaired_ids.json"))
    ap.add_argument("--d20", action="store_true", help="also count the 20 d20 poisoned keys as repaired")
    ap.add_argument("--detail-max", type=int, default=5000)
    ap.add_argument("--keys-out")
    a = ap.parse_args()
    t0 = time.time()
    rj = json.load(open(a.repaired))
    repaired = set(rj["repaired_2026_11_13"]) | (set(rj.get("poisoned_d20", [])) if a.d20 else set())
    do, no_, to, mo = pass1(a.old, repaired)
    t1 = time.time()
    dn, nn, tn, mn = pass1(a.new, repaired)
    t2 = time.time()
    ko, kn = set(do), set(dn)
    common = ko & kn
    changed = sorted(k for k in common if do[k] != dn[k])
    touched = to | tn
    untouched_changed = [k for k in changed if k not in touched]
    rep = {
        "old": os.path.abspath(a.old), "new": os.path.abspath(a.new),
        "n_repaired_set": len(repaired), "d20": a.d20,
        "n_old": len(do), "n_new": len(dn), "n_only_old": len(ko - kn), "n_only_new": len(kn - ko),
        "only_old_examples": sorted(ko - kn)[:20], "only_new_examples": sorted(kn - ko)[:20],
        "n_common": len(common), "n_identical": len(common) - len(changed), "n_changed": len(changed),
        "n_touched_old": len(to), "n_touched_new": len(tn), "n_touched_union": len(touched),
        "n_changed_touched": len(changed) - len(untouched_changed),
        "n_changed_untouched": len(untouched_changed),
        "changed_untouched_examples": untouched_changed[:30],
        "meta_diff_masked": meta_diff(mo, mn),
        "meta_summary_old": (mo or {}).get("summary"), "meta_summary_new": (mn or {}).get("summary"),
        "mean_n_selected_old": mean(no_.values()), "mean_n_selected_new": mean(nn.values()),
        "mean_n_selected_nonzero_old": mean([x for x in no_.values() if x]),
        "mean_n_selected_nonzero_new": mean([x for x in nn.values() if x]),
        "n_selected_changed": sum(1 for k in common if no_[k] != nn[k]),
        "seconds_pass1": [round(t1 - t0, 1), round(t2 - t1, 1)],
    }
    rep["records_identical"] = (rep["n_changed"] == 0 and rep["n_only_old"] == 0 and rep["n_only_new"] == 0)
    rep["changed_subset_of_touched"] = (rep["n_changed_untouched"] == 0 and rep["n_only_old"] == 0
                                       and rep["n_only_new"] == 0)
    if changed:
        # detail: all untouched changes first (they need explaining), then the rest up to the cap
        pick = untouched_changed + [k for k in changed if k in touched]
        pick = set(pick[: a.detail_max])
        ro, rn = pass2(a.old, pick), pass2(a.new, pick)
        sigs, paths = collections.Counter(), collections.Counter()
        stage_moves = collections.Counter()
        n_ids_ordered = n_ids_set = n_nsel = 0
        examples = []
        for k in sorted(pick):
            o, n = ro.get(k), rn.get(k)
            if o is None or n is None:
                continue
            ds = set()
            leaf_diffs(o, n, "", ds)
            for p in ds:
                paths[p] += 1
            sigs[" | ".join(sorted(ds))[:400]] += 1
            io, inn = sel_ids(o), sel_ids(n)
            n_ids_ordered += io != inn
            n_ids_set += set(io) != set(inn)
            n_nsel += o.get("n_selected") != n.get("n_selected")
            po, pn = prov_map(o), prov_map(n)
            for g in set(po) | set(pn):
                if po.get(g) != pn.get(g):
                    stage_moves[f"{po.get(g, ('absent', None))[1] or po.get(g, ('absent',))[0]} -> "
                                f"{pn.get(g, ('absent', None))[1] or pn.get(g, ('absent',))[0]}"] += 1
            if len(examples) < 15:
                examples.append({"key": k, "touched": k in touched, "selected_old": io, "selected_new": inn,
                                 "diff_paths": sorted(ds)[:20]})
        rep["detail"] = {
            "n_detailed": len(pick), "selected_ids_changed_ordered": n_ids_ordered,
            "selected_ids_changed_as_set": n_ids_set, "n_selected_changed": n_nsel,
            "top_paths": paths.most_common(40), "top_signatures": sigs.most_common(15),
            "provenance_stage_moves": stage_moves.most_common(30), "examples": examples,
        }
    if a.keys_out:
        with open(a.keys_out, "w") as fh:
            json.dump({"changed": changed, "touched_old": sorted(to), "touched_new": sorted(tn)}, fh)
    rep["seconds"] = round(time.time() - t0, 1)
    with open(a.out, "w") as fh:
        json.dump(rep, fh, indent=1, default=str)
    brief = {k: rep[k] for k in ("n_old", "n_new", "n_only_old", "n_only_new", "n_changed", "n_touched_union",
                                 "n_changed_untouched", "n_selected_changed", "mean_n_selected_old",
                                 "mean_n_selected_new", "records_identical", "changed_subset_of_touched",
                                 "seconds")}
    print(json.dumps(brief))
    print("meta_diff_masked:", json.dumps(rep["meta_diff_masked"])[:3000])
    print("RECORDS IDENTICAL" if rep["records_identical"] else
          ("CHANGED SUBSET OF TOUCHED" if rep["changed_subset_of_touched"] else "UNATTRIBUTED CHANGES"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
