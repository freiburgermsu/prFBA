#!/usr/bin/env python
"""Compare two shared-reference prFBA predictions.json files (run_prfba_fold.py output).

    cmp_prfba_pred.py OLD NEW [--out report.json] [--changed-genomes ids.json]

Content comparison: `pred` arm by arm, cell by cell (every field), and `meta` minus the run-specific
fields (generated, git_commit, wall_seconds, align_seconds). With --changed-genomes (a JSON list, or
a dict with key `changed_eligible`), every changed cell is attributed: reps (old or new) containing a
changed genome -> `rep_changed_genome`; else `unattributed`. Exit 0 iff identical (no --changed-genomes)
or iff 0 unattributed (with it).
"""
import argparse, json, sys
import orjson

RUN_META = {"generated", "git_commit", "wall_seconds", "align_seconds"}


def load(p):
    with open(p, "rb") as fh:
        return orjson.loads(fh.read())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("old"); ap.add_argument("new")
    ap.add_argument("--out"); ap.add_argument("--changed-genomes")
    a = ap.parse_args()
    o, n = load(a.old), load(a.new)
    ch = None
    if a.changed_genomes:
        x = json.load(open(a.changed_genomes))
        ch = set(map(str, x["changed_eligible"] if isinstance(x, dict) else x))
    mo = {k: v for k, v in o["meta"].items() if k not in RUN_META}
    mn = {k: v for k, v in n["meta"].items() if k not in RUN_META}
    meta_diff = sorted(k for k in set(mo) | set(mn) if mo.get(k) != mn.get(k))
    rep = {"old": a.old, "new": a.new, "meta_equal_minus_run_fields": not meta_diff,
           "meta_diff_keys": meta_diff,
           "meta_diff": {k: {"old": mo.get(k), "new": mn.get(k)} for k in meta_diff if k != "assertions"},
           "arms_old": sorted(o["pred"]), "arms_new": sorted(n["pred"]), "arms": {}}
    if "assertions" in meta_diff:
        ao, an = mo.get("assertions") or {}, mn.get("assertions") or {}
        rep["meta_diff"]["assertions"] = {k: {"old": ao.get(k), "new": an.get(k)}
                                          for k in set(ao) | set(an) if ao.get(k) != an.get(k)}
    unattr = 0; total_changed = 0
    for arm in sorted(set(o["pred"]) | set(n["pred"])):
        po, pn = o["pred"].get(arm, {}), n["pred"].get(arm, {})
        keys = set(po) | set(pn)
        one_sided = sorted(k for k in keys if (k in po) != (k in pn))
        changed = [k for k in keys if k in po and k in pn and po[k] != pn[k]]
        fields = {}
        reasons = {"rep_changed_genome": 0, "unattributed": 0}
        unattr_ex = []
        for k in changed:
            for f in set(po[k]) | set(pn[k]):
                if po[k].get(f) != pn[k].get(f):
                    fields[f] = fields.get(f, 0) + 1
            if ch is not None:
                reps = set(map(str, po[k].get("reps") or [])) | set(map(str, pn[k].get("reps") or []))
                if reps & ch:
                    reasons["rep_changed_genome"] += 1
                else:
                    reasons["unattributed"] += 1
                    if len(unattr_ex) < 10:
                        unattr_ex.append({"md5": k, "old": po[k], "new": pn[k]})
        ec_gain = sum(len(set(pn[k]["ec"]) - set(po[k]["ec"])) for k in changed)
        ec_loss = sum(len(set(po[k]["ec"]) - set(pn[k]["ec"])) for k in changed)
        rep["arms"][arm] = {"n_old": len(po), "n_new": len(pn), "one_sided": len(one_sided),
                            "n_changed": len(changed), "changed_fields": fields,
                            "ec_gained_total": ec_gain, "ec_lost_total": ec_loss,
                            **({"attribution": reasons, "unattributed_examples": unattr_ex} if ch is not None else {})}
        total_changed += len(changed) + len(one_sided)
        unattr += reasons["unattributed"] + len(one_sided)
    rep["identical_pred"] = total_changed == 0
    rep["n_changed_total"] = total_changed
    if ch is not None:
        rep["n_unattributed_total"] = unattr
    s = json.dumps(rep, indent=1)
    if a.out:
        open(a.out, "w").write(s)
    print(json.dumps({k: rep[k] for k in ("identical_pred", "n_changed_total", "meta_equal_minus_run_fields", "meta_diff_keys")}
                     | ({"n_unattributed_total": unattr} if ch is not None else {})))
    ok = (unattr == 0) if ch is not None else (total_changed == 0 and not meta_diff)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
