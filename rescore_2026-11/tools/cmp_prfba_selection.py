#!/usr/bin/env python
"""Content comparison of two prfba_selection.json files (re-score R2.4/R2.5 control and diff).

    cmp_prfba_selection.py --old OLD --new NEW --out report.json
        [--hits HITS.json] [--repaired $W/expected/repaired_ids.json] [--d20]
        [--cache-keys $W/stamps/R1.6.cache_keys.txt] [--max-examples 30]

Records: every amplicon key's record compared as a whole (dict equality; reps lists are ordered
as written). Differences are classified per regime (include / exclude_0.987): reps (as list, as
set), tier, n_selected, best_identity, n_removed_self; and best_nonself_identity (never expected to
move: it reads hits only). `_meta` is compared without `git_commit` and `inputs`.

With --hits, a record is "touched" when any genome of its top20 hits (the candidates of both
regimes) or of its old/new reps is in the repaired set: the only records a data-only change may
move. With --cache-keys (the key list of the frozen cache, in file = insertion order), each genome
named by a differing record is reported with its insertion index, so a difference caused by a
genome that entered the cache after the committed run (cache history) can be told from code.
Exit 0 iff records are identical (control) -- the verdict for a real old->new diff is
`changed_subset_of_touched`.
"""
import argparse
import json
import sys
from collections import Counter

import orjson

W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
REGIMES = ("include", "exclude_0.987")
META_MASK = ("git_commit", "inputs")


def load(path):
    d = orjson.loads(open(path, "rb").read())
    meta = d.pop("_meta", {})
    return d, meta


def repaired_set(path, d20):
    r = json.load(open(path))
    s = set(r["repaired_2026_11_13"])
    if d20:
        s |= set(r["poisoned_d20"])
    return s


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--hits")
    ap.add_argument("--repaired", default=f"{W}/expected/repaired_ids.json")
    ap.add_argument("--d20", action="store_true")
    ap.add_argument("--cache-keys")
    ap.add_argument("--max-examples", type=int, default=30)
    a = ap.parse_args(argv)

    old, mo = load(a.old)
    new, mn = load(a.new)
    rep = {"old": a.old, "new": a.new, "n_old": len(old), "n_new": len(new)}
    only_old = sorted(set(old) - set(new))
    only_new = sorted(set(new) - set(old))
    rep.update(n_only_old=len(only_old), n_only_new=len(only_new),
               only_old_examples=only_old[:10], only_new_examples=only_new[:10])
    diff = sorted(k for k in set(old) & set(new) if old[k] != new[k])
    rep["n_records_differ"] = len(diff)
    kinds = Counter()
    involved = Counter()
    for k in diff:
        o, n = old[k], new[k]
        if o.get("best_nonself_identity") != n.get("best_nonself_identity"):
            kinds["best_nonself_identity"] += 1
        for rg in REGIMES:
            ro, rn = o.get(rg, {}), n.get(rg, {})
            if ro == rn:
                continue
            kinds[f"{rg}|any"] += 1
            if ro.get("reps") != rn.get("reps"):
                kinds[f"{rg}|reps_list"] += 1
                if set(ro.get("reps", [])) != set(rn.get("reps", [])):
                    kinds[f"{rg}|reps_set"] += 1
                for g in set(ro.get("reps", [])) ^ set(rn.get("reps", [])):
                    involved[g] += 1
            for f in ("tier", "n_selected", "best_identity", "n_removed_self", "no_hits"):
                if ro.get(f) != rn.get(f):
                    kinds[f"{rg}|{f}"] += 1
            if not ro.get("n_selected") and rn.get("n_selected"):
                kinds[f"{rg}|abstain->answer"] += 1
            if ro.get("n_selected") and not rn.get("n_selected"):
                kinds[f"{rg}|answer->abstain"] += 1
    rep["differences_by_kind"] = dict(sorted(kinds.items()))
    rep["examples"] = [{"key": k, "old": old[k], "new": new[k]} for k in diff[:a.max_examples]]
    rep["n_genomes_entering_or_leaving_reps"] = len(involved)

    for name, d in (("old", old), ("new", new)):
        ns = [d[k][rg].get("n_selected", 0) for k in d for rg in REGIMES if rg in d[k]]
        rep[f"{name}_summary"] = {
            rg: {"n_abstain": sum(1 for k in d if not d[k].get(rg, {}).get("n_selected")),
                 "mean_n_selected_answered": round(
                     sum(d[k][rg]["n_selected"] for k in d if d[k].get(rg, {}).get("n_selected"))
                     / max(1, sum(1 for k in d if d[k].get(rg, {}).get("n_selected"))), 6)}
            for rg in REGIMES}
        del ns

    mmo = {k: v for k, v in mo.items() if k not in META_MASK}
    mmn = {k: v for k, v in mn.items() if k not in META_MASK}
    rep["meta_diff"] = {k: {"old": mmo.get(k), "new": mmn.get(k)}
                        for k in sorted(set(mmo) | set(mmn)) if mmo.get(k) != mmn.get(k)}
    rep["meta_old"], rep["meta_new"] = mo, mn

    repaired = repaired_set(a.repaired, a.d20)
    if a.hits:
        hits = orjson.loads(open(a.hits, "rb").read())
        untouched = []
        n_touched_total = 0
        for k in set(old) | set(new):
            hr = hits.get(k) or {}
            g = {str(h.get("genome_id")) for h in (hr.get("top20") or [])}
            for d in (old, new):
                for rg in REGIMES:
                    g.update((d.get(k) or {}).get(rg, {}).get("reps", []))
            t = bool(g & repaired)
            n_touched_total += t
            if k in diff and not t:
                untouched.append(k)
        del hits
        rep["n_touched_records"] = n_touched_total
        rep["n_changed_untouched"] = len(untouched)
        rep["changed_untouched_examples"] = [{"key": k, "old": old[k], "new": new[k]}
                                             for k in sorted(untouched)[:a.max_examples]]
        rep["changed_subset_of_touched"] = not untouched and not only_old and not only_new
    if a.cache_keys and involved:
        pos = {}
        with open(a.cache_keys) as fh:
            for i, line in enumerate(fh):
                g = line.rstrip("\n")
                if g in involved:
                    pos[g] = i
        rep["involved_genomes"] = {g: {"n_records": c, "cache_key_index": pos.get(g),
                                       "repaired": g in repaired}
                                   for g, c in involved.most_common()}
    rep["records_identical"] = not diff and not only_old and not only_new
    rep["meta_equal_masked"] = not rep["meta_diff"]
    json.dump(rep, open(a.out, "w"), indent=1, default=str)
    short = {k: rep[k] for k in ("n_old", "n_new", "n_only_old", "n_only_new", "n_records_differ",
                                 "differences_by_kind", "meta_diff", "records_identical")}
    for k in ("n_touched_records", "n_changed_untouched", "changed_subset_of_touched"):
        if k in rep:
            short[k] = rep[k]
    print(json.dumps(short, indent=1, default=str)[:4000])
    return 0 if rep["records_identical"] else 1


if __name__ == "__main__":
    sys.exit(main())
