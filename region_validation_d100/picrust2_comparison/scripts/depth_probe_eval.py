#!/usr/bin/env python
"""Does deeper retention recover the self-exclusion abstentions it caused?

Takes the sampled amplicons whose entire retained top-20 was deleted by the
self-exclusion filter, re-aligned with ``--k2 500``, and re-runs the SAME
post-processing the benchmark uses -- ``score_genecap.exclude_self_record`` then
``select_references.select_representatives`` -- on the deeper record.

Reports, per region and overall: how many now answer, at what tier, how deep the
first surviving (non-self, < tau) hit sits (i.e. what k2 would have sufficed), the
identity of that survivor, and whether the resulting anchor is taxonomically
correct -- because recovering an answer is only an improvement if the answer is
right.
"""
import argparse
import json
import os
import statistics as st
import sys
from collections import Counter, defaultdict

import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))
import _common as C  # noqa: E402
import score_genecap as SG  # noqa: E402
import score_taxacc as ST  # noqa: E402
import select_references as SR  # noqa: E402

PATHS = C.PATHS
DEPTH_BINS = [(20, "<=20 (shipped)"), (50, "21-50"), (100, "51-100"),
              (200, "101-200"), (500, "201-500")]


def depth_bin(rank):
    if rank is None:
        return "none in 500"
    for edge, label in DEPTH_BINS:
        if rank <= edge:
            return label
    return "none in 500"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    D = os.path.join(HERE, "..", "data", "depth_probe")
    ap.add_argument("--hits", default=os.path.join(D, "hits_k500",
                                                   "asv_top20_alignment_hits.json"))
    ap.add_argument("--keys", default=os.path.join(D, "wiped_sample.keys.json"))
    ap.add_argument("--out", default=os.path.join(D, "depth_probe_results.json"))
    args = ap.parse_args()

    keys = json.load(open(args.keys))
    truth = orjson.loads(open(PATHS.truth_json, "rb").read())
    fanout = orjson.loads(open(PATHS.amplicon_fanout, "rb").read())
    md5_of = {k: m for m, ks in fanout.items() for k in ks}
    hits = orjson.loads(open(args.hits, "rb").read())
    hits = {k: v for k, v in hits.items() if k != "_meta" and isinstance(v, dict)}

    cache = orjson.loads(open(PATHS.pgfam_cache, "rb").read())

    def provider(gid):
        v = cache.get(str(gid))
        return set(v) if v is not None else set()
    provider.cache = cache
    gcf = SR.gene_count_from_provider(provider)

    rows = []
    for key in keys["amplicons"]:
        t = truth[key]
        src, src_tid, region = t["src_genome_id"], t["src_taxon_id"], t["region"]
        hr = hits.get(md5_of[key])
        if hr is None:
            continue
        n_hits = len(hr.get("top20") or [])
        excl_rec, n_removed = SG.exclude_self_record(hr, src, src_tid, SG.TAU)
        survivors = excl_rec.get("top20") or []
        # depth of the first survivor in the ORIGINAL ranking = the k2 that suffices
        first_rank = None
        for i, h in enumerate(hr.get("top20") or [], 1):
            hid = str(h.get("genome_id")) if h.get("genome_id") is not None else None
            htid = h.get("taxon_id")
            ident = h.get("identity") or 0.0
            if hid == src or (src_tid is not None and htid is not None and int(htid) == src_tid) \
                    or ident >= SG.TAU:
                continue
            first_rank = i
            break
        sel = SR.select_representatives(excl_rec, gene_provider=provider, gene_count_fn=gcf)
        answered = (sel.get("n_selected") or 0) > 0
        row = {"amplicon": key, "region": region, "n_hits_retained": n_hits,
               "n_removed": n_removed, "n_survivors": len(survivors),
               "first_survivor_rank": first_rank,
               "first_survivor_identity": (survivors[0].get("identity") if survivors else None),
               "answered": answered, "tier": sel.get("tier"),
               "n_selected": sel.get("n_selected", 0),
               "reps": SG.selected_rep_gids(sel)}
        if answered:
            a_lin, a_gid, a_tid = ST.anchor_lineage(sel)
            truth_lin = t["truth_lineage"]
            for rank in ("Family", "Genus"):
                row[f"{rank.lower()}_correct"] = ST.rank_correct(truth_lin.get(rank),
                                                                 (a_lin or {}).get(rank))
        rows.append(row)

    # --- EC-space quality of the recovered answers -------------------------- #
    # Same projection and shared vocabulary as compare.py, so these numbers sit on
    # the same axis as the prFBA-vs-PICRUSt2 table.
    from annotation_ceiling import load_reference_ec
    maps = orjson.loads(open(os.path.join(HERE, "..", "data",
                                          "function_maps.json"), "rb").read())
    pg_ec = maps["pgfam_ec"]
    _ref_rows, ref_vocab = load_reference_ec()
    shared_ec = {e for v in pg_ec.values() for e in v} & ref_vocab

    def ec_of(gid):
        return {e for p in cache.get(str(gid), ()) for e in pg_ec.get(p, ())} & shared_ec

    ec_stats = {"recall": [], "precision": [], "f1": []}
    for r in rows:
        if not r["answered"]:
            continue
        T = ec_of(truth[r["amplicon"]]["src_genome_id"])
        if not T:
            continue
        Pr = set().union(*[ec_of(g) for g in r["reps"]]) if r["reps"] else set()
        tp = len(T & Pr)
        rec = tp / len(T)
        pre = tp / len(Pr) if Pr else 0.0
        f1 = 2 * rec * pre / (rec + pre) if (rec and pre) else 0.0
        r["ec_recall"], r["ec_precision"], r["ec_f1"] = round(rec, 4), round(pre, 4), round(f1, 4)
        ec_stats["recall"].append(rec)
        ec_stats["precision"].append(pre)
        ec_stats["f1"].append(f1)

    n = len(rows)
    answered = [r for r in rows if r["answered"]]
    out = {
        "_meta": {
            "question": "do the self-exclusion abstentions caused by k2=20 retention "
                        "recover when the same amplicons are re-aligned with k2=500?",
            "n_sampled": n, "n_wiped_total_in_benchmark": keys["n_wiped_total"],
            "hits": args.hits, "tau": SG.TAU,
            "selector": "shipped default (exact-tie union + species dedup + P3)",
        },
        "recovery": {
            "answered": len(answered), "rate": round(len(answered) / n, 4) if n else None,
            "still_empty": n - len(answered),
        },
        "depth_needed": dict(Counter(depth_bin(r["first_survivor_rank"]) for r in rows)),
        "tier": dict(Counter(r["tier"] for r in answered)),
        "first_survivor_identity": {},
        "anchor_accuracy_of_recovered": {},
        "ec_space_of_recovered": {
            m: {"n": len(v), "mean": round(st.mean(v), 4), "median": round(st.median(v), 4)}
            for m, v in ec_stats.items() if v},
        "by_region": {},
    }
    ids = [r["first_survivor_identity"] for r in answered if r["first_survivor_identity"]]
    if ids:
        out["first_survivor_identity"] = {
            "median": round(st.median(ids), 4), "mean": round(st.mean(ids), 4),
            "min": round(min(ids), 4), "max": round(max(ids), 4)}
    for rank in ("family", "genus"):
        cells = [r.get(f"{rank}_correct") for r in answered]
        scoreable = [c for c in cells if c is not None]
        out["anchor_accuracy_of_recovered"][rank] = {
            "n_scoreable": len(scoreable),
            "correct": round(sum(1 for c in scoreable if c) / len(scoreable), 4)
            if scoreable else None,
            "no_call_or_truth_undefined": len(cells) - len(scoreable)}
    per = defaultdict(list)
    for r in rows:
        per[r["region"]].append(r)
    for region, rs in sorted(per.items()):
        a = [r for r in rs if r["answered"]]
        g = [r.get("genus_correct") for r in a if r.get("genus_correct") is not None]
        out["by_region"][region] = {
            "n": len(rs), "answered": len(a),
            "rate": round(len(a) / len(rs), 4) if rs else None,
            "median_first_survivor_rank": st.median(
                [r["first_survivor_rank"] for r in rs if r["first_survivor_rank"]]) if any(
                r["first_survivor_rank"] for r in rs) else None,
            "genus_correct_of_recovered": round(sum(1 for c in g if c) / len(g), 4) if g else None,
        }
    with open(args.out, "w") as fh:
        json.dump({**out, "rows": rows}, fh, indent=1)
    print(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
