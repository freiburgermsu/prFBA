#!/usr/bin/env python
"""STUB Stage E weights for the Stage O smoke test (PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md).

NOT Stage E. Stage E (scripts/pfunc_estimators.py, a parallel task) selects the estimator
and writes the real file, data/pfunc_estimators/stage_o_weights.json. This script writes a
file with the same schema so the Stage O harness (scripts/objective_tuning.py) can be
smoke-tested before Stage E finishes. It refuses to write the real path, and it marks its
output `meta.stub = true`, which objective_tuning.py refuses unless --allow-stub is given.

The stub estimator is E1 at the calibration study's exploratory setting, with E3's fields
present:
  * records: each fold's own hits, data/prfba/f{fold:02d}_d100/hits/
    asv_top20_alignment_hits.json, streamed with ijson; each record truncated to its
    first --hits (100) hits, as the design says Stage E does (production records hold 100);
  * stage 0: calibrate_pfunc.stage0_ranked, imported (select_references stage 0, then
    identity / score ranking with genome dedup);
  * neighbours: the first --k (100) ranked survivors, w_j = exp(-(id_1 - id_j) / tau),
    tau = --tau (0.01);
  * clade: the fold's reference genomes (the keys of _cache_ec_subset.json) in the genus of
    the anchor (the first exact-tie genome, reps[0]); fewer than 3 backs off to family;
  * alpha = --alpha (0.0: the clade term is present but carries no weight);
  * meta.isotonic = null.

Checks, counted in meta.assertions (all should be 0 except the informational counts):
the source genome is never a neighbour or clade member; no neighbour or clade member shares
genus, species, taxon_id or holdout group with the source; every one is in the fold's
reference; the top-identity tie set among the stage-0 survivors equals the job's reps (an
exact-tie genome then has w = 1).

Usage:
  ~/Documents/py_venv/bin/python scripts/objective_tuning_stub_weights.py \
      --out data/objective_tuning/smoke/stage_o_weights.STUB.json
"""
import argparse
import hashlib
import json
import math
import os
import sys
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)
D = os.path.join(SR, "data")
REAL_PATH = os.path.abspath(os.path.join(D, "pfunc_estimators", "stage_o_weights.json"))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fold_job(arg):
    """One fold: stream its hits file once, rank and weight the population's records."""
    fold, jobs, k, tau, n_hits, alpha = arg
    import ijson
    import orjson
    sys.path.insert(0, HERE)
    from calibrate_pfunc import stage0_ranked          # noqa: E402  (select_references stage 0)
    t0 = time.time()
    with open(os.path.join(D, "prfba", "f{:02d}_d100".format(fold), "_cache_ec_subset.json"),
              "rb") as fh:
        ref = set(orjson.loads(fh.read()))
    with open(os.path.join(D, "panel_eligible.json"), "rb") as fh:
        panel = orjson.loads(fh.read())["genomes"]
    with open(os.path.join(D, "fold_plan.json"), "rb") as fh:
        ggroup = orjson.loads(fh.read())["genome_group"]
    by_md5 = {}
    for j in jobs:
        by_md5.setdefault(j["md5"], []).append(j)
    # the fold reference grouped by genus / family (the clade candidates)
    clade_idx = {"genus": {}, "family": {}}
    for g in ref:
        rec = panel.get(g) or {}
        for rank in ("genus", "family"):
            if rec.get(rank):
                clade_idx[rank].setdefault(rec[rank], []).append(g)
    A = Counter()
    out = {}
    path = os.path.join(D, "prfba", "f{:02d}_d100".format(fold), "hits",
                        "asv_top20_alignment_hits.json")
    with open(path, "rb") as fh:
        for md5, rec in ijson.kvitems(fh, "", use_float=True):
            if md5 not in by_md5 or not isinstance(rec, dict):
                continue
            rec = dict(rec)
            A["records_hits_total_gt_{}".format(n_hits)] += int(len(rec.get("top20") or [])
                                                               > n_hits)
            rec["top20"] = list(rec.get("top20") or [])[:n_hits]
            ranked, info = stage0_ranked(rec)
            for j in by_md5[md5]:
                src = j["src_genome_id"]
                srec = panel.get(src) or {}
                if not ranked:
                    A["no_stage0_survivors"] += 1
                    continue
                nb = ranked[:k]
                id1 = nb[0][1]
                neigh = [[g, math.exp(-(id1 - x) / tau)] for g, x in nb]
                tie = {g for g, x in ranked if x == id1}
                A["tie_set_ne_reps"] += int(tie != set(j["reps"]))
                A["reps_weight_ne_1"] += sum(1 for g, w in neigh if g in j["reps"] and w != 1.0)
                anchor = j["reps"][0]
                arec = panel.get(anchor) or {}
                rank_used, members = None, []
                for rank in ("genus", "family"):
                    cand = sorted(g for g in clade_idx[rank].get(arec.get(rank), [])
                                  if g != src)
                    if len(cand) >= 3:
                        rank_used, members = rank, cand
                        break
                A["anchor_genus_eq_source_genus"] += int(arec.get("genus") == srec.get("genus"))
                if rank_used is None:
                    A["clade_null_lt3_at_family"] += 1
                for role, gl in (("neighbour", [g for g, _ in neigh]), ("clade", members)):
                    A["leak_source_is_{}".format(role)] += int(src in gl)
                    A["{}_outside_fold_reference".format(role)] += sum(1 for g in gl
                                                                       if g not in ref)
                    for nm in ("genus", "species", "taxon_id"):
                        A["leak_{}_shares_{}".format(role, nm)] += sum(
                            1 for g in gl if (panel.get(g) or {}).get(nm) is not None and
                            (panel.get(g) or {}).get(nm) == srec.get(nm))
                    A["leak_{}_shares_group".format(role)] += sum(
                        1 for g in gl if ggroup.get(g) is not None and
                        ggroup.get(g) == ggroup.get(src))
                A["n_neighbours_lt_k"] += int(len(neigh) < k)
                out[j["key"]] = {
                    "md5": md5, "fold": fold, "neighbours": neigh,
                    "clade": ({"rank": rank_used, "members": members} if rank_used else None),
                    "alpha": float(alpha), "n_reps": j["n_reps"], "reps": list(j["reps"]),
                    "stub_info": {"n_survivors": info["n_survivors"], "n_hits_used":
                                  info["n_hits"], "anchor": anchor}}
    missing = sorted(j["key"] for j in jobs if j["key"] not in out)
    A["jobs_without_record"] = len(missing)
    return fold, out, dict(A), missing, round(time.time() - t0, 1)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--population", default=os.path.join(D, "weighted_objective",
                                                         "population.json"))
    ap.add_argument("--out", required=True)
    ap.add_argument("--k", type=int, default=100)
    ap.add_argument("--tau", type=float, default=0.01)
    ap.add_argument("--hits", type=int, default=100)
    ap.add_argument("--alpha", type=float, default=0.0)
    ap.add_argument("--only", default=None, metavar="KEY,KEY")
    args = ap.parse_args(argv)
    out_abs = os.path.abspath(args.out)
    if out_abs == REAL_PATH or os.path.dirname(out_abs) == os.path.dirname(REAL_PATH):
        raise SystemExit("[stub] refusing to write under data/pfunc_estimators/: that is "
                         "Stage E's output")
    if ".STUB." not in os.path.basename(out_abs):
        raise SystemExit("[stub] the output name must contain '.STUB.'")
    with open(args.population) as fh:
        pop = json.load(fh)
    jobs = pop["jobs"]
    if args.only:
        keep = set(args.only.split(","))
        jobs = [j for j in jobs if j["key"] in keep]
    by_fold = {}
    for j in jobs:
        by_fold.setdefault(int(j["fold"]), []).append(
            {k: j[k] for k in ("key", "md5", "fold", "reps", "n_reps", "src_genome_id")})
    import multiprocessing as mp
    t0 = time.time()
    tasks = [(f, by_fold[f], args.k, args.tau, args.hits, args.alpha) for f in sorted(by_fold)]
    with mp.get_context("fork").Pool(len(tasks)) as pool:
        res = pool.map(fold_job, tasks)
    amps, A, missing, secs = {}, Counter(), [], {}
    for fold, out, a, miss, s in res:
        amps.update(out)
        A.update(a)
        missing += miss
        secs[fold] = s
    meta = {"stub": True, "warning": "SMOKE-TEST STUB, NOT STAGE E OUTPUT",
            "generator": {"path": os.path.abspath(__file__), "sha256": sha256(__file__)},
            "created": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "estimator": "E1 hits (stub): K={}, tau={}, dedup=genome, first {} hits, "
                         "stage 0 = calibrate_pfunc.stage0_ranked; clade = anchor (reps[0]) "
                         "genus in the fold reference, family if < 3; alpha={}".format(
                             args.k, args.tau, args.hits, args.alpha),
            "K": args.k, "tau": args.tau, "dedup": "genome", "hits_used": args.hits,
            "alpha": args.alpha, "isotonic": None,
            "population": {"path": os.path.abspath(args.population),
                           "sha256": sha256(args.population)},
            "n_amplicons": len(amps), "n_jobs": len(jobs), "missing": missing,
            "assertions": dict(sorted(A.items())), "seconds_by_fold": secs,
            "seconds": round(time.time() - t0, 1)}
    os.makedirs(os.path.dirname(out_abs), exist_ok=True)
    with open(out_abs + ".tmp", "w") as fh:
        json.dump({"meta": meta, "amplicons": amps}, fh, sort_keys=True)
    os.replace(out_abs + ".tmp", out_abs)
    print(json.dumps({k: meta[k] for k in ("n_amplicons", "n_jobs", "assertions", "seconds",
                                           "seconds_by_fold")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
