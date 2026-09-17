#!/usr/bin/env python
"""Re-run the shipped prFBA reference selection for every 10k-benchmark amplicon.

For each ``genome__region`` amplicon this applies ``select_references.
select_representatives`` with the current production defaults (exact-tie
de-duplicated union, BV-BRC species dedup, P3 gene-richest tie-break fed from the
PGFam cache) to the amplicon's top-20 alignment hits, under two regimes:

* ``include``        the hits as-is (source genome may be in the database)
* ``exclude_0.987``  ``score_genecap.exclude_self_record`` first drops hits that are
                     the source genome/taxon or >= 0.987 identity, then re-selects

These are the reference sets whose union forms prFBA's synthetic genome. The
committed manuscript tables used the legacy reducer; ``--legacy`` reproduces that.
No network access: gene counts come from the local PGFam cache and the gene
provider is only consulted by the legacy reducer, where it is cache-only here.

Output: ``--out`` JSON  {amplicon_key: {regime: {tier, n_selected, reps[]}}, _meta}
"""
import argparse
import os
import subprocess
import sys
import time
from multiprocessing import get_context

import orjson

HERE = os.path.dirname(os.path.abspath(__file__))
RV = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RV, "scripts"))
import _common as C  # noqa: E402
import score_genecap as SG  # noqa: E402  (also puts prFBA root on sys.path)
import select_references as SR  # noqa: E402

PATHS = C.PATHS
_G = {}


def _cache_only_provider(cache):
    """gene_provider backed solely by the on-disk PGFam cache -- never fetches.

    The cache holds lists; the legacy reducer unions them into sets, so convert on
    first use (the union default never calls the provider at all).
    """
    def provider(gid):
        gid = str(gid)
        v = cache.get(gid)
        if v is None:
            return set()
        if not isinstance(v, set):
            v = set(v)
            cache[gid] = v
        return v
    provider.cache = cache
    return provider


def _reps(sel):
    return {"tier": sel.get("tier"), "n_selected": sel.get("n_selected", 0),
            "reps": SG.selected_rep_gids(sel)}


def _work(keys):
    hits, truth, provider, gcf, knobs = (_G["hits"], _G["truth"], _G["provider"],
                                         _G["gcf"], _G["knobs"])
    out = {}
    for key in keys:
        hr = hits.get(key)
        src_gid, _ = SG.split_key(key)
        src_tid = truth[key]["src_taxon_id"]
        if hr is None or "top20" not in hr:
            empty = {"tier": None, "n_selected": 0, "reps": [], "no_hits": True}
            out[key] = {SG.MODE_INCLUDE: empty, SG.MODE_EXCLUDE: dict(empty),
                        "best_nonself_identity": None}
            continue
        inc = SR.select_representatives(hr, gene_provider=provider, gene_count_fn=gcf,
                                        knobs=knobs)
        excl_rec, n_removed = SG.exclude_self_record(hr, src_gid, src_tid, SG.TAU)
        exc = SR.select_representatives(excl_rec, gene_provider=provider,
                                        gene_count_fn=gcf, knobs=knobs)
        out[key] = {
            SG.MODE_INCLUDE: {**_reps(inc), "best_identity": hr.get("best_identity")},
            SG.MODE_EXCLUDE: {**_reps(exc), "n_removed_self": n_removed,
                              "best_identity": excl_rec.get("best_identity")},
            # distance to the nearest genuinely different organism (the novelty axis
            # score_genecap bins on), independent of either regime's selection
            "best_nonself_identity": SG.best_nonself_identity(hr, src_gid, src_tid),
        }
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=32)
    ap.add_argument("--limit", type=int, default=0, help="first N amplicons (testing)")
    ap.add_argument("--legacy", action="store_true",
                    help="legacy ordered reducer (select_all=False), as in the manuscript")
    args = ap.parse_args()

    t0 = time.time()
    truth = orjson.loads(open(PATHS.truth_json, "rb").read())
    hits = orjson.loads(open(PATHS.hits_json, "rb").read())
    hits = {k: v for k, v in hits.items() if k != "_meta" and isinstance(v, dict)}
    print(f"[prfba_selection] truth={len(truth)} hits={len(hits)} ({time.time()-t0:.0f}s)")
    cache = orjson.loads(open(PATHS.pgfam_cache, "rb").read())
    print(f"[prfba_selection] PGFam cache: {len(cache)} genomes ({time.time()-t0:.0f}s)")
    provider = _cache_only_provider(cache)
    _G.update(hits=hits, truth=truth, provider=provider,
              gcf=SR.gene_count_from_provider(provider),
              knobs={"select_all": False} if args.legacy else None)

    keys = sorted(truth)
    if args.limit:
        keys = keys[:args.limit]
    chunks = [keys[i:i + 500] for i in range(0, len(keys), 500)]
    result = {}
    with get_context("fork").Pool(args.workers) as pool:
        for i, part in enumerate(pool.imap_unordered(_work, chunks), 1):
            result.update(part)
            if i % 20 == 0:
                print(f"[prfba_selection] {len(result)}/{len(keys)} ({time.time()-t0:.0f}s)")

    reps = {g for rec in result.values() for m in rec.values()
            if isinstance(m, dict) for g in m["reps"]}
    srcs = {SG.split_key(k)[0] for k in result}
    commit = subprocess.run(["git", "-C", RV, "rev-parse", "HEAD"], capture_output=True,
                            text=True).stdout.strip()
    result["_meta"] = {
        "select_mode": "legacy" if args.legacy else "default (exact-tie union)",
        "knobs": {**SR.DEFAULT_KNOBS, **(_G["knobs"] or {})},
        "p3_gene_count_fn": "PGFam cache family counts",
        "exclude_tau": SG.TAU,
        "n_amplicons": len(keys),
        "n_distinct_reps": len(reps),
        "n_reps_missing_from_pgfam_cache": len(reps - set(cache)),
        "n_sources_missing_from_pgfam_cache": len(srcs - set(cache)),
        "git_commit": commit,
        "inputs": {"hits": PATHS.hits_json, "truth": PATHS.truth_json,
                   "pgfam_cache": PATHS.pgfam_cache},
    }
    with open(args.out, "wb") as fh:
        fh.write(orjson.dumps(result))
    print(f"[prfba_selection] wrote {args.out}: {orjson.dumps(result['_meta']).decode()}"
          f" ({time.time()-t0:.0f}s)")


if __name__ == "__main__":
    main()
