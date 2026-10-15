#!/usr/bin/env python
"""Growth / no-growth across the DEFINED MEDIA PANEL, per arm, against the truth model.

Why this is a separate script from `phenotype_fold.py`
-----------------------------------------------------
SPEC D11 replaced the task's stated "growth/no-growth across a defined media panel"
readout with the 54-medium core ATP-yield vector, on the evidence that the media panel
returns an all-zero vector for every arm.  `phenotype_fold.py` therefore only evaluates
growth on the COMPLETE medium.  That makes the degeneracy a claim rather than a measured
result on this fixture, so this script measures it directly: it builds the same
complete-medium-gapfilled models by the same `phenotype_arm` calls and runs
`growth_vector` over the curated KBase/ModelSEED minimal media.

It is a READOUT, not a second pipeline: every model is built by
`phenotype_arm.build_gapfilled(ec_set, template, id, complete_media)` -- identical
procedure for truth, prFBA and PICRUSt2, differing only in the input EC set, with
`annotate_with_rast=False` inside `build_gapfilled`.  Gapfilling is on the COMPLETE
medium only; no medium that is tested is ever gapfilled on (D11).

Reported per (genome, arm):
  n_grow_minimal      media with bio1 > 1e-5 out of the panel
  growth_complete     bio1 on the complete medium (the positive reference: if this is
                      0 the model is dead and the panel's zeros say nothing)
  sanity_gate         phenotype_arm.sanity_gate (biomass and ATPM on water+protons)
and per genome the arm-vs-truth binary Hamming over the panel.

A panel that is all-zero for every arm INCLUDING truth is a finding about the metric
(it has no dynamic range on draft-quality reconstructions), not a result about the
methods, and is reported as such.
"""
import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
P2_SCRIPTS = os.path.join(REPO, "region_validation", "picrust2_comparison", "scripts")
if P2_SCRIPTS not in sys.path:
    sys.path.insert(0, P2_SCRIPTS)

GROWTH_THRESH = 1e-5
_PA = None


def pa():
    global _PA
    if _PA is None:
        import phenotype_arm as _m
        _PA = _m
    return _PA


def jload(path):
    import orjson
    with open(path, "rb") as fh:
        return orjson.loads(fh.read())


def worker_init(msd):
    import logging
    logging.disable(logging.WARNING)
    m = pa()
    m.preload_biochemistry(msd)
    for k in ("gramneg", "grampos", "archaea", "core"):
        m.load_template(k)


def run_one(job):
    """One genome, every arm: gapfill on complete, then the media panel."""
    m = pa()
    t0 = time.time()
    gf_media, gf_cpds = m.complete_media()
    panel = m.load_kbase_media(prefix=job["prefix"] or None, minimal_only=True)
    out = {"src_genome_id": job["genome_id"], "template": job["template"],
           "domain": job["domain"], "n_media": len(panel), "arms": {}}
    vectors = {}
    for arm, ec in sorted(job["arms"].items()):
        t1 = time.time()
        try:
            mdl, gapsol, _atpc, _tests = m.build_gapfilled(
                sorted(ec), job["template"], "{}_{}".format(
                    job["genome_id"].replace(".", "_"), arm), gf_media)
            gate = m.sanity_gate(mdl)
            gc = m.growth_vector(mdl, {"complete": gf_cpds}, thresh=GROWTH_THRESH)
            gv = m.growth_vector(mdl, panel, thresh=GROWTH_THRESH)
            vectors[arm] = gv
            out["arms"][arm] = {
                "n_ec_input": len(ec), "n_gapfilled_rxn": len(gapsol),
                "growth_complete": round(float(gc.get("complete", 0.0)), 6),
                "n_grow_minimal": sum(1 for v in gv.values() if v > GROWTH_THRESH),
                "grow_minimal": sorted(k for k, v in gv.items() if v > GROWTH_THRESH),
                "max_growth_minimal": round(max(gv.values()) if gv else 0.0, 9),
                "sanity_gate": {k: (round(v, 6) if isinstance(v, float) else v)
                                for k, v in gate.items()},
                "status": "ok", "wall_seconds": round(time.time() - t1, 2)}
        except Exception as exc:                        # noqa: BLE001
            out["arms"][arm] = {"n_ec_input": len(ec), "status": "exception",
                                "status_reason": "{}:{}".format(
                                    type(exc).__name__, exc)[:300],
                                "wall_seconds": round(time.time() - t1, 2)}
    tv = vectors.get("truth")
    for arm, gv in vectors.items():
        if tv is None or arm == "truth":
            continue
        d, n = m.hamming(tv, gv, thresh=GROWTH_THRESH)
        out["arms"][arm]["vs_truth"] = {
            "growth_hamming": int(d), "n_media_compared": int(n),
            "lost": sorted(k for k in tv
                           if tv[k] > GROWTH_THRESH and gv.get(k, 0.0) <= GROWTH_THRESH),
            "gained": sorted(k for k in gv
                             if gv[k] > GROWTH_THRESH and tv.get(k, 0.0) <= GROWTH_THRESH)}
    out["wall_seconds"] = round(time.time() - t0, 2)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--per-genome", required=True,
                    help="phenotype_fold.py's per_genome.jsonl (defines the genomes, "
                         "arms and templates actually used, so the panel readout is "
                         "measured on exactly the same models)")
    ap.add_argument("--panel", required=True, help="panel_eligible.json")
    ap.add_argument("--fixture-meta", default=None,
                    help="fixture_meta.json (only for provenance in the output)")
    ap.add_argument("--functions", required=True,
                    help="genome_id -> {ec: [...]} (truth EC sets)")
    ap.add_argument("--prfba", required=True)
    ap.add_argument("--picrust2", required=True)
    ap.add_argument("--select-mode", default="exact_tie")
    ap.add_argument("--n-genomes", type=int, default=3)
    ap.add_argument("--prefix", default="Carbon-",
                    help="KBase media name prefix; '' = every curated minimal medium")
    ap.add_argument("--out", required=True)
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--msd", default="/home/freiburger/Documents/ModelSEEDDatabase")
    args = ap.parse_args(argv)

    prf_arm = "prfba_" + args.select_mode
    panel = jload(args.panel)["genomes"]
    functions = jload(args.functions)

    # the EC sets, read from the SAME predictions.json files phenotype_fold.py consumed
    import glob
    prf = {}
    for p in sorted(glob.glob(os.path.join(args.prfba, "*", "predictions.json"))):
        d = jload(p)
        blk = d.get("pred", {}).get(args.select_mode, {})
        prf.update(blk)
    p2 = {}
    for p in sorted(glob.glob(os.path.join(args.picrust2, "*", "predictions.json"))):
        d = jload(p)
        dom = d.get("meta", {}).get("domain")
        for md5, rec in d.get("pred", {}).items():
            p2[(dom, md5)] = rec

    seen, jobs = set(), []
    with open(args.per_genome) as fh:
        for line in fh:
            r = json.loads(line)
            gid = r["src_genome_id"]
            if gid in seen or r["arm"] not in ("truth", prf_arm, "picrust2"):
                continue
            if r.get("status") != "ok":
                continue
            md5, dom = r.get("md5"), r.get("domain")
            dk = "arc" if (dom or "").lower().startswith("arch") else "bac"
            arms = {"truth": sorted(functions[gid]["ec"])}
            if md5 in prf:
                arms[prf_arm] = sorted(prf[md5]["ec"])
            if (dk, md5) in p2:
                arms["picrust2"] = sorted(p2[(dk, md5)]["ec"])
            if len(arms) < 3:
                continue
            seen.add(gid)
            jobs.append({"genome_id": gid, "domain": dom, "md5": md5,
                         "template": r.get("template") or panel[gid]["template"],
                         "arms": arms, "prefix": args.prefix})
            if len(jobs) >= args.n_genomes:
                break
    print("[panel] {} genomes x {} arms, media prefix {!r}".format(
        len(jobs), 3, args.prefix), flush=True)
    if not jobs:
        raise SystemExit("[panel] no eligible genomes in --per-genome")

    import multiprocessing as mp
    ctx = mp.get_context("spawn")
    nw = max(1, min(args.workers, len(jobs)))
    results = []
    if nw == 1:
        worker_init(args.msd)
        for j in jobs:
            results.append(run_one(j))
            print("[panel] {} done".format(results[-1]["src_genome_id"]), flush=True)
    else:
        pool = ctx.Pool(nw, initializer=worker_init, initargs=(args.msd,))
        try:
            for res in pool.imap_unordered(run_one, jobs):
                results.append(res)
                print("[panel] {} done in {}s".format(res["src_genome_id"],
                                                      res["wall_seconds"]), flush=True)
        finally:
            pool.close()
            pool.join()
    results.sort(key=lambda r: r["src_genome_id"])

    n_media = results[0]["n_media"] if results else 0
    arms = sorted({a for r in results for a in r["arms"]})
    summary = {}
    for a in arms:
        rows = [r["arms"][a] for r in results if r["arms"].get(a, {}).get("status") == "ok"]
        summary[a] = {
            "n_models": len(rows),
            "n_grow_minimal_total": sum(x["n_grow_minimal"] for x in rows),
            "n_models_with_any_minimal_growth": sum(1 for x in rows
                                                    if x["n_grow_minimal"] > 0),
            "growth_complete_mean": (round(sum(x["growth_complete"] for x in rows)
                                           / len(rows), 6) if rows else None),
            "n_models_growing_on_complete": sum(1 for x in rows
                                                if x["growth_complete"] > GROWTH_THRESH),
            "n_sanity_gate_pass": sum(1 for x in rows if x["sanity_gate"]["pass"]),
            "growth_hamming_vs_truth_total": sum(
                x.get("vs_truth", {}).get("growth_hamming", 0) for x in rows)
            if a != "truth" else 0}
    out = {"meta": {"n_genomes": len(results), "n_media_in_panel": n_media,
                    "media_prefix": args.prefix,
                    "growth_threshold": GROWTH_THRESH,
                    "gapfill_medium": "complete (347 compounds, KBaseMedia.cpd)",
                    "arms": arms, "select_mode": args.select_mode,
                    "per_genome_source": os.path.abspath(args.per_genome),
                    "fixture_meta": (jload(args.fixture_meta).get("prediction_source")
                                     if args.fixture_meta
                                     and os.path.exists(args.fixture_meta) else None),
                    "script": "phenotype_media_panel.py",
                    "total_growth_calls": len(results) * len(arms) * n_media},
           "summary": summary, "per_genome": results}
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(out, fh, indent=1, sort_keys=True)
    print("[panel] wrote {}".format(args.out), flush=True)
    for a in arms:
        s = summary[a]
        print("[panel] {:18s} grows on {}/{} minimal media across {} models; "
              "complete-medium growth {}/{}".format(
                  a, s["n_grow_minimal_total"], len(results) * n_media, s["n_models"],
                  s["n_models_growing_on_complete"], s["n_models"]), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
