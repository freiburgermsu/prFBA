#!/usr/bin/env python
"""vocab.py -- part 4: what dropping the RefSeq features does to the models themselves.

Every build up to prFBA 3ce12d0 merged two annotations of each source genome: BV-BRC's own
PATRIC (RASTtk) CDSs and the imported RefSeq (PGAP) CDSs. The BV-BRC-only rebuild keeps PATRIC
only. Parts 2 and 3 compare k/n with E3 within one vocabulary; this part compares the two
vocabularies under k/n, on the part-3 sample (80 n >= 2 genomes and 10 n = 1 controls).

Arms (both k/n; same draft rule, templates and gapfill):
  against (--against, default mixed)  the shipped k/n genome objects
  builds  (--builds, default bvbrc)   the BV-BRC-only k/n genome objects
Each arm's k/n draft, gapfilled model and per-genome record come from that build set's part-3
run (flux.py --builds <set>): build_prfba_community_models._build_one on ModelSEEDTemplates
21a8084 (c12run), then ONE gapfill on the complete medium the way Stage O gapfills (MSATPCorrection,
MSGapfill on the KBase complete medium at uptake 100, reversals applied). Both runs import the
same pinned ModelSEEDpy / MSCommunity and the same codiffusion files (hashes in each flux.json).

Checks
  * the two part-3 samples are the same genomes;
  * the template class of every genome is the same in both;
  * reproducibility of the reused drafts: the 'against' k/n draft is rebuilt here from its genome
    object with flux.py's draft rule (into <models>/part4/{dep}/{against}/) and must equal the
    saved one structurally (ids, stoichiometry, bounds, GPR as a DNF set).

Per genome
  features   per arm; 'against' features whose function string is absent from 'builds' (the
             RefSeq-only products), and the reverse (must be 0: BV-BRC-only is a subset)
  draft      reaction ids gained / lost; gene-backed reactions per arm; gene-backed gained / lost;
             GPRs changed on shared gene-backed reactions, compared as DNF sets over FUNCTION
             strings (feature ids are renumbered between the builds; one merged feature per
             function, so the function is the gene's identity), classed by gene set:
             fewer genes / more genes / same genes, other structure / other;
             reaction p (notes, k/n) changed on shared reactions
  growth     max bio1 on the complete medium (uptake 100, C <= 100) of the draft (before any
             gapfill; computed here) and of the gapfilled model (the part-3 record)
  burden     MSGapfill new reactions, reversals applied, reactions ATP correction added (in the
             gapfilled model, not in the draft, not in the MSGapfill solution), draft reactions
             absent after gapfilling, gene-less share of the gapfilled model's rxn reactions,
             and how many of the draft reactions the BV-BRC-only arm lost were added back
             gene-less by its gapfill

  run via systemd-run (README); outputs <results>/part4/{dep}/{asv}.json, <results>/vocab.json
"""
from __future__ import annotations

import argparse
import json
import os
import statistics
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402
import flux as F  # noqa: E402  (draft rule, gapfill rule, environment; paths of the --builds set)

EPS = 1e-9
OUT4 = f"{C.RESULTS}/part4"
MODELS4 = f"{C.MODELS_ROOT}/part4"


def arm_paths(name, dep, asv):
    spec = C.BUILD_SETS[name]
    return {"genome": f"{spec['builds'][dep]['kn']}/genome_objects/{asv}.json",
            "draft": f"{spec['models']}/part3/{dep}/kn/models_{F.TAG}/{asv}.xml",
            "gapfilled": f"{spec['models']}/part3/{dep}/gapfilled/{asv}.kn.xml",
            "record": f"{spec['results']}/part3/{dep}/{asv}.json",
            "assignment": f"{spec['results']}/part3/_work/{dep}_template_assignment.json"}


def fkey(feature):
    return " / ".join(feature.get("functions") or [])


def gpr_by_function(rule, fn_of):
    can = C.gpr_canonical(rule)
    if isinstance(can, tuple):
        return can
    return frozenset(frozenset(fn_of.get(g, f"?{g}") for g in clause) for clause in can)


# ------------------------------------------------------------------ worker
_V = {}


def _init(against):
    F._init_once()
    _V["against"] = against


def _rebuild_draft(dep, asv, against):
    """flux.py's draft rule (build_prfba_community_models._build_one) on the 'against' k/n
    genome object, into this part's own directory, with the 'against' template assignment."""
    gdir = C.BUILD_SETS[against]["builds"][dep]["kn"]
    B = F.import_B(dep, gdir)
    B.TPL_ASSIGN_PATH = arm_paths(against, dep, asv)["assignment"]
    d = f"{MODELS4}/{dep}/{against}"
    os.makedirs(f"{d}/models_{F.TAG}", exist_ok=True)
    os.chdir(d)
    B._W.clear()
    B._binit()
    _, res = B._build_one((asv, asv))
    return f"{d}/models_{F.TAG}/{asv}.xml", res


def run_one(job):
    from cobra.io import read_sbml_model
    dep, asv, against, builds = job["dep"], job["asv"], _V["against"], C.BUILD_SET
    t0 = time.time()
    rec = {"dep": dep, "asv": asv, "n": job["n"], "stratum": job["stratum"], "control": job["control"]}
    try:
        P = {against: arm_paths(against, dep, asv), builds: arm_paths(builds, dep, asv)}
        rec["inputs"] = P
        gen = {a: json.load(open(P[a]["genome"])) for a in P}
        recs = {a: json.load(open(P[a]["record"])) for a in P}
        for a in P:
            if recs[a].get("status") != "ok":
                raise RuntimeError(f"part-3 record of {a} is not ok: {recs[a].get('error')}")
        rec["template"] = {a: recs[a]["template"] for a in P}
        # features
        fn_of = {a: {f["id"]: fkey(f) for f in gen[a]["features"]} for a in P}
        fset = {a: set(fn_of[a].values()) for a in P}
        rec["features"] = {
            "n": {a: len(gen[a]["features"]) for a in P},
            "n_functions": {a: len(fset[a]) for a in P},
            "n_with_patric_id": {a: sum(1 for f in gen[a]["features"]
                                        if any(x[0] == "PATRIC_id" and x[1] for x in f.get("aliases", [])))
                                 for a in P},
            "n_against_only": len(fset[against] - fset[builds]),
            "n_builds_only": len(fset[builds] - fset[against])}
        # drafts: the saved ones, plus a fresh rebuild of the 'against' one
        t = time.time()
        fresh_path, fresh_res = _rebuild_draft(dep, asv, against)
        rec["fresh_draft"] = {"result": fresh_res, "secs": round(time.time() - t, 1)}
        dr = {a: read_sbml_model(P[a]["draft"]) for a in P}
        fresh = read_sbml_model(fresh_path)
        sa = {r.id: F._sig(r) for r in dr[against].reactions}
        sf = {r.id: F._sig(r) for r in fresh.reactions}
        rec["fresh_draft"].update({"same_reaction_ids": set(sa) == set(sf),
                                   "n_struct_diff": sum(1 for i in sa if i in sf and sa[i] != sf[i])})
        R = {a: {r.id: r for r in dr[a].reactions} for a in P}
        gb = {a: {i for i, r in R[a].items() if r.genes} for a in P}
        shared = set(R[against]) & set(R[builds])
        gpr_cls = {"fewer_genes": 0, "more_genes": 0, "same_genes_other_structure": 0, "other": 0}
        gpr_examples = []
        n_gpr_changed = 0
        for i in sorted(shared & gb[against] & gb[builds]):
            ga = gpr_by_function(R[against][i].gene_reaction_rule, fn_of[against])
            gbb = gpr_by_function(R[builds][i].gene_reaction_rule, fn_of[builds])
            if ga == gbb:
                continue
            n_gpr_changed += 1
            ua = set().union(*ga) if not isinstance(ga, tuple) else {ga}
            ub = set().union(*gbb) if not isinstance(gbb, tuple) else {gbb}
            cls = ("fewer_genes" if ub < ua else "more_genes" if ub > ua else
                   "same_genes_other_structure" if ub == ua else "other")
            gpr_cls[cls] += 1
            if len(gpr_examples) < 3:
                gpr_examples.append({"reaction": i, "n_genes": [len(ua), len(ub)],
                                     "functions_lost": sorted(ua - ub)[:4]})
        pa = {i: float(R[against][i].notes.get("probability", "nan")) for i in shared}
        pb = {i: float(R[builds][i].notes.get("probability", "nan")) for i in shared}
        dp = [pb[i] - pa[i] for i in shared if abs(pb[i] - pa[i]) > EPS]
        lost = set(R[against]) - set(R[builds])
        rec["draft"] = {
            "n_rxn": {a: len(R[a]) for a in P},
            "n_rxn_prefixed": {a: sum(1 for i in R[a] if i.startswith("rxn")) for a in P},
            "n_gene_backed": {a: len(gb[a]) for a in P},
            "n_genes": {a: len(dr[a].genes) for a in P},
            "n_gained": len(set(R[builds]) - set(R[against])),
            "n_lost": len(lost),
            "n_gene_backed_gained": len(gb[builds] - gb[against]),
            "n_gene_backed_lost": len(gb[against] - gb[builds]),
            "n_gene_backed_lost_reaction_kept_geneless": len((gb[against] - gb[builds]) & set(R[builds])),
            "n_shared_gene_backed": len(shared & gb[against] & gb[builds]),
            "n_gpr_changed": n_gpr_changed, "gpr_change_class": gpr_cls, "gpr_examples": gpr_examples,
            "n_p_changed_shared": len(dp), "n_p_down_shared": sum(1 for x in dp if x < 0),
            "mean_dp_changed_shared": statistics.mean(dp) if dp else 0.0,
            "lost_reactions": sorted(lost)}
        rec["draft_growth_complete_c100"] = {a: F._max_growth(dr[a]) for a in P}
        # gapfilled models
        gf = {a: read_sbml_model(P[a]["gapfilled"]) for a in P}
        G = {a: {r.id for r in gf[a].reactions} for a in P}
        burden = {}
        for a in P:
            new = {f"{x}0" for x in recs[a]["gapfill"]["new_rxn"]}       # template id -> _c0
            added = G[a] - set(R[a])
            rx = [r for r in gf[a].reactions if r.id.startswith("rxn")]
            burden[a] = {"growth_complete_c100": recs[a]["max_growth_complete_c100"],
                         "infeasible": recs[a]["gapfill"]["infeasible"],
                         "n_new_rxn": recs[a]["gapfill"]["n_new_rxn"],
                         "n_reversed_applied": recs[a]["gapfill"]["n_reversed_applied"],
                         "n_added_total": len(added),
                         "n_added_not_in_gapfill_solution": len(added - new),
                         "n_new_rxn_not_in_model": len(new - G[a]),
                         "n_draft_rxn_absent_after": len(set(R[a]) - G[a]),
                         "n_rxn_prefixed": len(rx),
                         "geneless_share_rxn": (sum(1 for r in rx if not r.genes) / len(rx)) if rx else None}
        burden["lost_draft_rxn_readded_by_gapfill"] = len(lost & G[builds])
        burden["lost_draft_rxn_readded_in_gapfill_solution"] = len(
            lost & {f"{x}0" for x in recs[builds]["gapfill"]["new_rxn"]})
        burden["lost_draft_rxn_in_against_gapfilled_model"] = len(lost & G[against])
        rec["gapfilled"] = burden
        rec["status"] = "ok"
    except Exception as exc:  # noqa: BLE001
        import traceback
        rec["status"] = "error"
        rec["error"] = f"{type(exc).__name__}: {str(exc)[:400]}"
        rec["traceback"] = traceback.format_exc()[-2000:]
    rec["secs"] = round(time.time() - t0, 1)
    C.dump(rec, f"{OUT4}/{dep}/{asv}.json")
    return {k: rec.get(k) for k in ("dep", "asv", "status", "error", "secs")}


# ------------------------------------------------------------------ summary
def summarise(rows, against, builds):
    d = F.describe
    ok = [r for r in rows if r["status"] == "ok"]
    out = {"n_genomes": len(rows), "n_ok": len(ok),
           "errors": [{k: r.get(k) for k in ("dep", "asv", "error")} for r in rows if r["status"] != "ok"],
           "checks": {
               "template_class_equal": sum(1 for r in ok if r["template"][against] == r["template"][builds]),
               "fresh_against_draft_equal": sum(1 for r in ok if r["fresh_draft"]["same_reaction_ids"]
                                                and r["fresh_draft"]["n_struct_diff"] == 0),
               "builds_functions_subset_of_against": sum(1 for r in ok if r["features"]["n_builds_only"] == 0),
               "n_ok": len(ok)}}
    groups = (("n>=2", lambda r: not r["control"]), ("n=1 controls", lambda r: r["control"]),
              ("all", lambda r: True))
    for gname, sel in groups:
        for dep in list(C.BUILDS) + ["pooled"]:
            R = [r for r in ok if sel(r) and (dep == "pooled" or r["dep"] == dep)]
            if not R:
                continue
            ratio = lambda f: [f(r) for r in R]
            blk = {"n_genomes": len(R)}
            blk["features"] = {a: d([r["features"]["n"][a] for r in R]) for a in (against, builds)}
            blk["features_ratio"] = d(ratio(lambda r: r["features"]["n"][builds] / r["features"]["n"][against]))
            blk["features_against_only_share"] = d(ratio(
                lambda r: r["features"]["n_against_only"] / r["features"]["n_functions"][against]))
            blk["n_genomes_features_unchanged"] = sum(1 for r in R if r["features"]["n_against_only"] == 0)
            dr = lambda k: [r["draft"][k] for r in R]
            blk["draft"] = {
                "n_rxn": {a: d([r["draft"]["n_rxn"][a] for r in R]) for a in (against, builds)},
                "n_gene_backed": {a: d([r["draft"]["n_gene_backed"][a] for r in R]) for a in (against, builds)},
                "gene_backed_ratio": d(ratio(lambda r: r["draft"]["n_gene_backed"][builds]
                                             / r["draft"]["n_gene_backed"][against]
                                             if r["draft"]["n_gene_backed"][against] else None)),
                "n_gained": d(dr("n_gained")), "n_lost": d(dr("n_lost")),
                "n_gene_backed_lost": d(dr("n_gene_backed_lost")),
                "n_gene_backed_gained": d(dr("n_gene_backed_gained")),
                "n_gene_backed_lost_reaction_kept_geneless": d(dr("n_gene_backed_lost_reaction_kept_geneless")),
                "n_gpr_changed": d(dr("n_gpr_changed")),
                "gpr_changed_share_of_shared_gene_backed": d(ratio(
                    lambda r: r["draft"]["n_gpr_changed"] / r["draft"]["n_shared_gene_backed"]
                    if r["draft"]["n_shared_gene_backed"] else None)),
                "gpr_change_class_total": {c: sum(r["draft"]["gpr_change_class"][c] for r in R)
                                           for c in R[0]["draft"]["gpr_change_class"]},
                "n_p_changed_shared": d(dr("n_p_changed_shared")),
                "mean_dp_changed_shared": d([r["draft"]["mean_dp_changed_shared"] for r in R
                                             if r["draft"]["n_p_changed_shared"]]),
                "totals": {"gained": sum(dr("n_gained")), "lost": sum(dr("n_lost")),
                           "gene_backed_against": sum(r["draft"]["n_gene_backed"][against] for r in R),
                           "gene_backed_builds": sum(r["draft"]["n_gene_backed"][builds] for r in R),
                           "gpr_changed": sum(dr("n_gpr_changed")),
                           "shared_gene_backed": sum(dr("n_shared_gene_backed"))},
                "n_genomes_identical_reaction_sets": sum(1 for r in R if not r["draft"]["n_lost"]
                                                         and not r["draft"]["n_gained"]),
                "n_genomes_no_gpr_change": sum(1 for r in R if not r["draft"]["n_gpr_changed"])}
            dg = {a: [r["draft_growth_complete_c100"][a] for r in R] for a in (against, builds)}
            blk["draft_growth"] = {a: {"n_grows": sum(1 for x in dg[a] if x and x > 1e-6),
                                       "growth": d(dg[a])} for a in dg}
            gf = {a: [r["gapfilled"][a] for r in R] for a in (against, builds)}
            blk["gapfilled"] = {
                a: {"n_grows": sum(1 for x in gf[a] if (x["growth_complete_c100"] or 0) > 1e-6),
                    "n_infeasible": sum(1 for x in gf[a] if x["infeasible"]),
                    "growth": d([x["growth_complete_c100"] for x in gf[a]]),
                    "n_new_rxn": d([x["n_new_rxn"] for x in gf[a]]),
                    "n_new_rxn_sum": sum(x["n_new_rxn"] for x in gf[a]),
                    "n_reversed_applied_sum": sum(x["n_reversed_applied"] for x in gf[a]),
                    "n_added_total": d([x["n_added_total"] for x in gf[a]]),
                    "n_added_total_sum": sum(x["n_added_total"] for x in gf[a]),
                    "n_added_not_in_gapfill_solution": d([x["n_added_not_in_gapfill_solution"] for x in gf[a]]),
                    "n_draft_rxn_absent_after": d([x["n_draft_rxn_absent_after"] for x in gf[a]]),
                    "geneless_share_rxn": d([x["geneless_share_rxn"] for x in gf[a]])}
                for a in gf}
            blk["gapfilled"]["growth_relative_change"] = d(
                [(r["gapfilled"][builds]["growth_complete_c100"] - r["gapfilled"][against]["growth_complete_c100"])
                 / r["gapfilled"][against]["growth_complete_c100"] for r in R
                 if r["gapfilled"][against]["growth_complete_c100"]
                 and r["gapfilled"][builds]["growth_complete_c100"] is not None])
            blk["gapfilled"]["n_new_rxn_change"] = d([r["gapfilled"][builds]["n_new_rxn"]
                                                      - r["gapfilled"][against]["n_new_rxn"] for r in R])
            blk["gapfilled"]["n_genomes_more_gapfill"] = sum(
                1 for r in R if r["gapfilled"][builds]["n_new_rxn"] > r["gapfilled"][against]["n_new_rxn"])
            blk["gapfilled"]["n_genomes_less_gapfill"] = sum(
                1 for r in R if r["gapfilled"][builds]["n_new_rxn"] < r["gapfilled"][against]["n_new_rxn"])
            blk["gapfilled"]["lost_draft_rxn_readded_by_gapfill_sum"] = sum(
                r["gapfilled"]["lost_draft_rxn_readded_by_gapfill"] for r in R)
            blk["gapfilled"]["lost_draft_rxn_readded_in_gapfill_solution_sum"] = sum(
                r["gapfilled"]["lost_draft_rxn_readded_in_gapfill_solution"] for r in R)
            out.setdefault(gname, {})[dep] = blk
    out["by_stratum"] = {}
    for st in F.STRATA + ["1"]:
        R = [r for r in ok if r["stratum"] == st]
        out["by_stratum"][st] = {
            "n_genomes": len(R),
            "features_ratio": d([r["features"]["n"][builds] / r["features"]["n"][against] for r in R]),
            "gene_backed_ratio": d([r["draft"]["n_gene_backed"][builds] / r["draft"]["n_gene_backed"][against]
                                    for r in R if r["draft"]["n_gene_backed"][against]]),
            "n_lost": d([r["draft"]["n_lost"] for r in R]),
            "n_gpr_changed": d([r["draft"]["n_gpr_changed"] for r in R]),
            "n_new_rxn": {a: d([r["gapfilled"][a]["n_new_rxn"] for r in R]) for a in (against, builds)}}
    out["per_genome"] = []
    for r in sorted(rows, key=lambda r: (r["dep"], r["control"], r["stratum"], r["asv"])):
        row = {k: r.get(k) for k in ("dep", "asv", "n", "stratum", "control", "status")}
        if r["status"] == "ok":
            row.update({"template": r["template"][builds],
                        "features": [r["features"]["n"][against], r["features"]["n"][builds]],
                        "gene_backed": [r["draft"]["n_gene_backed"][against], r["draft"]["n_gene_backed"][builds]],
                        "lost": r["draft"]["n_lost"], "gained": r["draft"]["n_gained"],
                        "gpr_changed": r["draft"]["n_gpr_changed"],
                        "draft_growth": [r["draft_growth_complete_c100"][against],
                                         r["draft_growth_complete_c100"][builds]],
                        "growth": [r["gapfilled"][against]["growth_complete_c100"],
                                   r["gapfilled"][builds]["growth_complete_c100"]],
                        "n_new_rxn": [r["gapfilled"][against]["n_new_rxn"], r["gapfilled"][builds]["n_new_rxn"]],
                        "n_added_total": [r["gapfilled"][against]["n_added_total"],
                                          r["gapfilled"][builds]["n_added_total"]],
                        "lost_readded": r["gapfilled"]["lost_draft_rxn_readded_by_gapfill"]})
        else:
            row["error"] = r.get("error")
        out["per_genome"].append(row)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--builds", required=True, choices=sorted(C.BUILD_SETS),
                    help="the build set measured (its part-3 run must be complete); results go to its dir")
    ap.add_argument("--against", default="mixed", choices=sorted(C.BUILD_SETS),
                    help="the reference build set (its part-3 run must be complete)")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--summary-only", action="store_true")
    a = ap.parse_args()
    assert a.builds == C.BUILD_SET and a.against != a.builds, (a.builds, C.BUILD_SET, a.against)
    F.setup_env()
    C.ensure_pins()
    t_all = time.time()
    samples = {n: json.load(open(f"{C.BUILD_SETS[n]['results']}/part3/_work/sample.json"))
               for n in (a.against, a.builds)}
    key = lambda smp: {d: [(g["asv"], g["n"], g["stratum"], g["control"]) for g in v["genomes"]]
                       for d, v in smp.items()}
    if key(samples[a.against]) != key(samples[a.builds]):
        raise SystemExit("the two part-3 samples differ")
    sample = samples[a.builds]
    jobs = [{"dep": dep, **g} for dep, v in sample.items() for g in v["genomes"]]
    if not a.summary_only:
        todo = [j for j in jobs if not (os.path.exists(f"{OUT4}/{j['dep']}/{j['asv']}.json") and
                                        json.load(open(f"{OUT4}/{j['dep']}/{j['asv']}.json")).get("status") == "ok")]
        C.log(f"part 4: {len(jobs)} genomes, {len(todo)} to run, {a.workers} workers")
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(a.workers, initializer=_init, initargs=(a.against,),
                                          maxtasksperchild=10) as pool:
            for i, r in enumerate(pool.imap_unordered(run_one, todo, chunksize=1), 1):
                C.log(f"  [{i}/{len(todo)}] {r['dep']} {r['asv']} {r['status']} {r.get('error') or ''} ({r['secs']} s)")
    rows = [json.load(open(f"{OUT4}/{j['dep']}/{j['asv']}.json")) for j in jobs
            if os.path.exists(f"{OUT4}/{j['dep']}/{j['asv']}.json")]
    summ = summarise(rows, a.against, a.builds)
    flux_meta = {n: json.load(open(f"{C.BUILD_SETS[n]['results']}/flux.json"))["_meta"]["code"]
                 for n in (a.against, a.builds) if os.path.exists(f"{C.BUILD_SETS[n]['results']}/flux.json")}
    prev = None
    if a.summary_only and os.path.exists(f"{C.RESULTS}/vocab.json"):
        prev = json.load(open(f"{C.RESULTS}/vocab.json"))["_meta"].get("secs_total_run")
    res = {"_meta": {"script": "vocab.py", "code": C.code_record([os.path.abspath(__file__),
                                                                  os.path.join(HERE, "flux.py")]),
                     "pins": C.pin_record(), "templates_root": C.TPL_ROOT["c12run"],
                     "arms": {"against": {"build_set": a.against, **{k: v for k, v in C.BUILD_SETS[a.against].items()}},
                              "builds": {"build_set": a.builds, **{k: v for k, v in C.BUILD_SETS[a.builds].items()}}},
                     "part3_code_of_each_arm": flux_meta,
                     "sample": "the part-3 sample (identical in both build sets)",
                     "definitions": {
                         "gpr_changed": "a reaction gene-backed in both drafts whose GPR, as a DNF set over the "
                                        "genes' function strings, differs",
                         "gpr_change_class": "fewer_genes: the BV-BRC-only gene set is a strict subset; more_genes: "
                                             "strict superset; same_genes_other_structure; other",
                         "growth_complete_c100": "max bio1, KBase complete medium at uptake 100, C uptake <= 100 "
                                                 "(flux._max_growth), on the gapfilled model (part-3 record)",
                         "n_added_total": "reactions in the gapfilled model that are not in the draft (MSGapfill "
                                          "solution + ATP-correction gapfilling)",
                         "lost_draft_rxn_readded_by_gapfill": "draft reactions of the reference arm that the "
                                                              "BV-BRC-only draft lacks and its gapfilled model has"},
                     "secs_total_run": prev if a.summary_only else round(time.time() - t_all, 1),
                     "summary_written": time.strftime("%Y-%m-%dT%H:%M:%S")},
           **summ}
    C.dump(res, f"{C.RESULTS}/vocab.json")
    C.log(f"part 4 done: {summ['n_ok']}/{summ['n_genomes']} ok -> {C.RESULTS}/vocab.json")


if __name__ == "__main__":
    main()
