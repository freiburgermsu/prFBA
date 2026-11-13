#!/usr/bin/env python
"""flux.py -- part 3: what E3 changes in prFBA fluxes of models built from the deployed depth-100
genome objects. DESCRIPTIVE: there is no truth here, only the k/n-vs-E3 difference.

Sample (deterministic). Per deployment, the built ASVs whose E3 re-weighting was applied
(n >= 2 source genomes), stratified by n: 2 | 3-4 | 5-9 | >= 10, 10 per stratum, ranked by
sha256("downstream_e3|" + asv); plus 5 n = 1 controls ranked the same way (their genome objects
are byte-identical between the builds, so every result must be identical).

Per genome (one worker process per genome, recycled after it):
  1. draft from the k/n genome object AND from the E3 genome object with
     build_prfba_community_models._build_one (codiffusion's production draft rule:
     MSBuilder.build_base_model(rid, "0", True, False) + ATPM, notes['probability'] = max over GPR
     genes); template class by that module's stage_templates() rule on the genome lineage
     (+ MiDAS kingdom for codiffusion ASVs; EmilyKin has no MiDAS table keyed by ASV md5, so its
     rule is lineage-only). Templates: ModelSEEDTemplates 21a8084 (c12run, the Stage O pin).
     The two drafts are compared reaction by reaction (ids, stoichiometry, bounds, GPR).
  2. ONE gapfill, on the k/n draft, the way Stage O gapfills (phenotype_fold.gapfill_block,
     pipeline as_shipped, A-C12 reversals applied): MSATPCorrection on the ModelSEED ATP test
     media against Core-V6 -> tests -> MSGapfill on the KBase complete medium
     (ModelSEEDDatabase Media/KBaseMedia.cpd, uptake 100) -> integrate + reversals. Gapfilling
     reads no probabilities (no reaction_scores), and the drafts are identical apart from the
     notes, so one gapfill serves both arms.
  3. annotate two copies with prfba_objective.set_reaction_probabilities (after gapfilling:
     gapfilled and other gene-less reactions p = 0.0) from the k/n and from the E3 genome
     object (prfba_objective.feature_probabilities_from_genome, the production loader rule).
  4. prfba_objective.run_prfba on each copy, tuned (ex_weight 1, missing_probability 0,
     biomass bio1) and released (ex_weight 100; objective set to bio1 first), environment =
     the complete medium's compounds at uptake 100 (Stage O's environment). The k/n tuned solve
     is repeated on a fresh copy (solver determinism).

Build set (--builds, required; common.BUILD_SETS): which k/n and E3 builds, and where the
results and models go -- mixed: results/, models/; bvbrc: results_bvbrc/, models_bvbrc/. The
sample rule reads only the E3 manifest (status, n_genomes, pfunc.applied), so a build set whose
manifests agree on those draws the same sample; main() checks it against the other sets' saved
samples.

  run via systemd-run (README); outputs <results>/part3/{deployment}/{asv}.json, <results>/flux.json
"""
from __future__ import annotations

import argparse
import json
import math
import os
import statistics
import sys
import time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

TAG = "v7c12"                       # any tag starting "v7" makes _load_template read TPL_V7
MODELS = f"{C.MODELS_ROOT}/part3"         # per build set (common.BUILD_SETS, --builds)
RES = f"{C.RESULTS}/part3"
STRATA = ["2", "3-4", "5-9", ">=10"]
PER_STRATUM = 10
N_CONTROLS = 5
EPS = 1e-9


def setup_env():
    os.environ.setdefault("PYTHONHASHSEED", "0")
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["CODIF_TAG"] = TAG
    os.environ["CODIF_TPL_ROOT"] = C.TPL_ROOT["c12run"]
    os.environ["PYTHONPATH"] = C.PYTHONPATH
    for p in reversed(C.PYTHONPATH.split(":")):
        if p not in sys.path:
            sys.path.insert(0, p)
    for p in (C.CODIF, C.PRFBA):
        if p not in sys.path:
            sys.path.append(p)


def import_B(dep, gdir):
    setup_env()
    import build_prfba_community_models as B
    B.TAG = TAG
    B.GDIR = gdir
    B.SETS_PATH = f"{RES}/_work/{dep}_sets.json"
    B.TPL_ASSIGN_PATH = f"{RES}/_work/{dep}_template_assignment.json"
    B.TAXONOMY_CSV = (f"{C.CODIF}/new_data/taxonomy_codif_all.csv" if dep == "codiffusion"
                      else f"{RES}/_work/empty_taxonomy.csv")
    assert B.TPL_ROOT == C.TPL_ROOT["c12run"], B.TPL_ROOT
    return B


# ------------------------------------------------------------------ sample
def choose_sample():
    out = {}
    for dep, b in C.BUILDS.items():
        man = json.load(open(f"{b['e3']}/synthetic_genome_manifest.json"))
        strata = defaultdict(list)
        for asv, v in man.items():
            if asv.startswith("_") or v.get("status") != "built":
                continue
            n = v["n_genomes"]
            applied = bool((v.get("pfunc") or {}).get("applied"))
            if n >= 2 and not applied:
                continue
            strata[C.n_stratum(n)].append(asv)
        pick = []
        for s in STRATA:
            pool = sorted(strata[s], key=C.sample_rank)
            pick += [{"asv": a, "n": man[a]["n_genomes"], "stratum": s, "control": False}
                     for a in pool[:PER_STRATUM]]
        pool = sorted(strata["1"], key=C.sample_rank)
        pick += [{"asv": a, "n": 1, "stratum": "1", "control": True} for a in pool[:N_CONTROLS]]
        out[dep] = {"pool_sizes": {s: len(strata[s]) for s in ["1"] + STRATA}, "genomes": pick}
    return out


def sample_agreement(sample):
    """For every other build set with a saved sample: is it the same genomes, n, strata, controls?"""
    key = lambda smp: {d: [(g["asv"], g["n"], g["stratum"], g["control"]) for g in v["genomes"]]
                       for d, v in smp.items()}
    out = {}
    for name, spec in C.BUILD_SETS.items():
        p = f"{spec['results']}/part3/_work/sample.json"
        if name != C.BUILD_SET and os.path.exists(p):
            out[name] = key(json.load(open(p))) == key(sample)
    return out


# ------------------------------------------------------------------ worker
_W = {}


def _init():
    setup_env()
    C.check_imports()
    import logging
    for name in ("modelseedpy", "cobra", "optlang", "mscommunity"):
        logging.getLogger(name).setLevel(logging.ERROR)
    try:
        import gurobipy
        gurobipy.setParam("Threads", 1)
        gurobipy.setParam("OutputFlag", 0)
    except Exception:
        pass
    from modelseedpy.biochem import modelseed_biochem as MB
    if MB.ModelSEEDBiochem.default_biochemistry is None:      # phenotype_arm.preload_biochemistry
        MB.ModelSEEDBiochem.default_biochemistry = MB.from_local(C.MSD)
    from modelseedpy.core.msmedia import MSMedia
    cpds = [l.strip() for l in open(f"{C.MSD}/Media/KBaseMedia.cpd") if l.strip()]
    media = MSMedia.from_dict({c: 100.0 for c in cpds})       # phenotype_arm.complete_media
    media.id = "complete"
    _W["media"], _W["env"] = media, {c: 100.0 for c in cpds}


def _draft(dep, arm, asv):
    B = import_B(dep, C.BUILDS[dep][arm])
    d = f"{MODELS}/{dep}/{arm}"
    os.makedirs(f"{d}/models_{TAG}", exist_ok=True)
    os.chdir(d)
    B._W.clear()
    B._binit()
    rid, res = B._build_one((asv, asv))
    return f"{d}/models_{TAG}/{asv}.xml", res, B


def _sig(r):
    """stoichiometry, bounds and the GPR as a DNF set (GPR string order is not reproducible)."""
    return (tuple(sorted((m.id, round(c, 12)) for m, c in r.metabolites.items())),
            round(r.lower_bound, 9), round(r.upper_bound, 9), C.gpr_canonical(r.gene_reaction_rule))


def _gapfill(B, draft, cls):
    """phenotype_fold.gapfill_block (pipeline as_shipped, integrate_reversals) on a draft model."""
    from modelseedpy import MSGapfill
    from modelseedpy.core.msatpcorrection import MSATPCorrection, load_default_medias
    from modelseedpy.core.msbuilder import MSBuilder
    from modelseedpy.core.msmodel import get_reaction_constraints_from_direction
    from modelseedpy.core.mstemplate import MSTemplateBuilder
    tmpl = B._load_template(cls)
    core_t = MSTemplateBuilder.from_dict(json.load(open(f"{C.TPL_ROOT['c12run']}/Core-V6.json"))).build()
    base = draft
    if "ATPM_c0" not in base.reactions:
        MSBuilder.add_atpm(base)
    atpc = MSATPCorrection(base, core_t, load_default_medias(), compartment="c0",
                           atp_hydrolysis_id="ATPM_c0", load_default_medias=False)
    atpc.evaluate_growth_media()
    atpc.determine_growth_media()
    atpc.apply_growth_media_gapfilling()
    atpc.expand_model_to_genome_scale()
    tests = atpc.build_tests()
    gf = MSGapfill(base, default_gapfill_templates=[tmpl], test_conditions=tests, default_target="bio1")
    infeasible = ""
    try:
        res = gf.run_gapfilling(_W["media"])
    except Exception as exc:  # noqa: BLE001
        res, infeasible = None, f"gapfill_exception:{type(exc).__name__}:{str(exc)[:120]}"
    gap_sol = {}
    if res and res.get("new"):
        for rid, d in res["new"].items():
            if rid[:-1] in tmpl.reactions:
                gap_sol[rid[:-1]] = get_reaction_constraints_from_direction(d)
    elif res is None and not infeasible:
        infeasible = "gapfill_infeasible"
    proposed_rev = dict((res or {}).get("reversed") or {})
    mdl = base.copy()
    MSBuilder.integrate_gapfill_solution(tmpl, mdl, gap_sol)
    applied = []
    for rid, d in sorted(proposed_rev.items()):
        if rid in mdl.reactions:
            r = mdl.reactions.get_by_id(rid)
            if d in (">", "="):
                r.upper_bound = max(r.upper_bound, 1000.0)
            if d in ("<", "="):
                r.lower_bound = min(r.lower_bound, -1000.0)
            applied.append(rid)
    return mdl, {"n_new_rxn": len(gap_sol), "new_rxn": sorted(gap_sol), "n_reversed_applied": len(applied),
                 "infeasible": infeasible, "selected_atp_media": [m.id for m in (atpc.selected_media or [])],
                 "n_tests": len(tests)}


def _max_growth(model):
    from modelseedpy.core.msmodelutl import MSModelUtil
    from modelseedpy.fbapkg.elementuptakepkg import ElementUptakePkg
    m = model.copy()
    u = MSModelUtil(m)
    u.add_medium(_W["env"])
    ElementUptakePkg(m).build_package({"C": 100})
    m.objective = "bio1"
    v = m.slim_optimize()
    return None if v != v else float(v)


def _metrics(model, sol, probs, changed, cfg):
    import prfba_objective as PO
    kw = PO.prfba_kwargs(cfg == "tuned")
    fl = sol.fluxes
    rx = [r for r in model.reactions if r.id.startswith("rxn")]
    ex = [r for r in model.reactions if r.id.startswith("EX_")]
    gb = {r.id for r in rx if r.genes}
    cost = {r.id: max(kw["min_prob"], 1 - probs[r.id] ** kw["prob_exp"]) for r in rx}
    a = {r.id: abs(float(fl[r.id])) for r in rx}
    e = {r.id: float(fl[r.id]) for r in ex}
    tot = sum(a.values())
    gbf = sum(v for k, v in a.items() if k in gb)
    return {
        "status": sol.status, "objective": float(sol.objective_value),
        "growth": float(fl["bio1"]) if "bio1" in fl.index else None,
        "internal_term": sum(cost[k] * v for k, v in a.items()),
        "exchange_term": kw["ex_weight"] * sum(abs(v) for v in e.values()),
        "flux_rxn_total": tot, "flux_gene_backed": gbf, "flux_geneless_rxn": tot - gbf,
        "flux_p_changed": sum(v for k, v in a.items() if k in changed),
        "flux_exchange_total": sum(abs(v) for v in e.values()),
        "uptake_total": sum(-v for v in e.values() if v < 0),
        "secretion_total": sum(v for v in e.values() if v > 0),
        "flux_weighted_mean_p_gene_backed": (sum(probs[k] * v for k, v in a.items() if k in gb) / gbf) if gbf > EPS else None,
        "n_active_rxn": sum(1 for v in a.values() if v > EPS),
        "n_active_p_changed": sum(1 for k, v in a.items() if v > EPS and k in changed),
    }


def _cost(fl, probs, rx, ex, kw):
    return (sum(max(kw["min_prob"], 1 - probs[r] ** kw["prob_exp"]) * abs(float(fl[r])) for r in rx)
            + kw["ex_weight"] * sum(abs(float(fl[r])) for r in ex))


def _pair(model, s_kn, s_e3, probs, cfg):
    import prfba_objective as PO
    kw = PO.prfba_kwargs(cfg == "tuned")
    rx = [r.id for r in model.reactions if r.id.startswith("rxn")]
    ex = [r.id for r in model.reactions if r.id.startswith("EX_")]
    fk, fe = s_kn.fluxes, s_e3.fluxes
    # each arm's flux priced under both probability sets: how far is the E3 optimum from the
    # k/n optimum under k/n prices, and vice versa (>= 0 up to solver tolerance)
    cross = {"kn_price_of_kn_flux": _cost(fk, probs["kn"], rx, ex, kw),
             "kn_price_of_e3_flux": _cost(fe, probs["kn"], rx, ex, kw),
             "e3_price_of_e3_flux": _cost(fe, probs["e3"], rx, ex, kw),
             "e3_price_of_kn_flux": _cost(fk, probs["e3"], rx, ex, kw)}
    cross["e3_flux_excess_under_kn_prices"] = cross["kn_price_of_e3_flux"] - cross["kn_price_of_kn_flux"]
    cross["kn_flux_excess_under_e3_prices"] = cross["e3_price_of_kn_flux"] - cross["e3_price_of_e3_flux"]
    cross["objective_change_from_repricing_only"] = cross["e3_price_of_kn_flux"] - cross["kn_price_of_kn_flux"]
    l1 = sum(abs(float(fe[r]) - float(fk[r])) for r in rx)
    tk = sum(abs(float(fk[r])) for r in rx)
    te = sum(abs(float(fe[r])) for r in rx)
    ak = {r for r in rx if abs(float(fk[r])) > EPS}
    ae = {r for r in rx if abs(float(fe[r])) > EPS}
    dex = sorted(((r, float(fe[r]) - float(fk[r])) for r in ex if abs(float(fe[r]) - float(fk[r])) > 1e-6),
                 key=lambda x: -abs(x[1]))
    return {**cross, "l1_rxn": l1, "rerouted_share": (l1 / (tk + te)) if tk + te > EPS else 0.0,
            "max_abs_dflux": max((abs(float(fe[r]) - float(fk[r])) for r in rx + ex), default=0.0),
            "jaccard_active": (len(ak & ae) / len(ak | ae)) if ak | ae else 1.0,
            "n_active_only_kn": len(ak - ae), "n_active_only_e3": len(ae - ak),
            "l1_exchange": sum(abs(d) for _, d in dex), "n_exchange_changed": len(dex),
            "top_exchange_changes": [[r, round(d, 6), round(float(fk[r]), 6), round(float(fe[r]), 6)]
                                     for r, d in dex[:8]]}


def run_one(job):
    _init_once()
    dep, asv = job["dep"], job["asv"]
    out_path = f"{RES}/{dep}/{asv}.json"
    t0 = time.time()
    rec = {"dep": dep, "asv": asv, "n": job["n"], "stratum": job["stratum"], "control": job["control"],
           "timing": {}}
    try:
        from cobra.io import read_sbml_model
        import prfba_objective as PO
        paths, B = {}, None
        for arm in ("kn", "e3"):
            t = time.time()
            paths[arm], res, B = _draft(dep, arm, asv)
            rec[f"draft_{arm}"] = res
            rec["timing"][f"draft_{arm}"] = round(time.time() - t, 1)
            if res.startswith("ERROR"):
                raise RuntimeError(f"draft {arm}: {res}")
        rec["template"] = B._W["assign"][asv]
        dk, de = read_sbml_model(paths["kn"]), read_sbml_model(paths["e3"])
        rk, re_ = {r.id: r for r in dk.reactions}, {r.id: r for r in de.reactions}
        notes_p = {i: float(r.notes["probability"]) for i, r in rk.items()}
        notes_p_e3 = {i: float(r.notes["probability"]) for i, r in re_.items()}
        rec["drafts"] = {"n_rxn": len(rk), "same_reaction_ids": set(rk) == set(re_),
                         "n_struct_diff": sum(1 for i in rk if i in re_ and _sig(rk[i]) != _sig(re_[i])),
                         "n_gpr_string_order_only": sum(1 for i in rk if i in re_ and rk[i].gene_reaction_rule
                                                        != re_[i].gene_reaction_rule and _sig(rk[i]) == _sig(re_[i])),
                         "n_rxn_p_changed": sum(1 for i in rk if i in re_ and abs(notes_p[i] - notes_p_e3[i]) > EPS),
                         "n_gene_bearing": sum(1 for r in dk.reactions if r.genes)}
        if not rec["drafts"]["same_reaction_ids"] or rec["drafts"]["n_struct_diff"]:
            raise RuntimeError("k/n and E3 drafts differ in structure")
        del de
        t = time.time()
        dk.id = asv
        gfm, gfrec = _gapfill(B, dk, rec["template"])
        rec["gapfill"] = gfrec
        rec["timing"]["gapfill"] = round(time.time() - t, 1)
        fp = {arm: PO.feature_probabilities_from_genome(f"{C.BUILDS[dep][arm]}/genome_objects/{asv}.json")
              for arm in ("kn", "e3")}
        annotated, probs = {}, {}
        for arm in ("kn", "e3"):
            m = gfm.copy()
            probs[arm] = PO.set_reaction_probabilities(m, fp[arm])
            annotated[arm] = m
        # the production notes of the drafts == set_reaction_probabilities on the gapfilled model?
        rec["notes_vs_rule_mismatch"] = sum(1 for i, p in notes_p.items()
                                            if i in probs["kn"] and abs(p - probs["kn"][i]) > EPS)
        rx_ids = [r.id for r in gfm.reactions if r.id.startswith("rxn")]
        changed = {i for i in rx_ids if abs(probs["kn"][i] - probs["e3"][i]) > EPS}
        dp = [probs["e3"][i] - probs["kn"][i] for i in changed]
        rec["model"] = {"n_rxn": len(gfm.reactions), "n_rxn_prefixed": len(rx_ids),
                        "n_gene_backed": sum(1 for r in gfm.reactions if r.id.startswith("rxn") and r.genes),
                        "n_p_changed": len(changed),
                        "n_p_up": sum(1 for x in dp if x > 0), "n_p_down": sum(1 for x in dp if x < 0),
                        "mean_dp_changed": statistics.mean(dp) if dp else 0.0,
                        "n_p_zero_e3_gene_backed": sum(1 for r in gfm.reactions if r.id.startswith("rxn")
                                                       and r.genes and probs["e3"][r.id] == 0)}
        rec["max_growth_complete_c100"] = _max_growth(gfm)
        sols, rec["prfba"], rec["pairs"] = {}, {}, {}
        t = time.time()
        for cfg in ("tuned", "released"):
            for arm in ("kn", "e3"):
                M = annotated[arm].copy()
                if cfg == "released":
                    M.objective = "bio1"
                try:
                    s = PO.run_prfba(M, environment=_W["env"], tuned=(cfg == "tuned"))
                    sols[(cfg, arm)] = s
                    rec["prfba"][f"{cfg}/{arm}"] = _metrics(annotated[arm], s, probs[arm], changed, cfg)
                except Exception as exc:  # noqa: BLE001
                    rec["prfba"][f"{cfg}/{arm}"] = {"error": f"{type(exc).__name__}: {str(exc)[:300]}"}
            if (cfg, "kn") in sols and (cfg, "e3") in sols:
                rec["pairs"][cfg] = _pair(gfm, sols[(cfg, "kn")], sols[(cfg, "e3")], probs, cfg)
        # determinism: the k/n tuned solve again on a fresh copy
        if ("tuned", "kn") in sols:
            M = annotated["kn"].copy()
            s2 = PO.run_prfba(M, environment=_W["env"], tuned=True)
            f1, f2 = sols[("tuned", "kn")].fluxes, s2.fluxes
            rec["determinism_tuned_kn"] = {
                "objective_diff": float(s2.objective_value) - float(sols[("tuned", "kn")].objective_value),
                "max_abs_flux_diff": float(max(abs(float(f1[r]) - float(f2[r])) for r in f1.index))}
        rec["timing"]["prfba"] = round(time.time() - t, 1)
        rec["fluxes"] = {f"{cfg}/{arm}": {r: round(float(v), 9) for r, v in s.fluxes.items() if abs(v) > EPS}
                         for (cfg, arm), s in sols.items()}
        rec["p_changed_reactions"] = {i: [round(probs["kn"][i], 6), round(probs["e3"][i], 6)]
                                      for i in sorted(changed)
                                      if any(abs(float(s.fluxes[i])) > EPS for s in sols.values())}
        rec["status"] = "ok"
        try:
            from cobra.io import write_sbml_model
            os.makedirs(f"{MODELS}/{dep}/gapfilled", exist_ok=True)
            write_sbml_model(annotated["kn"], f"{MODELS}/{dep}/gapfilled/{asv}.kn.xml")
        except Exception as exc:  # noqa: BLE001
            rec["write_warning"] = str(exc)[:200]
    except Exception as exc:  # noqa: BLE001
        import traceback
        rec["status"] = "error"
        rec["error"] = f"{type(exc).__name__}: {str(exc)[:400]}"
        rec["traceback"] = traceback.format_exc()[-2000:]
    rec["timing"]["total"] = round(time.time() - t0, 1)
    try:
        import resource
        rec["maxrss_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576, 2)
    except Exception:
        pass
    C.dump(rec, out_path)
    return {k: rec.get(k) for k in ("dep", "asv", "status", "error")} | {"secs": rec["timing"]["total"]}


_INIT = {"done": False}


def _init_once():
    if not _INIT["done"]:
        _init()
        _INIT["done"] = True


# ------------------------------------------------------------------ summary
def q(xs, f):
    s = sorted(xs)
    return s[min(len(s) - 1, int(round(f * (len(s) - 1))))] if s else None


def describe(xs):
    xs = [x for x in xs if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if not xs:
        return {"n": 0}
    return {"n": len(xs), "mean": statistics.mean(xs), "median": q(xs, 0.5), "p10": q(xs, 0.1),
            "p25": q(xs, 0.25), "p75": q(xs, 0.75), "p90": q(xs, 0.9), "min": min(xs), "max": max(xs)}


def summarise(sample):
    rows = []
    for dep, v in sample.items():
        for g in v["genomes"]:
            p = f"{RES}/{dep}/{g['asv']}.json"
            if os.path.exists(p):
                r = json.load(open(p))
                r.pop("fluxes", None)
                rows.append(r)
    ok = [r for r in rows if r["status"] == "ok"]
    out = {"n_genomes": len(rows), "n_ok": len(ok),
           "errors": [{k: r.get(k) for k in ("dep", "asv", "n", "error")} for r in rows if r["status"] != "ok"]}
    metrics = ["objective", "growth", "flux_rxn_total", "flux_gene_backed", "flux_p_changed",
               "flux_exchange_total", "uptake_total", "flux_weighted_mean_p_gene_backed",
               "internal_term", "exchange_term", "n_active_rxn"]
    for group, sel in (("n>=2", lambda r: not r["control"]), ("n=1 controls", lambda r: r["control"])):
        for dep in list(C.BUILDS) + ["pooled"]:
            R = [r for r in ok if sel(r) and (dep == "pooled" or r["dep"] == dep)]
            if not R:
                continue
            blk = {"n_genomes": len(R)}
            for cfg in ("tuned", "released"):
                P = [r for r in R if f"{cfg}/kn" in r["prfba"] and "error" not in r["prfba"][f"{cfg}/kn"]
                     and "error" not in r["prfba"].get(f"{cfg}/e3", {"error": 1})]
                c = {"n_solved_both": len(P),
                     "n_errors": sum(1 for r in R if r not in P)}
                for m in metrics:
                    kn = [r["prfba"][f"{cfg}/kn"][m] for r in P]
                    e3 = [r["prfba"][f"{cfg}/e3"][m] for r in P]
                    d = [b - a if (a is not None and b is not None) else None for a, b in zip(kn, e3)]
                    rel = [(b - a) / abs(a) if (a not in (None, 0) and b is not None) else None for a, b in zip(kn, e3)]
                    c[m] = {"kn": describe(kn), "e3": describe(e3), "e3_minus_kn": describe(d),
                            "relative": describe(rel),
                            "sum_kn": sum(x for x in kn if x is not None), "sum_e3": sum(x for x in e3 if x is not None),
                            "n_e3_lower": sum(1 for x in d if x is not None and x < -1e-9),
                            "n_e3_higher": sum(1 for x in d if x is not None and x > 1e-9)}
                pr = [r["pairs"][cfg] for r in P if cfg in r.get("pairs", {})]
                c["rerouted_share"] = describe([x["rerouted_share"] for x in pr])
                for k in ("e3_flux_excess_under_kn_prices", "kn_flux_excess_under_e3_prices",
                          "objective_change_from_repricing_only"):
                    c[k] = describe([x[k] for x in pr])
                c["relative_e3_flux_excess_under_kn_prices"] = describe(
                    [x["e3_flux_excess_under_kn_prices"] / x["kn_price_of_kn_flux"] for x in pr if x["kn_price_of_kn_flux"] > EPS])
                c["relative_kn_flux_excess_under_e3_prices"] = describe(
                    [x["kn_flux_excess_under_e3_prices"] / x["e3_price_of_e3_flux"] for x in pr if x["e3_price_of_e3_flux"] > EPS])
                c["max_rel_objective_decomposition_error"] = max(
                    (abs(r["prfba"][f"{cfg}/{arm}"]["objective"] - r["prfba"][f"{cfg}/{arm}"]["internal_term"]
                         - r["prfba"][f"{cfg}/{arm}"]["exchange_term"]) / max(1.0, abs(r["prfba"][f"{cfg}/{arm}"]["objective"]))
                     for r in P for arm in ("kn", "e3")), default=None)
                c["jaccard_active"] = describe([x["jaccard_active"] for x in pr])
                c["l1_exchange"] = describe([x["l1_exchange"] for x in pr])
                c["n_genomes_flux_identical"] = sum(1 for x in pr if x["max_abs_dflux"] <= 1e-9)
                c["n_genomes_any_exchange_change"] = sum(1 for x in pr if x["n_exchange_changed"])
                pk = [r["prfba"][f"{cfg}/kn"] for r in P]
                pe = [r["prfba"][f"{cfg}/e3"] for r in P]
                c["share_of_flux_through_p_changed"] = {
                    "kn": describe([a["flux_p_changed"] / a["flux_rxn_total"] for a in pk if a["flux_rxn_total"] > EPS]),
                    "e3": describe([a["flux_p_changed"] / a["flux_rxn_total"] for a in pe if a["flux_rxn_total"] > EPS])}
                c["pooled_share_rerouted"] = (sum(x["l1_rxn"] for x in pr) /
                                              sum(a["flux_rxn_total"] + b["flux_rxn_total"] for a, b in zip(pk, pe))
                                              if pr else None)
                blk[cfg] = c
            blk["model"] = {k: describe([r["model"][k] for r in R]) for k in
                            ("n_rxn_prefixed", "n_gene_backed", "n_p_changed", "n_p_up", "n_p_down", "mean_dp_changed",
                             "n_p_zero_e3_gene_backed")}
            blk["gapfill_n_new_rxn"] = describe([r["gapfill"]["n_new_rxn"] for r in R])
            # draft reactions absent from the gapfilled model (removed by ATP correction, not re-added)
            blk["draft_rxn_net_removed"] = describe([r["drafts"]["n_rxn"] + r["gapfill"]["n_new_rxn"] - r["model"]["n_rxn"]
                                                     for r in R])
            blk["determinism_max_abs_flux_diff"] = max((r.get("determinism_tuned_kn", {}).get("max_abs_flux_diff", 0)
                                                       for r in R), default=None)
            out.setdefault(group, {})[dep] = blk
    # n >= 2, both deployments, by n stratum
    out["by_stratum"] = {}
    for st in STRATA:
        R = [r for r in ok if r["stratum"] == st]
        blk = {"n_genomes": len(R)}
        for cfg in ("tuned", "released"):
            pr = [r["pairs"][cfg] for r in R if cfg in r.get("pairs", {})]
            pk = [r["prfba"][f"{cfg}/kn"] for r in R if cfg in r.get("pairs", {})]
            pe = [r["prfba"][f"{cfg}/e3"] for r in R if cfg in r.get("pairs", {})]
            blk[cfg] = {"rerouted_share": describe([x["rerouted_share"] for x in pr]),
                        "pooled_share_rerouted": (sum(x["l1_rxn"] for x in pr) /
                                                  sum(a["flux_rxn_total"] + b["flux_rxn_total"] for a, b in zip(pk, pe))
                                                  if pr else None),
                        "n_exchange_changed_genomes": sum(1 for x in pr if x["n_exchange_changed"]),
                        "objective_relative": describe([(b["objective"] - a["objective"]) / a["objective"]
                                                        for a, b in zip(pk, pe)]),
                        "relative_e3_flux_excess_under_kn_prices": describe(
                            [x["e3_flux_excess_under_kn_prices"] / x["kn_price_of_kn_flux"] for x in pr])}
        blk["n_p_changed"] = describe([r["model"]["n_p_changed"] for r in R])
        blk["mean_dp_changed"] = describe([r["model"]["mean_dp_changed"] for r in R])
        out["by_stratum"][st] = blk
    # per-genome table
    out["per_genome"] = []
    for r in sorted(rows, key=lambda r: (r["dep"], r["control"], STRATA.index(r["stratum"]) if r["stratum"] in STRATA else -1, r["asv"])):
        row = {k: r.get(k) for k in ("dep", "asv", "n", "stratum", "control", "status", "template")}
        if r["status"] == "ok":
            row["n_p_changed"] = r["model"]["n_p_changed"]
            row["gapfill_n_new"] = r["gapfill"]["n_new_rxn"]
            for cfg in ("tuned", "released"):
                a, b = r["prfba"].get(f"{cfg}/kn", {}), r["prfba"].get(f"{cfg}/e3", {})
                if "error" in a or "error" in b:
                    row[cfg] = {"error": a.get("error") or b.get("error")}
                    continue
                pr = r["pairs"][cfg]
                row[cfg] = {"objective": [a["objective"], b["objective"]], "growth": [a["growth"], b["growth"]],
                            "flux_gene_backed": [a["flux_gene_backed"], b["flux_gene_backed"]],
                            "flux_p_changed": [a["flux_p_changed"], b["flux_p_changed"]],
                            "flux_exchange_total": [a["flux_exchange_total"], b["flux_exchange_total"]],
                            "rerouted_share": pr["rerouted_share"], "jaccard_active": pr["jaccard_active"],
                            "max_abs_dflux": pr["max_abs_dflux"]}
        else:
            row["error"] = r.get("error")
        out["per_genome"].append(row)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--builds", required=True, choices=sorted(C.BUILD_SETS),
                    help="build set (common.BUILD_SETS): its k/n and E3 builds, results and models dirs")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--limit", type=int, default=None, help="smoke: first N jobs only")
    ap.add_argument("--summary-only", action="store_true",
                    help="re-summarise the saved per-genome files (uses the saved sample)")
    a = ap.parse_args()
    assert a.builds == C.BUILD_SET, (a.builds, C.BUILD_SET)
    setup_env()
    C.ensure_pins()
    os.makedirs(f"{RES}/_work", exist_ok=True)
    t_all = time.time()
    if a.summary_only and os.path.exists(f"{RES}/_work/sample.json"):
        sample = json.load(open(f"{RES}/_work/sample.json"))
    else:
        C.require_builds()
        sample = choose_sample()
        C.dump(sample, f"{RES}/_work/sample.json")
    same = sample_agreement(sample)
    C.log(f"sample vs the other build sets' saved samples: {same}")
    jobs = [{"dep": dep, **g} for dep, v in sample.items() for g in v["genomes"]]
    if not a.summary_only:
        # template classes: stage_templates() on the sampled k/n genome objects (own output files)
        with open(f"{RES}/_work/empty_taxonomy.csv", "w") as fh:
            fh.write("seq,Kingdom,Phylum\n")
        for dep in C.BUILDS:
            if os.path.exists(f"{RES}/_work/{dep}_template_assignment.json"):
                continue
            B = import_B(dep, C.BUILDS[dep]["kn"])
            sets = {g["asv"]: {"rep_asv": g["asv"]} for g in sample[dep]["genomes"]}
            C.dump({"_meta": {"note": "downstream_e3 part 3 sample"}, **sets}, B.SETS_PATH)
            B.stage_templates()
        todo = [j for j in jobs if not (os.path.exists(f"{RES}/{j['dep']}/{j['asv']}.json") and
                                        json.load(open(f"{RES}/{j['dep']}/{j['asv']}.json")).get("status") == "ok")]
        if a.limit is not None:
            todo = todo[:a.limit]
        C.log(f"part 3: {len(jobs)} genomes in the sample, {len(todo)} to run, {a.workers} workers")
        import multiprocessing as mp
        # one task per process (spawn): cobra copies leak solver objects; templates hold gapfill
        # state. multiprocessing.Pool: ProcessPoolExecutor(max_tasks_per_child) deadlocked under
        # Python 3.12 once workers retired (penalties.py, 2026-11-12)
        with mp.get_context("spawn").Pool(a.workers, maxtasksperchild=1) as pool:
            for i, r in enumerate(pool.imap_unordered(run_one, todo, chunksize=1), 1):
                C.log(f"  [{i}/{len(todo)}] {r['dep']} {r['asv']} {r['status']} {r.get('error') or ''} ({r['secs']} s)")
    summ = summarise(sample)
    prev_secs = None
    if a.summary_only:             # keep the run's wall time: from the last flux.json, else its log
        if os.path.exists(f"{C.RESULTS}/flux.json"):
            prev_secs = json.load(open(f"{C.RESULTS}/flux.json"))["_meta"].get("secs_total_run")
        log_path = f"{C.OUT}/logs/flux.log" if C.BUILD_SET == "mixed" else f"{C.OUT}/logs/flux_{C.BUILD_SET}.log"
        if prev_secs is None and os.path.exists(log_path):
            import re
            m = re.findall(r"part 3 done in ([0-9.]+) s", open(log_path).read())
            prev_secs = float(m[-1]) if m else None
    res = {"_meta": {"script": "flux.py", "code": C.code_record([os.path.abspath(__file__)]),
                     "pins": C.pin_record(), "templates_root": C.TPL_ROOT["c12run"],
                     "sample_rule": f"per deployment: E3-applied built ASVs, strata {STRATA}, {PER_STRATUM} per "
                                    f"stratum by sha256('downstream_e3|'+asv); {N_CONTROLS} n=1 controls likewise",
                     "pool_sizes": {d: v["pool_sizes"] for d, v in sample.items()},
                     "sample_equals_saved_sample_of": same,
                     "environment": "KBase complete medium (ModelSEEDDatabase Media/KBaseMedia.cpd), every "
                                    "compound at uptake 100; prFBA's carbon cap C <= 100",
                     "configs": {"tuned": "prfba_objective.TUNED_PRFBA_KWARGS", "released": "RELEASED_PRFBA_KWARGS "
                                 "(objective bio1 set first)"},
                     "definitions": {
                         "flux_gene_backed": "sum |v| over rxn reactions with a GPR in the gapfilled model",
                         "flux_p_changed": "sum |v| over rxn reactions whose p differs between k/n and E3",
                         "rerouted_share": "sum_rxn |v_e3 - v_kn| / (sum_rxn |v_kn| + sum_rxn |v_e3|): 0 identical, "
                                           "1 disjoint supports",
                         "pooled_share_rerouted": "the same ratio with numerator and denominator summed over genomes",
                         "objective": "value of prFBA's weighted minimisation (each arm priced with its own p)"},
                     "secs_total_run": prev_secs if a.summary_only else round(time.time() - t_all, 1),
                     "summary_written": time.strftime("%Y-%m-%dT%H:%M:%S")},
           **summ}
    C.dump(res, f"{C.RESULTS}/flux.json" if a.limit is None else f"{RES}/_work/flux_smoke.json")
    C.log(f"part 3 done ({res['_meta']['secs_total_run']} s run)")


if __name__ == "__main__":
    main()
