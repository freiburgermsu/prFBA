#!/usr/bin/env python
"""penalties.py -- part 2: what E3 changes in compartment-aware community gapfill penalties
(the "2 - p" factor of codiffusion_bioreactor/compartment_aware_gapfill.py, a PROTOTYPE whose
only caller is demo_compartment_scoring.py).

The demo (codiffusion_bioreactor/demo_compartment_scoring.py, output
modeling_files/compartment_scoring_demo_day174_4members.json) prices every gapfill candidate of a
day-174 subcommunity twice: stock (GapfillingPkg strips _cN; no reaction_scores -> factor 2 for
all) and compartment-aware (candidate rxn_ci priced base x (2 - max p of member i's genes for
rxn)). Here the member drafts are rebuilt from the DEPTH-100 genome objects, once from the
shipped k/n build and once from the E3 build, with build_prfba_community_models._build_one itself
(imported, not edited; module globals pointed at our paths), and the penalty vectors compared.

Member -> genome object. A codiffusion member is an ASVset of the 20-hit build
(modeling_files/ASV_sets_54d1150.json: built ASVs grouped by identical selected-genome lists,
rep = highest summed abundance). Its d100 genome object is the file named by the SAME
representative ASV md5: synthetic_genomes_d100[_e3]/genome_objects/{rep_asv}.json. Unmapped =
the rep has no built d100 genome. The d100 selection of a rep can differ from the 20-hit one, so a
d100 draft is not the production (54d1150) draft; the comparison is d100-k/n vs d100-E3 only.

Sets
  A  the demo's own 4 members (day 174, top 4 by abundance; Wolfe medium + H2/CO2 feed),
     community assembled with mscommunity.build_from_species_models, gram_neg gapfill template
     for the community (as the demo)
  B  day 174's 4 most abundant members whose d100 rep was reweighted (n >= 2), same procedure
  C  every day-174 member alone (single-member GapfillingPkg, member scores at index 0, same
     gram_neg template): the per-member picture for the whole 88-member community. Set A/B
     members are checked against their community compartments (equivalence block).
For A and B also: does community gapfilling solve at all (stock, compartment-aware k/n,
compartment-aware E3)? One MSGapfill.run_gapfilling each, 30-min solver limit.

Templates: codiffusion's production set ModelSEEDTemplates b6fe86b (codif-curated,
rebuild_pipeline.sh), from the pinned export (common.py). Template class per member:
build_prfba_community_models.stage_templates() run on the d100 k/n genome objects (MiDAS kingdom
+ genome lineage), written to our own file.

  run via systemd-run (see README); outputs: results/penalties.json, models/part2/{kn,e3}/
"""
from __future__ import annotations

import argparse
import json
import os
import re
import statistics
import sys
import time
from collections import defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import common as C  # noqa: E402

DAY = "174"
TAG = "v7e"
MODELS = f"{C.OUT}/models/part2"
WORK = f"{C.OUT}/results/part2_work"
SETS_PATH = f"{C.CODIF}/modeling_files/ASV_sets_54d1150.json"
AB_PATH = f"{C.CODIF}/modeling_files/ASVset_day_abundances_54d1150.json"
DEMO_OUT = f"{C.CODIF}/modeling_files/compartment_scoring_demo_day174_4members.json"
DEMO_PROV = f"{C.CODIF}/comm_models_v7/day174_Jb_comm.provenance.json"
PROD_ASSIGN = f"{C.CODIF}/modeling_files/template_assignment_v7.json"
MY_SETS = f"{WORK}/ASV_sets_day{DAY}.json"
MY_ASSIGN = f"{WORK}/template_assignment_d100.json"
TOL = 1e-9


def setup_env():
    os.environ.setdefault("PYTHONHASHSEED", "0")
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    os.environ["CODIF_TAG"] = TAG
    os.environ["CODIF_TPL_ROOT"] = C.TPL_ROOT["codif_curated"]
    for p in reversed(C.PYTHONPATH.split(":")):
        if p not in sys.path:
            sys.path.insert(0, p)
    if C.CODIF not in sys.path:
        sys.path.insert(1, C.CODIF)
    os.environ["PYTHONPATH"] = C.PYTHONPATH


def import_B(gdir=None):
    """build_prfba_community_models with its relative-path globals pointed at absolute paths."""
    setup_env()
    import build_prfba_community_models as B
    B.TAG = TAG
    B.SETS_PATH = MY_SETS
    B.TPL_ASSIGN_PATH = MY_ASSIGN
    B.TAXONOMY_CSV = f"{C.CODIF}/new_data/taxonomy_codif_all.csv"
    B.MEDIA_JSON = f"{C.CODIF}/modeling_files/kbase_media_207617.json"
    B.AB_PATH = AB_PATH
    if gdir:
        B.GDIR = gdir
    assert B.TPL_ROOT == C.TPL_ROOT["codif_curated"], B.TPL_ROOT
    return B


# ------------------------------------------------------------------ draft building (workers)
_ARM = {}


def _draft_init(arm):
    B = import_B(C.BUILDS["codiffusion"][arm])
    C.check_imports()
    import logging
    for name in ("modelseedpy", "cobra", "optlang"):
        logging.getLogger(name).setLevel(logging.ERROR)
    os.makedirs(f"{MODELS}/{arm}", exist_ok=True)
    os.chdir(f"{MODELS}/{arm}")          # _build_one writes models_<TAG>/{rid}.xml under cwd
    os.makedirs(f"models_{TAG}", exist_ok=True)
    B._binit()
    _ARM["B"], _ARM["arm"] = B, arm


def _draft_one(task):
    B = _ARM["B"]
    t0 = time.time()
    rid, res = B._build_one(task)
    return {"arm": _ARM["arm"], "rid": rid, "result": res, "secs": round(time.time() - t0, 1)}


def build_drafts(tasks, workers):
    import multiprocessing as mp
    out = []
    for arm in ("kn", "e3"):
        todo = [t for t in tasks if not os.path.exists(f"{MODELS}/{arm}/models_{TAG}/{t[0]}.xml")]
        C.log(f"drafts {arm}: {len(todo)} of {len(tasks)} to build")
        if not todo:
            continue
        # multiprocessing.Pool, not ProcessPoolExecutor(max_tasks_per_child=...): the latter
        # deadlocked here (Python 3.12) once workers began to retire
        with mp.get_context("spawn").Pool(workers, initializer=_draft_init, initargs=(arm,),
                                          maxtasksperchild=10) as pool:
            for r in pool.imap_unordered(_draft_one, todo, chunksize=1):
                out.append(r)
                if r["result"].startswith("ERROR"):
                    C.log(f"  draft ERROR {arm} {r['rid']}: {r['result']}")
    return out


# ------------------------------------------------------------------ draft comparison
def _rxn_sig(r):
    """stoichiometry, bounds and the GPR as a DNF set (GPR string order is not reproducible)."""
    return (tuple(sorted((m.id, round(c, 12)) for m, c in r.metabolites.items())),
            round(r.lower_bound, 9), round(r.upper_bound, 9), C.gpr_canonical(r.gene_reaction_rule))


def compare_drafts(rid):
    setup_env()
    from cobra.io import read_sbml_model
    a = read_sbml_model(f"{MODELS}/kn/models_{TAG}/{rid}.xml")
    b = read_sbml_model(f"{MODELS}/e3/models_{TAG}/{rid}.xml")
    ra, rb = {r.id: r for r in a.reactions}, {r.id: r for r in b.reactions}
    same_ids = set(ra) == set(rb)
    struct_diff = [i for i in ra if i in rb and _rxn_sig(ra[i]) != _rxn_sig(rb[i])]
    gpr_string_order_only = sum(1 for i in ra if i in rb and ra[i].gene_reaction_rule != rb[i].gene_reaction_rule
                                and i not in struct_diff)
    pa = {i: float(r.notes.get("probability", "nan")) for i, r in ra.items()}
    pb = {i: float(r.notes.get("probability", "nan")) for i, r in rb.items()}
    changed = {i: (pa[i], pb[i]) for i in pa if i in pb and abs(pa[i] - pb[i]) > TOL}
    notes_other = [i for i in ra if i in rb and
                   {k: v for k, v in ra[i].notes.items() if k != "probability"} !=
                   {k: v for k, v in rb[i].notes.items() if k != "probability"}]
    gb = [i for i, r in ra.items() if r.genes]
    return {"rid": rid, "n_rxn": len(ra), "n_gene_bearing": len(gb), "same_reaction_ids": same_ids,
            "n_struct_diff": len(struct_diff), "struct_diff_examples": struct_diff[:5],
            "n_gpr_string_order_only": gpr_string_order_only,
            "n_other_notes_diff": len(notes_other),
            "n_rxn_p_changed": len(changed),
            "mean_p_gene_bearing_kn": statistics.mean(pa[i] for i in gb) if gb else None,
            "mean_p_gene_bearing_e3": statistics.mean(pb[i] for i in gb) if gb else None,
            "rxn_p_changed": {i: [round(x, 6), round(y, 6)] for i, (x, y) in sorted(changed.items())}}


# ------------------------------------------------------------------ penalties
def scores_for(rids, arm, single_index=None):
    """compartment_aware_gapfill.member_scores_from_drafts over OUR draft dir."""
    setup_env()
    import compartment_aware_gapfill as CAG
    sc = CAG.member_scores_from_drafts({"abundances": {r: 1.0 for r in rids}}, tag=TAG,
                                       repo=f"{MODELS}/{arm}")
    if single_index is not None:
        assert list(sc) == ["1"], list(sc)
        sc = {str(single_index): sc["1"]}
    return sc


def penalty_vectors(model, template_factory, media, scores_by_arm):
    """stock + compartment-aware per arm, with demo_compartment_scoring.penalties_for itself."""
    import demo_compartment_scoring as D
    out, stats = {}, {}
    t0 = time.time()
    out["stock"], stats["stock"] = D.penalties_for(model.copy(), template_factory(), media, {}, patched=False)
    for arm, sc in scores_by_arm.items():
        out[arm], stats[arm] = D.penalties_for(model.copy(), template_factory(), media, sc, patched=True)
    return out, stats, round(time.time() - t0, 1)


IDX = re.compile(r"_([a-z])(\d+)$")


def base_of(rid):
    return re.sub(r"_[a-z]\d+$", "", rid)


def analyse(pen, scores_by_arm, index_to_member, p_lookup):
    """kn vs e3 compartment-aware penalties. p_lookup[arm][index][base] -> member max p."""
    stock, kn, e3 = pen["stock"], pen["kn"], pen["e3"]
    cand = sorted(set(stock) & set(kn) & set(e3))
    rows, per_member = [], defaultdict(lambda: {"n_candidates": 0, "n_member_scored": 0, "n_changed": 0,
                                                "dfactor": [], "rel_change": []})
    for r in cand:
        m = IDX.search(r)
        idx = m.group(2) if m else None
        b = base_of(r)
        mem = index_to_member.get(idx, f"index{idx}")
        pm = per_member[mem]
        pm["n_candidates"] += 1
        scored = idx in scores_by_arm["kn"] and b in scores_by_arm["kn"][idx]
        pm["n_member_scored"] += scored
        base = stock[r] / 2.0
        f_kn, f_e3 = kn[r] / base, e3[r] / base
        if abs(e3[r] - kn[r]) > TOL:
            pm["n_changed"] += 1
            pm["dfactor"].append(f_e3 - f_kn)
            pm["rel_change"].append(e3[r] / kn[r] - 1)
            rows.append({"candidate": r, "member": mem, "base_penalty": round(base, 6),
                         "p_kn": round(p_lookup["kn"].get(idx, {}).get(b, float("nan")), 6),
                         "p_e3": round(p_lookup["e3"].get(idx, {}).get(b, float("nan")), 6),
                         "penalty_kn": round(kn[r], 6), "penalty_e3": round(e3[r], 6),
                         "factor_kn": round(f_kn, 6), "factor_e3": round(f_e3, 6)})
    # stock factor must be exactly 2 everywhere; compartment-aware factor = 2 - p
    stock_bad = 0
    # cheapest-member flips among bases with candidates in >= 2 member compartments
    by_base = defaultdict(dict)
    for r in cand:
        m = IDX.search(r)
        if m and m.group(2) in index_to_member:
            by_base[base_of(r)][m.group(2)] = r
    flips, multi = [], 0
    for b, per in by_base.items():
        if len(per) < 2:
            continue
        multi += 1

        def argmin(vec):
            lo = min(vec[r] for r in per.values())
            return sorted(i for i, r in per.items() if abs(vec[r] - lo) <= TOL)
        a_kn, a_e3, a_st = argmin(kn), argmin(e3), argmin(stock)
        if a_kn != a_e3:
            flips.append({"reaction": b,
                          "cheapest_kn": [index_to_member[i] for i in a_kn],
                          "cheapest_e3": [index_to_member[i] for i in a_e3],
                          "cheapest_stock": [index_to_member[i] for i in a_st],
                          "penalties_kn": {index_to_member[i]: round(kn[r], 4) for i, r in per.items()},
                          "penalties_e3": {index_to_member[i]: round(e3[r], 4) for i, r in per.items()}})
    pm_out = {}
    for mem, v in per_member.items():
        d = v.pop("dfactor")
        rc = v.pop("rel_change")
        v["dfactor_e3_minus_kn"] = describe(d)
        v["relative_penalty_change"] = describe(rc)
        v["n_penalty_up"] = sum(1 for x in d if x > 0)
        v["n_penalty_down"] = sum(1 for x in d if x < 0)
        pm_out[mem] = v
    rows.sort(key=lambda x: -abs(x["penalty_e3"] - x["penalty_kn"]))
    all_d = [x["factor_e3"] - x["factor_kn"] for x in rows]
    return {"n_candidates": len(cand), "n_changed": len(rows),
            "dfactor_e3_minus_kn": describe(all_d),
            "n_penalty_up": sum(1 for x in all_d if x > 0), "n_penalty_down": sum(1 for x in all_d if x < 0),
            "relative_penalty_change": describe([x["penalty_e3"] / x["penalty_kn"] - 1 for x in rows]),
            "reactions_differentiated_between_members": {
                arm: sum(1 for b, per in by_base.items() if len(per) >= 2 and
                         len({round(vec[r], 9) for r in per.values()}) > 1)
                for arm, vec in (("stock", stock), ("kn", kn), ("e3", e3))},
            "n_bases_in_multiple_member_compartments": multi,
            "n_cheapest_member_flips_kn_vs_e3": len(flips), "cheapest_member_flips": flips[:25],
            "per_member": pm_out, "largest_changes": rows[:15]}


def describe(xs):
    if not xs:
        return {"n": 0}
    s = sorted(xs)
    q = lambda f: s[min(len(s) - 1, int(round(f * (len(s) - 1))))]
    return {"n": len(s), "mean": round(statistics.mean(s), 6), "median": round(q(0.5), 6),
            "p05": round(q(0.05), 6), "p25": round(q(0.25), 6), "p75": round(q(0.75), 6),
            "p95": round(q(0.95), 6), "min": round(s[0], 6), "max": round(s[-1], 6)}


def p_lookup_from_scores(scores_by_arm):
    out = {}
    for arm, sc in scores_by_arm.items():
        out[arm] = {idx: {b: max(g.values()) if g else 0.0 for b, g in tab.items()} for idx, tab in sc.items()}
    return out


# ------------------------------------------------------------------ community sets (A, B)
def community(rids, ab):
    setup_env()
    from cobra.io import read_sbml_model
    from mscommunity import build_from_species_models
    import demo_compartment_scoring as D
    members = []
    for rid in rids:
        m = read_sbml_model(f"{MODELS}/kn/models_{TAG}/{rid}.xml")
        m.id = rid
        members.append(m)
    tot = sum(ab[r] for r in rids)
    norm = {r: ab[r] / tot for r in rids}
    comm = build_from_species_models(members, f"day{DAY}_ds", f"day{DAY}_ds", abundances=norm, printing=False)
    comm.id = f"day{DAY}_ds"
    comm = D.sanitize(comm)
    return comm, norm


def try_gapfill(comm, template_factory, media, scores, timeout_s):
    """One MSGapfill.run_gapfilling on the community: stock (scores None) or compartment-aware."""
    import io
    import logging
    import compartment_aware_gapfill as CAG
    from modelseedpy.core.msgapfill import MSGapfill
    buf = io.StringIO()
    h = logging.StreamHandler(buf)
    h.setLevel(logging.INFO)
    lg = logging.getLogger("modelseedpy")
    old = lg.level
    lg.setLevel(logging.INFO)
    lg.addHandler(h)
    rec = {}
    t0 = time.time()
    try:
        gf = MSGapfill(comm.copy(), default_gapfill_templates=[template_factory()], test_conditions=[])
        gf.gfmodel.solver.configuration.timeout = timeout_s
        if scores is not None:
            CAG.install()
            pkg = CAG.attach_scores(gf, scores)
            pkg.compute_gapfilling_penalties()
            pkg.build_gapfilling_objective_function()
            rec["scoring_stats"] = getattr(pkg, "compartment_scoring_stats", None)
        sol = gf.run_gapfilling(media=media, target="bio1", minimum_obj=0.01, binary_check=False,
                                prefilter=True)
        rec["solved"] = bool(sol)
        if sol:
            rec["n_new"] = len(sol.get("new", {}))
            rec["n_reversed"] = len(sol.get("reversed", {}))
        ts = getattr(gf.gfpkgmgr.getpkg("GapfillingPkg"), "test_solution", None)
        rec["test_solution_status"] = getattr(ts, "status", None)
        rec["test_solution_objective"] = getattr(ts, "objective_value", None)
    except Exception as exc:
        rec["error"] = f"{type(exc).__name__}: {str(exc)[:300]}"
    finally:
        CAG.uninstall()
        lg.removeHandler(h)
        lg.setLevel(old)
    rec["secs"] = round(time.time() - t0, 1)
    msgs = [l for l in buf.getvalue().splitlines() if l.strip()]
    rec["log_tail"] = msgs[-8:]
    return rec


def run_community_set(name, rids, ab, idx_assign, gapfill_timeout, attempt_gapfill=True):
    B = import_B(C.BUILDS["codiffusion"]["kn"])
    t0 = time.time()
    comm, norm = community(rids, ab)
    media = B._day_media()[B.MEDIA_ASSIGN.get(DAY, B.DEFAULT_MEDIA)[0]]
    tf = lambda: B._load_template("gram_neg")
    scores = {arm: scores_for(rids, arm) for arm in ("kn", "e3")}
    pen, stats, secs = penalty_vectors(comm, tf, media, scores)
    index_to_member = {str(i): r for i, r in enumerate(rids, 1)}
    res = {"members": rids, "abundance_normalized": {r: round(v, 6) for r, v in norm.items()},
           "media": media.id, "community_n_rxn": len(comm.reactions), "community_n_met": len(comm.metabolites),
           "n_member_scored_reactions": {arm: {index_to_member[k]: len(v) for k, v in sc.items()}
                                         for arm, sc in scores.items()},
           "scoring_stats": stats, "penalty_secs": secs,
           "kn_vs_e3": analyse(pen, scores, index_to_member, p_lookup_from_scores(scores)),
           "demo_comparison": None}
    # stock factor is 2 for every candidate (no reaction_scores): check it on the stock vector
    res["stock_equal_kn_for_unscored"] = sum(
        1 for r in pen["stock"] if r in pen["kn"] and abs(pen["stock"][r] - pen["kn"][r]) <= TOL)
    gf = {}
    for label, sc in ((("stock", None), ("kn", scores["kn"]), ("e3", scores["e3"])) if attempt_gapfill else ()):
        C.log(f"  set {name}: community gapfill attempt [{label}]")
        gf[label] = try_gapfill(comm, tf, media, sc, gapfill_timeout)
        C.log(f"    -> {gf[label].get('solved')} {gf[label].get('error', '')} ({gf[label]['secs']} s)")
    res["community_gapfill"] = gf
    res["secs_total"] = round(time.time() - t0, 1)
    # per-compartment penalties kept for the single-member equivalence check
    res["_pen"] = {arm: {r: v for r, v in vec.items() if IDX.search(r) and IDX.search(r).group(2) in index_to_member}
                   for arm, vec in pen.items()}
    return res


# ------------------------------------------------------------------ single-member (set C)
_S = {}


def _single_init():
    B = import_B(C.BUILDS["codiffusion"]["kn"])
    C.check_imports()
    import logging
    for name in ("modelseedpy", "cobra", "optlang"):
        logging.getLogger(name).setLevel(logging.ERROR)
    _S["B"] = B
    _S["media"] = B._day_media()[B.MEDIA_ASSIGN.get(DAY, B.DEFAULT_MEDIA)[0]]


def _single_one(rid):
    from cobra.io import read_sbml_model
    B = _S["B"]
    t0 = time.time()
    try:
        model = read_sbml_model(f"{MODELS}/kn/models_{TAG}/{rid}.xml")
        model.id = rid
        scores = {arm: scores_for([rid], arm, single_index=0) for arm in ("kn", "e3")}
        pen, stats, _ = penalty_vectors(model, lambda: B._load_template("gram_neg"), _S["media"], scores)
        res = analyse(pen, scores, {"0": rid}, p_lookup_from_scores(scores))
        res.pop("cheapest_member_flips", None)
        res["scoring_stats"] = stats
        res["n_member_scored_reactions"] = len(scores["kn"]["0"])
        out = {"rid": rid, "ok": True, "result": res, "secs": round(time.time() - t0, 1),
               "_pen": {arm: {r: v for r, v in vec.items() if IDX.search(r) and IDX.search(r).group(2) == "0"}
                        for arm, vec in pen.items()}}
    except Exception as exc:
        out = {"rid": rid, "ok": False, "error": f"{type(exc).__name__}: {str(exc)[:300]}",
               "secs": round(time.time() - t0, 1)}
    path = f"{WORK}/single/{rid}.json"
    C.dump(out, path)
    return {k: v for k, v in out.items() if k != "_pen"}


def run_singles(rids, workers):
    import multiprocessing as mp
    todo = [r for r in rids if not os.path.exists(f"{WORK}/single/{r}.json")]
    C.log(f"single-member penalties: {len(todo)} of {len(rids)} to run")
    if todo:
        with mp.get_context("spawn").Pool(workers, initializer=_single_init, maxtasksperchild=1) as pool:
            for r in pool.imap_unordered(_single_one, todo, chunksize=1):
                if not r["ok"]:
                    C.log(f"  single ERROR {r['rid']}: {r['error']}")
    return {r: json.load(open(f"{WORK}/single/{r}.json")) for r in rids}


def equivalence(comm_res, singles):
    """community compartment ci penalties == the member's single-model c0 penalties?"""
    out = {}
    for i, rid in enumerate(comm_res["members"], 1):
        s = singles.get(rid)
        if not s or not s.get("ok"):
            continue
        rec = {}
        for arm in ("stock", "kn", "e3"):
            cv = {base_of(r): v for r, v in comm_res["_pen"][arm].items() if IDX.search(r).group(2) == str(i)
                  and IDX.search(r).group(1) == "c"}
            sv = {base_of(r): v for r, v in s["_pen"][arm].items() if IDX.search(r).group(1) == "c"}
            common = set(cv) & set(sv)
            rec[arm] = {"n_comm": len(cv), "n_single": len(sv), "n_common": len(common),
                        "n_equal": sum(1 for b in common if abs(cv[b] - sv[b]) <= 1e-6)}
        out[rid] = rec
    return out


# ------------------------------------------------------------------ cheapest-member analysis
def flips_from_singles(members, singles):
    """Which member is cheapest for each reaction offered in >= 2 members' compartments.

    Uses the single-member penalty vectors (c0), which equal the members' community
    compartments exactly (equivalence block). Classes of a k/n -> E3 change of the argmin set:
      unique_changed  one cheapest member under each and they differ
      tie_broken      E3's cheapest set is a strict subset of k/n's (E3 breaks an exact tie)
      tie_created     k/n's cheapest set is a strict subset of E3's
      other           any other change of the set
    """
    pen = {}
    for rid in members:
        s = singles.get(rid)
        if s and s.get("ok"):
            pen[rid] = {arm: {base_of(r): v for r, v in s["_pen"][arm].items() if IDX.search(r).group(1) == "c"}
                        for arm in ("stock", "kn", "e3")}
    bases = defaultdict(list)
    for rid, p in pen.items():
        for b in p["kn"]:
            bases[b].append(rid)
    counts = defaultdict(int)
    examples = defaultdict(list)
    n_multi = 0
    for b, rids in bases.items():
        if len(rids) < 2:
            continue
        n_multi += 1
        am = {}
        for arm in ("stock", "kn", "e3"):
            lo = min(pen[r][arm][b] for r in rids)
            am[arm] = sorted(r for r in rids if abs(pen[r][arm][b] - lo) <= TOL)
        for ref in ("kn", "stock"):
            a, e = set(am[ref]), set(am["e3"])
            if a == e:
                continue
            if len(a) == 1 and len(e) == 1:
                cls = "unique_changed"
            elif e < a:
                cls = "tie_broken"
            elif a < e:
                cls = "tie_created"
            else:
                cls = "other"
            counts[f"{ref}_vs_e3:{cls}"] += 1
            if ref == "kn" and len(examples[cls]) < 8:
                shown = sorted(set(am["kn"]) ^ set(am["e3"]))[:6] + sorted(set(am["kn"]) & set(am["e3"]))[:2]
                examples[cls].append({"reaction": b, "n_cheapest_kn": len(am["kn"]), "n_cheapest_e3": len(am["e3"]),
                                      "n_cheapest_stock": len(am["stock"]),
                                      "cheapest_kn": am["kn"][:10], "cheapest_e3": am["e3"][:10],
                                      "penalty_kn": {r: round(pen[r]["kn"][b], 4) for r in shown},
                                      "penalty_e3": {r: round(pen[r]["e3"][b], 4) for r in shown}})
        if set(am["kn"]) != set(am["stock"]):
            counts["stock_vs_kn:changed"] += 1
    return {"n_members": len(pen), "n_reactions_offered_in_2plus_members": n_multi,
            "counts": dict(sorted(counts.items())), "examples_kn_vs_e3": dict(examples),
            "note": "an exact tie is broken arbitrarily by the MILP; unique_changed is the clean count"}


def post_flips(res, ab):
    singles = {}
    for f in os.listdir(f"{WORK}/single"):
        if f.endswith(".json"):
            singles[f[:-5]] = json.load(open(f"{WORK}/single/{f}"))
    out = {name: flips_from_singles(c["members"], singles) for name, c in res["communities"].items() if c}
    out["day174_all_members"] = flips_from_singles(list(ab), singles)
    return out


# ------------------------------------------------------------------ main
def mapping(sets, man_kn, man_e3, rids):
    rows = {}
    for rid in rids:
        s = sets[rid]
        a = s["rep_asv"]
        ek, ee = man_kn.get(a), man_e3.get(a)
        built = bool(ek and ek.get("status") == "built" and
                     os.path.exists(f"{C.BUILDS['codiffusion']['kn']}/genome_objects/{a}.json") and
                     os.path.exists(f"{C.BUILDS['codiffusion']['e3']}/genome_objects/{a}.json"))
        same_sel = [m for m in s["member_asvs"]
                    if man_kn.get(m, {}).get("status") == "built" and
                    sorted(man_kn[m]["genomes"]) == sorted(ek["genomes"])] if built else []
        rows[rid] = {"rep_asv": a, "mapped": built,
                     "n_genomes_20hit": len(s["genomes"]),
                     "n_genomes_d100": ek.get("n_genomes") if ek else None,
                     "d100_selection_equals_20hit": bool(ek) and sorted(ek.get("genomes", [])) == sorted(s["genomes"]),
                     "e3_applied": bool(ee and (ee.get("pfunc") or {}).get("applied")),
                     "e3_clade_rank": ((ee or {}).get("pfunc") or {}).get("clade_rank"),
                     "n_member_asvs": len(s["member_asvs"]),
                     "n_member_asvs_with_rep_d100_selection": len(same_sel)}
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--gapfill-timeout", type=int, default=1800)
    ap.add_argument("--skip-community", action="store_true", help="smoke: no sets A/B")
    ap.add_argument("--no-gapfill-attempt", action="store_true")
    ap.add_argument("--post-only", action="store_true",
                    help="recompute the cheapest-member block of results/penalties.json from the saved "
                         "single-member penalty vectors")
    ap.add_argument("--limit", type=int, default=None,
                    help="smoke: drafts/single-member runs for sets A+B plus the first N other members")
    a = ap.parse_args()
    if a.post_only:
        res = json.load(open(f"{C.OUT}/results/penalties.json"))
        res["cheapest_member"] = post_flips(res, json.load(open(AB_PATH))[DAY]["abundances"])
        res["cheapest_member"]["_written_by"] = {"mode": "penalties.py --post-only (saved single-member vectors)",
                                                 "time": time.strftime("%Y-%m-%dT%H:%M:%S"),
                                                 "penalties_py_sha256": C.sha256_file(os.path.abspath(__file__), 16)}
        C.dump(res, f"{C.OUT}/results/penalties.json")
        C.log("cheapest-member block rewritten")
        return
    setup_env()
    C.ensure_pins()
    os.makedirs(WORK, exist_ok=True)
    t_all = time.time()
    sets = json.load(open(SETS_PATH)); sets.pop("_meta", None)
    ab = json.load(open(AB_PATH))[DAY]["abundances"]
    man_kn = json.load(open(f"{C.BUILDS['codiffusion']['kn']}/synthetic_genome_manifest.json"))
    man_e3 = json.load(open(f"{C.BUILDS['codiffusion']['e3']}/synthetic_genome_manifest.json"))
    day_members = list(ab)                                   # abundance order (as the demo)
    mp_all = mapping(sets, man_kn, man_e3, list(sets))
    mp_day = {r: mp_all[r] for r in day_members}
    unmapped = [r for r in day_members if not mp_day[r]["mapped"]]
    demo = json.load(open(DEMO_OUT))
    demo_prov = list(json.load(open(DEMO_PROV))["abundances"])
    set_a = demo["members"]
    assert set_a == demo_prov[:4] == day_members[:4], (set_a, demo_prov[:4], day_members[:4])
    set_b = [r for r in day_members if mp_day[r]["e3_applied"] and mp_day[r]["mapped"]][:4]
    C.log(f"day {DAY}: {len(day_members)} members, {len(unmapped)} unmapped, "
          f"{sum(v['e3_applied'] for v in mp_day.values())} reweighted | A={set_a} | B={set_b}")

    # template classes on the d100 genome objects (production rule, our own output file)
    C.dump({"_meta": {"source": "subset of " + SETS_PATH}, **{r: sets[r] for r in day_members}}, MY_SETS)
    B = import_B(C.BUILDS["codiffusion"]["kn"])
    imports = C.check_imports()
    if not os.path.exists(MY_ASSIGN):
        B.stage_templates()
    assign = json.load(open(MY_ASSIGN)); assign_meta = assign.pop("_meta")
    prod = json.load(open(PROD_ASSIGN)); prod.pop("_meta", None)
    tpl_diff = {r: [prod.get(r, {}).get("template"), assign[r]["template"]] for r in day_members
                if prod.get(r, {}).get("template") != assign[r]["template"]}

    # drafts
    tasks = [(r, sets[r]["rep_asv"]) for r in day_members if mp_day[r]["mapped"]]
    if a.limit is not None:
        ab_sets = set(set_a) | set(set_b)
        keep = ab_sets | {r for r in [r for r, _ in tasks if r not in ab_sets][:a.limit]}
        tasks = [t for t in tasks if t[0] in keep]
    t0 = time.time()
    draft_log = build_drafts(tasks, a.workers)
    draft_secs = round(time.time() - t0, 1)
    C.dump(draft_log, f"{WORK}/draft_log.json")
    bad = [r for r in draft_log if r["result"].startswith("ERROR")]
    C.log(f"drafts done in {draft_secs} s, {len(bad)} errors")
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    both = [rid for rid, _ in tasks
            if all(os.path.exists(f"{MODELS}/{arm}/models_{TAG}/{rid}.xml") for arm in ("kn", "e3"))]
    with ProcessPoolExecutor(max_workers=a.workers, mp_context=mp.get_context("spawn")) as ex:
        cmp = dict(zip(both, ex.map(compare_drafts, both, chunksize=2)))
    C.dump(cmp, f"{WORK}/draft_comparison.json")
    C.log("draft comparison: " + json.dumps({
        "n": len(cmp), "same_ids": sum(v["same_reaction_ids"] for v in cmp.values()),
        "struct_diff": sum(v["n_struct_diff"] for v in cmp.values()),
        "rxn_p_changed": sum(v["n_rxn_p_changed"] for v in cmp.values())}))

    # set C first (parallel), then the two communities (serial, larger)
    singles = run_singles([r for r, _ in tasks], a.workers)
    comm = {}
    for name, rids in (("A_demo_top4", set_a), ("B_top4_reweighted", set_b)):
        C.log(f"set {name}: {rids}")
        r = (run_community_set(name, rids, ab, None, a.gapfill_timeout, not a.no_gapfill_attempt)
             if not a.skip_community else None)
        comm[name] = r
    if "A_demo_top4" in comm and comm["A_demo_top4"]:
        A = comm["A_demo_top4"]
        A["demo_comparison"] = {
            "demo_n_candidates_20hit_v7_drafts": demo["n_candidates"],
            "demo_scoring_stats": demo["scoring_stats"],
            "demo_reactions_differentiated": [demo["reactions_differentiated_stock"],
                                              demo["reactions_differentiated_patched"]],
            "note": "the demo priced 54d1150 (20-hit) v7 drafts on the ModelSEEDTemplates tree of "
                    "2026-09; this run prices d100 drafts on codif-curated b6fe86b, so the "
                    "candidate counts differ by construction"}
    equiv = {name: equivalence(r, singles) for name, r in comm.items() if r}
    for r in comm.values():
        if r:
            r.pop("_pen", None)

    # set C summary over the whole community
    per = {}
    pooled_d, pooled_rel = [], []
    for rid, _ in tasks:
        s = singles.get(rid)
        if not s or not s.get("ok"):
            per[rid] = {"error": (s or {}).get("error", "missing")}
            continue
        res = s["result"]
        pm = res["per_member"].get(rid, {})
        per[rid] = {"n_genomes_d100": mp_day[rid]["n_genomes_d100"], "e3_applied": mp_day[rid]["e3_applied"],
                    "abundance": round(ab[rid], 6), "n_candidates": res["n_candidates"],
                    "n_member_scored_candidates": pm.get("n_member_scored"),
                    "n_member_scored_reactions": res["n_member_scored_reactions"],
                    "n_changed": res["n_changed"], "n_penalty_up": res["n_penalty_up"],
                    "n_penalty_down": res["n_penalty_down"],
                    "dfactor_e3_minus_kn": res["dfactor_e3_minus_kn"],
                    "relative_penalty_change": res["relative_penalty_change"],
                    "largest_changes": res["largest_changes"][:3],
                    "draft_rxn_p_changed": cmp.get(rid, {}).get("n_rxn_p_changed")}
        sx = json.load(open(f"{WORK}/single/{rid}.json"))
        pen = sx["_pen"]
        for r in pen["kn"]:
            if r in pen["e3"] and abs(pen["e3"][r] - pen["kn"][r]) > TOL:
                base = pen["stock"][r] / 2
                pooled_d.append(pen["e3"][r] / base - pen["kn"][r] / base)
                pooled_rel.append(pen["e3"][r] / pen["kn"][r] - 1)
    n1 = [r for r in per if not mp_day[r]["e3_applied"]]
    summary_c = {
        "n_members": len(per), "n_reweighted": sum(1 for r in per if mp_day[r]["e3_applied"]),
        "n_members_with_changed_candidates": sum(1 for r in per if per[r].get("n_changed")),
        "n1_members_with_any_change": sum(1 for r in n1 if per[r].get("n_changed")),
        "n_changed_candidates_total": sum(per[r].get("n_changed") or 0 for r in per),
        "n_member_scored_candidates_total": sum(per[r].get("n_member_scored_candidates") or 0 for r in per),
        "n_candidates_total": sum(per[r].get("n_candidates") or 0 for r in per),
        "dfactor_e3_minus_kn_pooled": describe(pooled_d),
        "relative_penalty_change_pooled": describe(pooled_rel),
        "n_penalty_up": sum(1 for x in pooled_d if x > 0), "n_penalty_down": sum(1 for x in pooled_d if x < 0),
        "abundance_share_of_members_with_changes": round(sum(ab[r] for r in per if per[r].get("n_changed")), 6)}

    res = {"_meta": {"script": "penalties.py", "day": DAY, "tag": TAG,
                     "code": C.code_record([os.path.abspath(__file__)]), "pins": C.pin_record(),
                     "templates_root": C.TPL_ROOT["codif_curated"], "template_assignment": assign_meta,
                     "imports": imports,
                     "secs_total": round(time.time() - t_all, 1), "draft_secs": draft_secs,
                     "gapfill_timeout_s": a.gapfill_timeout,
                     "definitions": {
                         "penalty": "GapfillingPkg base penalty x factor; stock factor = 2 for every candidate "
                                    "(no reaction_scores); compartment-aware factor = 2 - max p of member i's "
                                    "genes for the candidate's base reaction (unscored: 2)",
                         "dfactor_e3_minus_kn": "factor(E3) - factor(k/n) = p_kn - p_e3, over candidates whose "
                                                "compartment-aware penalty differs between the arms",
                         "relative_penalty_change": "penalty(E3) / penalty(k/n) - 1"}},
           "mapping": {"rule": "member (20-hit ASVset) -> d100 genome object of the same representative ASV md5",
                       "n_asvsets_all_days": len(mp_all),
                       "n_unmapped_all_days": sum(1 for v in mp_all.values() if not v["mapped"]),
                       "n_reweighted_all_days": sum(1 for v in mp_all.values() if v["e3_applied"]),
                       "n_rep_selection_changed_20hit_to_d100_all_days":
                           sum(1 for v in mp_all.values() if not v["d100_selection_equals_20hit"]),
                       "day174": mp_day, "day174_unmapped": unmapped,
                       "template_class_differs_from_production_assignment": tpl_diff},
           "drafts": {"n_built": len(cmp), "errors": bad,
                      "all_same_reaction_ids": all(v["same_reaction_ids"] for v in cmp.values()),
                      "n_struct_diff_total": sum(v["n_struct_diff"] for v in cmp.values()),
                      "n_gpr_string_order_only_total": sum(v["n_gpr_string_order_only"] for v in cmp.values()),
                      "n_other_notes_diff_total": sum(v["n_other_notes_diff"] for v in cmp.values()),
                      "per_member": {r: {k: v for k, v in x.items() if k != "rxn_p_changed"}
                                     for r, x in cmp.items()}},
           "communities": comm, "equivalence_single_vs_community": equiv,
           "single_member_day174": {"summary": summary_c, "per_member": per}}
    res["cheapest_member"] = post_flips(res, ab)
    path = f"{C.OUT}/results/penalties.json" if a.limit is None else f"{WORK}/penalties_smoke.json"
    C.dump(res, path)
    C.log(f"done in {res['_meta']['secs_total']} s -> {path}")


if __name__ == "__main__":
    main()
