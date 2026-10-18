#!/usr/bin/env python
"""Phenotype arm of the prFBA / PICRUSt2 comparison: EC set -> ModelSEED model -> phenotype.

Three models per held-out organism, built by an IDENTICAL procedure that differs only
in the input EC set:

  truth      the held-out genome's own BV-BRC EC set (genome_function_sets.json["ec"])
  prfba      prFBA's predicted EC set for the same amplicon
  picrust2   PICRUSt2's predicted EC set for the same amplicon

EC -> model input (the one mapping, applied identically to all three):
  a ModelSEED template role is granted iff EVERY complete EC number appearing in the
  role's name is present in the input EC set (65 / 7463 EC-bearing roles name more than
  one EC; requiring all of them never credits an activity the input lacks).  Each granted
  role becomes one synthetic MSFeature whose id is the template role id and whose
  ontology_terms["RAST"] is the role name verbatim -- MSBuilder's normalize_role() keeps
  EC digits, so the match into the template is exact.  MSBuilder then applies the normal
  ModelSEED complex/GPR rules.  Measured coverage: 96.4% (bacteria) / 94.5% (archaea) of a
  genome's ECs name at least one template role.

Template choice is a pure function of the held-out genome's TRUE lineage, so it is
identical for all three models and is never under test.  auto_select_template() is
deliberately NOT used (it downloads a classifier and would make the template depend on
the input gene set).

Un-gapfilled drafts are NOT usable for flux phenotypes: bio1 is 0 even with every
exchange open, and naive ATP-yield tests are corrupted by template artifact loops (see
module docstring of run_smoke()).  Phenotypes are therefore read from the standard
ModelSEED pipeline (MSATPCorrection + MSGapfill, i.e. MSBuilder.build_and_gapfill),
run with identical settings on all three inputs, with the pre-gapfill draft reaction
sets and the gapfill burden reported alongside.
"""
import json
import os
import re
import time

import cobra
import pandas as pd
from cobra.core import Metabolite, Reaction

from modelseedpy.core.msbuilder import MSBuilder
from modelseedpy.core.msgenome import MSFeature, MSGenome
from modelseedpy.core.msmedia import MSMedia
from modelseedpy.core.mstemplate import MSTemplateBuilder

MSD = "/home/freiburger/Documents/ModelSEEDDatabase"
MST = "/home/freiburger/Documents/ModelSEEDTemplates/templates/v7.0/"
TEMPLATE_FILES = {"gramneg": "GramNegModelTemplateV7.json",
                  "grampos": "GramPosModelTemplateV7.json",
                  "archaea": "ArchaeaTemplateV6.json",
                  "core": "Core-V6.json"}
ATP_MEDIAS = "/home/freiburger/Documents/ModelSEEDpy/modelseedpy/data/atp_medias.tsv"
KB_MEDIA = f"{MSD}/Media/media_list_with_meta.txt"
KB_COMPLETE_CPD = f"{MSD}/Media/KBaseMedia.cpd"

EC_RE = re.compile(r"\(EC:?\s*([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)\)", re.I)
# GTDB/BV-BRC phyla reconstructed with the Gram-positive template.
GRAMPOS_PHYLA = {"Bacillota", "Bacillota_A", "Bacillota_B", "Bacillota_C", "Bacillota_G",
                 "Firmicutes", "Actinomycetota", "Actinobacteria", "Chloroflexota",
                 "Deinococcota", "Thermotogota", "Mycoplasmatota", "Tenericutes"}

_tcache = {}


def preload_biochemistry(path=MSD):
    """modelseedpy's config.cfg points biochem at the RELATIVE path ../ModelSEEDDatabase,
    so ModelSEEDBiochem.get() raises FileNotFoundError unless cwd happens to be
    ~/Documents/ModelSEEDpy.  MSATPCorrection.expand_model_to_genome_scale() needs it.
    Preload the singleton instead of editing the user's installed config (~13 s)."""
    from modelseedpy.biochem import modelseed_biochem as MB
    if MB.ModelSEEDBiochem.default_biochemistry is None:
        MB.ModelSEEDBiochem.default_biochemistry = MB.from_local(path)
    return MB.ModelSEEDBiochem.default_biochemistry


def load_template(key):
    """-> (MSTemplate, {role_id: (role_name, frozenset(ECs in name))})"""
    if key not in _tcache:
        with open(MST + TEMPLATE_FILES[key]) as fh:
            d = json.load(fh)
        t = MSTemplateBuilder.from_dict(d).build()
        role_ecs = {}
        for r in d["roles"]:
            ecs = frozenset(EC_RE.findall(r["name"]))
            if ecs:
                role_ecs[r["id"]] = (r["name"], ecs)
        _tcache[key] = (t, role_ecs)
    return _tcache[key]


def template_for(rec):
    """Template key from the held-out genome's true lineage (selection_10k.json record)."""
    if (rec.get("domain") or "").lower().startswith("arch"):
        return "archaea"
    return "grampos" if rec.get("phylum") in GRAMPOS_PHYLA else "gramneg"


def genome_from_ec(ec_set, role_ecs, genome_id="g"):
    ec_set = set(ec_set)
    g = MSGenome()
    g.id = genome_id
    feats = []
    for rid, (name, ecs) in sorted(role_ecs.items()):
        if ecs <= ec_set:
            f = MSFeature(rid, "")
            f.add_ontology_term("RAST", name)
            feats.append(f)
    if feats:
        g.add_features(feats)
    return g


def build_draft(ec_set, tmpl_key, model_id):
    tmpl, role_ecs = load_template(tmpl_key)
    g = genome_from_ec(ec_set, role_ecs, model_id)
    b = MSBuilder(g, tmpl, ontology_term="RAST",
                  template_core=load_template("core")[0])
    m = b.build(model_id, "0", allow_all_non_grp_reactions=False,
                annotate_with_rast=False, biomass_classic=False, biomass_gc=0.5,
                add_reaction_from_rast_annotation=True)
    return m, b, len(g.features)


def complete_media():
    cpds = [l.strip() for l in open(KB_COMPLETE_CPD) if l.strip()]
    m = MSMedia.from_dict({c: 100.0 for c in cpds})
    m.id = "complete"
    return m, {c: (-100.0, 1000.0) for c in cpds}


def build_gapfilled(ec_set, tmpl_key, model_id, gf_media):
    """Standard ModelSEED reconstruction: draft -> ATP correction -> gapfill on gf_media."""
    tmpl, role_ecs = load_template(tmpl_key)
    g = genome_from_ec(ec_set, role_ecs, model_id)
    b = MSBuilder(g, tmpl, ontology_term="RAST",
                  template_core=load_template("core")[0])
    mdl, atpc, tests, gapsol = b.build_and_gapfill(model_id, gf_media,
                                                   annotate_with_rast=False)
    return mdl, gapsol, atpc, tests


# ---------------------------------------------------------------- media / phenotypes
def load_atp_media_panel(path=ATP_MEDIAS):
    df = pd.read_csv(path, sep="\t", index_col=0)
    panel = {}
    for mid in df.columns:
        cpds = {ex.split("_")[1]: (-1000.0, 1000.0) for ex, v in df[mid].items() if v > 0}
        cpds["cpd00001"] = (-1000.0, 1000.0)
        cpds["cpd00067"] = (-1000.0, 1000.0)
        panel[mid] = cpds
    return panel


def load_kbase_media(path=KB_MEDIA, prefix=None, minimal_only=True):
    """KBase/ModelSEED curated media: 329 defined minimal Biolog media
    (168 Carbon-*, plus Nitrogen-*/Phosphate-*/Sulfate-*)."""
    out = {}
    with open(path) as fh:
        fh.readline()
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) < 6:
                continue
            name, minimal = f[0], f[3]
            if minimal_only and minimal != "1":
                continue
            if prefix and not name.startswith(prefix):
                continue
            cpds = {}
            for tok in f[5:]:
                p = tok.split(";")
                if len(p) >= 3 and p[0]:
                    # KBase media columns are minFlux;maxFlux in KBase sign convention
                    # (positive = uptake). cobra exchange bounds are the negation:
                    # lower = -maxFlux (uptake cap), upper = -minFlux (secretion cap).
                    # Cf. MSMedia.from_kbase_object and KBaseMediaPkg.build_package.
                    cpds[p[0]] = (-float(p[2]), -float(p[1]))
            if cpds:
                out[name] = cpds
    return out


def ensure_exchanges(model, cpd_ids):
    """Add the e0 metabolite + EX_ reaction for media compounds the model lacks.
    Annotation-independent (same media list for all three models); uptake into c0 still
    needs a gene-associated transport reaction, so the phenotype stays gene-gated."""
    added = []
    for c in sorted(cpd_ids):
        mid = f"{c}_e0"
        if mid not in model.metabolites:
            model.add_metabolites([Metabolite(mid, None, mid, 0, "e0")])
        ex = f"EX_{mid}"
        if ex not in model.reactions:
            r = Reaction(ex, "Exchange for " + mid, "exchanges", -1000, 1000)
            r.add_metabolites({model.metabolites.get_by_id(mid): -1})
            model.add_reactions([r])
            added.append(ex)
    return added


def set_media_bounds(model, cpds, default_excretion=100.0):
    """KBaseMediaPkg semantics: every exchange defaults to no uptake / free excretion,
    then the media compounds get their own bounds."""
    for r in model.reactions:
        if r.id.startswith("EX_"):
            r.bounds = (0.0, default_excretion)
    for c, b in cpds.items():
        ex = f"EX_{c}_e0"
        if ex in model.reactions:
            lb, ub = b if isinstance(b, tuple) else (-float(b), 1000.0)
            model.reactions.get_by_id(ex).bounds = (lb, max(ub, 0.0))


def _val(s, thresh):
    return 0.0 if (s is None or s != s or s < thresh) else float(s)


def growth_vector(model, media, target="bio1", thresh=1e-5):
    ensure_exchanges(model, {c for m in media.values() for c in m})
    out = {}
    with model:
        model.objective = target
        for mid, cpds in media.items():
            set_media_bounds(model, cpds)
            out[mid] = _val(model.slim_optimize(), thresh)
    return out


def atp_yield_vector(model, panel=None, thresh=1e-6):
    """ATP produced per medium (maximize ATPM_c0).  WARNING: on un-ATP-corrected models
    this is corrupted by template artifact loops -- measured 1000 ATP on the 'empty'
    medium via rxn15249/rxn20583 (duplicate alpha-glucose stereochemistry).  Only
    meaningful after MSATPCorrection."""
    panel = panel or load_atp_media_panel()
    MSBuilder.add_atpm(model)
    ensure_exchanges(model, {c for m in panel.values() for c in m})
    out = {}
    with model:
        model.objective = "ATPM_c0"
        for mid, cpds in panel.items():
            set_media_bounds(model, cpds)
            out[mid] = _val(model.slim_optimize(), thresh)
    return out


def precursor_producibility(model, cpds, bio="bio1", thresh=1e-6):
    """Can each biomass precursor be made from this medium? (demand-reaction probe)"""
    precursors = sorted(m.id for m, v in model.reactions.get_by_id(bio).metabolites.items()
                        if v < 0)
    ensure_exchanges(model, set(cpds))
    out = {}
    with model:
        set_media_bounds(model, cpds)
        for mid in precursors:
            if mid not in model.metabolites:
                out[mid] = 0.0
                continue
            with model:
                dm = Reaction("DM_probe", "demand", "demand", 0, 1000)
                dm.add_metabolites({model.metabolites.get_by_id(mid): -1})
                model.add_reactions([dm])
                model.objective = "DM_probe"
                out[mid] = _val(model.slim_optimize(), thresh)
    return out


def essentiality(model, cpds, target="bio1", thresh=1e-5, frac=0.01):
    """Single-reaction deletion on one medium -> {rxn_id: essential?}, plus wild-type flux.
    Returns (None, 0.0) if the model does not grow on the medium."""
    ensure_exchanges(model, set(cpds))
    ess = {}
    with model:
        model.objective = target
        set_media_bounds(model, cpds)
        wt = model.slim_optimize()
        if wt is None or wt != wt or wt < thresh:
            return None, 0.0
        for r in [x for x in model.reactions
                  if not x.id.startswith(("EX_", "SK_", "DM_", "bio"))]:
            b = r.bounds
            r.bounds = (0, 0)
            s = model.slim_optimize()
            r.bounds = b
            ess[r.id] = bool(s is None or s != s or s < frac * wt)
    return ess, float(wt)


def core_atp_vector(ec_set, tmpl_key, model_id, no_gapfilling=True):
    """The cheap, loop-controlled, gapfill-free phenotype: ATP yield on each of the 54
    ModelSEEDpy ATP media, computed by MSATPCorrection on the model restricted to the
    core template (252 reactions).  disable_noncore_reactions() removes the artifact
    loops that make a raw draft produce 1000 ATP from nothing -- measured ATP on the
    'empty' medium is 0.0 here.  2-3 s per model after templates are cached.

    Measured on 2866714.3: 11/54 media positive for the truth EC set, 1/54 after
    dropping 20% of the ECs (all pyruvate and acetate respiration lost), 11/54 with 30%
    false-positive ECs added (binary Hamming 0, mean |dyield| 0.02).
    """
    from modelseedpy.core.msatpcorrection import MSATPCorrection, load_default_medias
    tmpl, role_ecs = load_template(tmpl_key)
    core = load_template("core")[0]
    g = genome_from_ec(ec_set, role_ecs, model_id)
    b = MSBuilder(g, tmpl, ontology_term="RAST", template_core=core)
    mb = b.build_base_model(model_id, "0", False, False)
    MSBuilder.add_atpm(mb)
    atpc = MSATPCorrection(mb, core, load_default_medias(), compartment="c0",
                           atp_hydrolysis_id="ATPM_c0", load_default_medias=False)
    v = atpc.evaluate_growth_media(no_gapfilling=no_gapfilling)
    return {k: (0.0 if (x is None or x != x) else float(x)) for k, x in v.items()}


def sanity_gate(model, bio="bio1"):
    """Necessary but NOT sufficient: no biomass and no ATP from water and protons alone.
    Measured: the complete-medium-gapfilled truth model passes (0.0 / 0.0) and grows
    (4.779).  The same EC set gapfilled on glucose minimal FAILS (ATPM = 100 from water
    and protons) and does not grow on the very medium it was gapfilled for (biomass 0.0)
    -- gapfill on complete, not on a minimal medium.  The test is still only necessary,
    not sufficient: it cannot see a mass-generating loop that needs a substrate to prime
    it, so pair it with yield_check()."""
    ensure_exchanges(model, {"cpd00001", "cpd00067"})
    empty = {"cpd00001": (-100.0, 100.0), "cpd00067": (-100.0, 100.0)}
    with model:
        model.objective = bio
        set_media_bounds(model, empty)
        b = model.slim_optimize()
    with model:
        MSBuilder.add_atpm(model)
        model.objective = "ATPM_c0"
        set_media_bounds(model, empty)
        a = model.slim_optimize()
    b = 0.0 if (b is None or b != b) else float(b)
    a = 0.0 if (a is None or a != a) else float(a)
    return {"empty_biomass": b, "empty_atp": a, "pass": b < 1e-6 and a < 1e-6}


def yield_check(model, cpds, bio="bio1", max_ratio=0.5):
    """Second half of the validity gate: biomass flux must be small relative to the
    CARBON it is made from, to catch substrate-primed mass-generating loops that
    sanity_gate() cannot see.  Total uptake is reported too but must not be used: media
    ions are allowed in at 100 each and would mask the loop.
    Measured on the complete-medium-gapfilled truth model: biomass 4.779,
    total uptake 345.2, ratio_total 0.0138 (passes).  max_ratio still wants calibration
    on a few dozen models before it is used to discard anything."""
    ensure_exchanges(model, set(cpds))
    with model:
        model.objective = bio
        set_media_bounds(model, cpds)
        sol = model.optimize()
        if sol.status != "optimal":
            return {"biomass": 0.0, "total_uptake": 0.0, "carbon_uptake": 0.0,
                    "ratio_total": 0.0, "ratio_carbon": 0.0, "pass": True}
        tot, carb = 0.0, 0.0
        for rid, f in sol.fluxes.items():
            if not rid.startswith("EX_") or f >= 0:
                continue
            tot -= f
            mets = list(model.reactions.get_by_id(rid).metabolites)
            formula = mets[0].formula if mets else None
            n_c = 0
            if formula:
                m = re.search(r"C(\d*)(?![a-z])", formula)
                if m:
                    n_c = int(m.group(1) or 1)
            carb += -f * n_c
        return {"biomass": float(sol.objective_value), "total_uptake": tot,
                "carbon_uptake": carb,
                "ratio_total": (sol.objective_value / tot) if tot > 1e-9 else float("inf"),
                "ratio_carbon": (sol.objective_value / carb) if carb > 1e-9 else float("inf"),
                "pass": bool(carb > 1e-9 and sol.objective_value / carb < max_ratio)}


def rxn_sets(model):
    """(metabolic reactions, gene-associated reactions): the structural readout."""
    core = {r.id for r in model.reactions
            if not r.id.startswith(("EX_", "SK_", "DM_", "bio")) and r.id != "ATPM_c0"}
    gpr = {r.id for r in model.reactions if r.gene_reaction_rule}
    return core, gpr


def jaccard(a, b):
    a, b = set(a), set(b)
    return (len(a & b) / len(a | b)) if (a | b) else 1.0


def hamming(a, b, thresh=0.0):
    keys = set(a) | set(b)
    return sum(1 for k in keys if (a.get(k, 0) > thresh) != (b.get(k, 0) > thresh)), len(keys)


# ------------------------------------------------------------------ smoke test / CLI
def run_smoke(gid, rec, ec_sets, gf_media=None, growth_media=None, out=None):
    """Build one model per entry of ec_sets ({tag: [EC,...]}) for genome gid and compare
    every model against ec_sets['truth'].

    Measured on 2866714.3 (Oscillochloris sp. ZM17-4, Chloroflexota -> grampos,
    587 ECs), gapfill media = complete (347 compounds), 64-core box, Gurobi 13.0.1 via optlang:
      draft build (no gapfill)          1.5 s   980 metabolic rxns, 866 genes
      build_and_gapfill                 42-70 s 78-100 gapfilled rxns
      draft growth on complete media    bio1 = 0.0  (every exchange open)
      draft ATP yield on 'empty' media  1000 (artifact loop rxn15249/rxn20583)
      gapfilled growth on complete      4.78 (truth) vs 3.57 (20% ECs dropped)
      essentiality (complete)           266 vs 232 essential; 74/908 shared disagree
      ATP-correction media selected     11 (truth) vs 1 (20% ECs dropped)
      growth on 167 Carbon-* minimal    0 for every model, gapfilled on complete OR on
                                        glucose minimal -- the Biolog panel is dead here
    """
    preload_biochemistry()
    tk = template_for(rec)
    if gf_media is None:
        gf_media, gf_cpds = complete_media()
    else:
        gf_cpds = {c.id: (c.lower_bound, c.upper_bound) for c in gf_media.mediacompounds}
    res = {}
    for tag, ecs in ec_sets.items():
        t0 = time.time()
        draft, _, n_feat = build_draft(ecs, tk, f"{gid}_{tag}_draft")
        d_core, _ = rxn_sets(draft)
        mdl, gapsol, atpc, tests = build_gapfilled(ecs, tk, f"{gid}_{tag}", gf_media)
        g_core, _ = rxn_sets(mdl)
        ess, wt = essentiality(mdl, gf_cpds)
        res[tag] = {"n_ec": len(ecs), "n_feat": n_feat, "template": tk,
                    "draft_core": sorted(d_core), "gf_core": sorted(g_core),
                    "n_gapfilled": len(gapsol), "gapfilled": sorted(gapsol),
                    "growth_gf_media": wt, "ess": ess or {},
                    "atp_media": [m.id for m in atpc.selected_media],
                    "growth": growth_vector(mdl, growth_media) if growth_media else None,
                    "secs": round(time.time() - t0, 1)}
    T = res["truth"]
    for tag, R in res.items():
        if tag == "truth":
            continue
        both = set(T["ess"]) & set(R["ess"])
        R["vs_truth"] = {
            "draft_jaccard": round(jaccard(T["draft_core"], R["draft_core"]), 4),
            "gf_jaccard": round(jaccard(T["gf_core"], R["gf_core"]), 4),
            "growth_ratio": (R["growth_gf_media"] / T["growth_gf_media"]
                             if T["growth_gf_media"] else None),
            "ess_disagree": sum(1 for k in both if T["ess"][k] != R["ess"][k]),
            "ess_shared": len(both),
            "ess_truth_only_rxn": len(set(T["ess"]) - set(R["ess"])),
            "atp_media_jaccard": round(jaccard(T["atp_media"], R["atp_media"]), 4)}
        if growth_media:
            d, n = hamming(T["growth"], R["growth"])
            R["vs_truth"]["growth_disagree"] = [d, n]
    if out:
        with open(out, "w") as fh:
            json.dump(res, fh, indent=1)
    return res


if __name__ == "__main__":
    import argparse
    import random
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome", required=True, help="genome_id in genome_function_sets.json")
    ap.add_argument("--function-sets", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "data",
        "genome_function_sets.json"))
    ap.add_argument("--panel", default=os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "..", "data",
        "selection_10k.json"))
    ap.add_argument("--drop-frac", type=float, default=0.2,
                    help="synthetic 'prediction': drop this fraction of the truth ECs")
    ap.add_argument("--growth-media-prefix", default=None,
                    help="e.g. Carbon- to add a KBase minimal-media growth panel")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    import orjson
    with open(args.function_sets, "rb") as fh:
        gfs = orjson.loads(fh.read())
    with open(args.panel) as fh:
        rec = {r["genome_id"]: r for r in json.load(fh)["genomes"]}[args.genome]
    ec = sorted(gfs[args.genome]["ec"])
    random.seed(7)
    mut = sorted(random.sample(ec, int((1 - args.drop_frac) * len(ec))))
    gm = load_kbase_media(prefix=args.growth_media_prefix) if args.growth_media_prefix else None
    res = run_smoke(args.genome, rec, {"truth": ec, f"drop{int(args.drop_frac*100)}": mut},
                    growth_media=gm, out=args.out)
    for tag, r in res.items():
        print(tag, {k: v for k, v in r.items()
                    if k not in ("draft_core", "gf_core", "ess", "growth", "gapfilled")})
