#!/usr/bin/env python
"""Stage E: function-probability estimators for prFBA against the shipped k/n.

DESIGN (binding)
  ../PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md, section "Stage E" (committed b46598b before
  any Stage E number was computed). Declared POST HOC, after the calibration study's
  exploratory top-k analysis had been seen. Every deviation and every choice the design
  leaves open is listed in DEVIATIONS below and in the README this script writes.

POPULATION
  The controlled genus-blocked shared-reference benchmark exactly as calibrate_pfunc
  scored it: cells.csv.gz arm prfba_exact_tie, space ec, outcome scored, density 1.0
  (57,021 cells). Truth = the cell's source genome's EC set on the 2,235-EC universe
  (_cache_panel_ec.json through norm_ec, calibrate_pfunc.load_context verbatim). Reps =
  predictions.json pred.exact_tie[md5].reps of the cell's fold run. Holdout groups =
  cells.csv.gz group_id (== fold_plan genome_group, asserted).

RECORDS
  Each fold's own hits (data/prfba/f{fold}_d100/hits/asv_top20_alignment_hits.json),
  searched against that fold's reference only. The records hold 200 hits in align-score
  rank order (asserted); only the FIRST 100 are used (production records hold 100).
  select_references stage 0 (coverage >= 0.90, |SW - edlib| <= 0.02, identity >= 0.865),
  then ranked by identity, then align_score (stable), then deduplicated by genome (first
  occurrence) or additionally by BV-BRC species (select_references' species key on the
  hit lineage: _norm(Species) -> 'taxon:<taxon_id>' -> genome_id; the highest-identity
  genome of each species is kept). "Top K" = the first K of those survivors (K = 100 =
  every survivor of the first 100 hits). The 200-hit, genome-dedup list (first 100
  survivors, calibrate_pfunc.topk_job verbatim) is kept ONLY for the reproduction check.

ESTIMATORS (x_jf = 1 if neighbour j carries EC f; EC sets = the fold's
_cache_ec_subset.json, asserted equal to the panel truth matrix row by row)
  E0  k/n over the exact-tie reps (shipped).
  E1  q = sum_j w_j x_jf / sum_j w_j over the top K, w_j = exp(-(id_1 - id_j)/tau).
      K in {20, 50, 100}, tau in {0.0025, 0.005, 0.01, 0.02, 0.05}, dedup genome|species.
  E2  q = (k + alpha pi_c) / (n + alpha); pi_c = share of the fold's reference genomes in
      the anchor's clade carrying f; anchor = reps[0] (prFBA's own anchor order); clade
      start c in {genus, family}, back-off genus -> family -> order when the clade has
      < 3 reference genomes or its label is a panel placeholder ('unclassified_*');
      no admissible clade up to order -> alpha term dropped (q = k/n).
      alpha in {0.5, 1, 2, 4, 8}.
  E3  q = (sum_j w_j x_jf + alpha pi_c) / (sum_j w_j + alpha), K = 100,
      tau in {0.005, 0.01, 0.02}, both dedups, c in {genus, family}, alpha in {0.5,1,2,4}.
  E4  out-of-fold isotonic recalibration of each family's selected variant: for outer
      fold f, PAVA fitted on the four other folds' universe points of the variant
      selected for f (binned: exact 0, 4,000 uniform bins, exact 1; x = bin mean q,
      weight = points), applied to fold f's exact q by the fitted piecewise-linear map.

TUNING (nested)
  For each outer fold, each family's variant minimises the pooled universe log-loss
  (eps 0.01) of the other four folds; the outer fold is scored with it; reported numbers
  pool the outer folds. Exact ties: smaller K, larger tau, genome before species, genus
  before family, smaller alpha. Carried forward: the configuration chosen in most outer
  folds (same tie order).

METRICS
  universe (all 2,235 ECs): Brier, log-loss eps 0.01 (every point clipped to
    [eps, 1-eps]), ECE over deciles (calibrate's bins: exact 0 | 10 deciles | exact 1);
    OLS calibration slope (extra).
  union support = the exact-tie reps' union (E0's support) for EVERY estimator: Brier,
    WLS calibration slope (= point-level OLS of y on q), AUROC (exact for the estimate;
    1,001-bin histograms for the bootstrap; ties count 1/2).
  decision: macro and micro F1 of q >= 0.5 (over the universe) minus the union's,
    read against the SESOI 0.02 with score_folds' rule.
  strata: all; n = 1, 2, 3, >= 4; n >= 2 (multi).

UNCERTAINTY
  calibrate_pfunc.make_schemes / boot_weights verbatim (score_folds' draw:
  default_rng(seed).integers(0, K, size=(B, K)) over the stratum's sorted cluster keys),
  B = 2,000, seed 1729; primary unit holdout group, sensitivities source genome and
  genus. All estimators of a stratum share the draws, so contrasts vs E0 are paired.

DEVIATIONS / INTERPRETATIONS: see write_readme() (README section "Deviations").

USAGE
  PY=~/Documents/py_venv/bin/python
  $PY scripts/pfunc_estimators.py --workers 12 --bootstrap 2000 --seed 1729 \\
      --densities 1.0,0.5,0.25 --outdir data/pfunc_estimators --cache-dir <scratch>

Interpreter: /home/freiburger/Documents/py_venv/bin/python (the ONLY interpreter).
Reuses calibrate_pfunc's stateless helpers by import (never modifies it); the EC-universe,
truth-matrix and taxonomy-code block of calibrate_pfunc.load_context and the stage-0
ranking of calibrate_pfunc.stage0_ranked are copied verbatim (attributed inline).
"""
from __future__ import annotations

import os

# single-threaded BLAS per process: the worker count is the CPU budget (as calibrate_pfunc)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse            # noqa: E402
import collections         # noqa: E402
import hashlib             # noqa: E402
import math                # noqa: E402
import multiprocessing as mp  # noqa: E402
import pickle              # noqa: E402
import resource            # noqa: E402
import sys                 # noqa: E402
import time                # noqa: E402

import numpy as np         # noqa: E402
import orjson              # noqa: E402
import pandas as pd        # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)                                   # .../shared_reference
D = os.path.join(SR, "data")
REPO = os.path.dirname(os.path.dirname(SR))                  # .../prFBA
for _p in (HERE, REPO):
    if _p not in sys.path:
        sys.path.insert(0, _p)

import calibrate_pfunc as CP      # noqa: E402  stateless helpers only (CP.G is never used)
import select_references as SRF   # noqa: E402  the species-dedup key (_norm)

CELLS = os.path.join(D, "cells.csv.gz")
FOLD_PLAN = os.path.join(D, "fold_plan.json")
PANEL = os.path.join(D, "panel_eligible.json")
EC_UNI = os.path.join(D, "panel_ec_universe.json")
PANEL_EC = os.path.join(D, "_cache_panel_ec.json")
POPULATION = os.path.join(D, "weighted_objective", "population.json")
CALIB_SCORES = os.path.join(D, "calibration", "scores.json")
DESIGN = os.path.join(SR, "PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md")

DTAG = {"1.0": "d100", "0.5": "d050", "0.25": "d025"}
STAGE0 = CP.STAGE0                     # coverage 0.90, |SW-edlib| 0.02, floor 0.865
HITS_USED = 100
K1 = (20, 50, 100)
TAU1 = (0.0025, 0.005, 0.01, 0.02, 0.05)
DEDUPS = ("genome", "species")
CLADES = ("genus", "family")
ALPHA2 = (0.5, 1.0, 2.0, 4.0, 8.0)
TAU3 = (0.005, 0.01, 0.02)
ALPHA3 = (0.5, 1.0, 2.0, 4.0)
K3 = 100
REP_K = (5, 10, 20, 50, 100)           # 200-hit reproduction of calibrate's exploratory top-k
REP_TAU = 0.01
EPS = 0.01
NBF = 4000                             # isotonic fit bins (plus exact 0 and exact 1)
NBA = 1001                             # AUROC bootstrap bins: floor(1000 q) for q < 1; q == 1
BACKOFF = ("genus", "family", "order")
MIN_CLADE = 3
PLACEHOLDER = "unclassified_"
STRATA = ("all", "n1", "n2", "n3", "n4+", "multi")
SCHEMES = ("group", "src", "genus")
FAMILIES = ("E1", "E2", "E3")
TARGETS = {"E0_universe_brier": 0.08891442, "E0_union_macro_F1": 0.7672148}

G: dict = {}            # global context, populated in the parent BEFORE forking workers


def log(msg):
    print(f"[pfunc {time.strftime('%H:%M:%S')}] {msg}", flush=True)


jload = CP.jload


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ============================================================================ variants
def build_variants():
    """The fixed variant table. Row 0 = E0; then E1 (dedup, tau, K); E2 (clade, alpha);
    E3 (dedup, tau, clade, alpha); then the 200-hit reproduction rows (R200)."""
    V = [dict(name="E0", family="E0", K=None, tau=None, dedup=None, clade=None, alpha=None,
              hits=200)]
    ix = {"E1": {}, "E2": {}, "E3": {}, "R200": {}}
    for d in DEDUPS:
        for tau in TAU1:
            for K in K1:
                ix["E1"][(d, tau, K)] = len(V)
                V.append(dict(name=f"E1|K={K}|tau={tau}|{d}", family="E1", K=K, tau=tau,
                              dedup=d, clade=None, alpha=None, hits=HITS_USED))
    for c in CLADES:
        for a in ALPHA2:
            ix["E2"][(c, a)] = len(V)
            V.append(dict(name=f"E2|{c}|alpha={a:g}", family="E2", K=None, tau=None,
                          dedup=None, clade=c, alpha=a, hits=200))
    for d in DEDUPS:
        for tau in TAU3:
            for c in CLADES:
                for a in ALPHA3:
                    ix["E3"][(d, tau, c, a)] = len(V)
                    V.append(dict(name=f"E3|K={K3}|tau={tau}|{d}|{c}|alpha={a:g}",
                                  family="E3", K=K3, tau=tau, dedup=d, clade=c, alpha=a,
                                  hits=HITS_USED))
    for K in REP_K:
        ix["R200"][K] = len(V)
        V.append(dict(name=f"R200|K={K}|tau={REP_TAU}|genome", family="R200", K=K,
                      tau=REP_TAU, dedup="genome", clade=None, alpha=None, hits=200))
    return V, ix


def tie_key(v):
    """Exact-tie order (smaller K, larger tau, genome before species, genus before family,
    smaller alpha)."""
    return ((v["K"] if v["K"] is not None else 0),
            -(v["tau"] if v["tau"] is not None else 0.0),
            (0 if v["dedup"] in (None, "genome") else 1),
            (0 if v["clade"] in (None, "genus") else 1),
            (v["alpha"] if v["alpha"] is not None else 0.0))


# ============================================================================ context
def load_context():
    """Panel, taxonomy codes, holdout groups, EC universe and the truth matrix E.
    The panel / taxonomy / EC-universe / truth block is calibrate_pfunc.load_context
    (l.218-281) verbatim; the clade ranks (genus, family, order) are added for E2/E3."""
    t0 = time.time()
    A = collections.OrderedDict()
    panel = jload(PANEL)["genomes"]
    gids = sorted(panel)
    gidx = {g: i for i, g in enumerate(gids)}
    NG = len(gids)
    G.update(gids=gids, gidx=gidx, NG=NG)

    fp = jload(FOLD_PLAN)
    gg = fp["genome_group"]
    G["SESOI"] = float(fp["meta"]["frozen"]["sesoi_delta_f1"])
    A["fold_plan_sesoi_delta_f1"] = G["SESOI"]
    A["panel_genomes_without_fold_plan_group"] = sum(1 for g in gids if g not in gg)

    def codes(vals):
        u = {v: i for i, v in enumerate(sorted(set(vals), key=str))}
        return np.array([u[v] for v in vals], np.int64)
    G["group_name"] = [gg.get(g, f"nogroup:{g}") for g in gids]
    G["group_code"] = codes(G["group_name"])
    G["genus_code"] = codes([panel[g]["genus"] or f"nogenus:{g}" for g in gids])
    G["species_code"] = codes([panel[g]["species"] or f"nospecies:{g}" for g in gids])
    G["taxon_code"] = codes([str(panel[g]["taxon_id"]) for g in gids])
    A["panel_taxonomy"] = {
        "genomes": NG, "groups": int(G["group_code"].max() + 1),
        "genera": int(G["genus_code"].max() + 1),
        "species": int(G["species_code"].max() + 1),
        "taxon_ids": int(G["taxon_code"].max() + 1),
        "genera_spanning_more_than_one_group": int(sum(
            1 for n in collections.Counter(
                gen for gen, _ in set(zip(G["genus_code"].tolist(),
                                          G["group_code"].tolist()))).values() if n > 1))}
    del fp, gg

    # ---- clade ranks for E2/E3 (panel lineage; 'unclassified_*' labels are placeholders)
    rank_code, rank_labels, rank_ph = {}, {}, {}
    for rank in BACKOFF:
        labels = [panel[g][rank] or f"no{rank}:{g}" for g in gids]
        u = sorted(set(labels))
        ui = {v: i for i, v in enumerate(u)}
        rank_code[rank] = np.array([ui[v] for v in labels], np.int64)
        rank_labels[rank] = u
        rank_ph[rank] = np.array([v.startswith(PLACEHOLDER) or v.startswith(f"no{rank}:")
                                  for v in u], bool)
        A[f"clade_rank_{rank}"] = {
            "labels": len(u), "placeholder_labels": int(rank_ph[rank].sum()),
            "genomes_with_placeholder_label": int(rank_ph[rank][rank_code[rank]].sum())}
    G.update(rank_code=rank_code, rank_labels=rank_labels, rank_ph=rank_ph)
    G["family_code"] = rank_code["family"]
    G["order_code"] = rank_code["order"]

    # ---- EC universe + truth (score_folds: panel_ec_universe after norm_ec) -----------
    uni, seen = [], set()
    for e in jload(EC_UNI)["ec"]:
        n = CP.norm_ec(e)
        if n and n not in seen:
            seen.add(n)
            uni.append(n)
    ec_index = {e: i for i, e in enumerate(uni)}
    NU = len(uni)
    pec = jload(PANEL_EC)
    E = np.zeros((NG, NU), np.uint8)
    n_changed = 0
    for g, v in pec.items():
        i = gidx.get(g)
        if i is None:
            continue
        ids = []
        for e in v["ec"]:
            ne = CP.norm_ec(e)
            n_changed += int(ne != e)
            if ne is not None and ne in ec_index:
                ids.append(ec_index[ne])
        E[i, ids] = 1
    A["panel_ec_entries_changed_by_norm_ec"] = n_changed
    A["panel_genomes_without_panel_ec_record"] = sum(1 for g in gids if g not in pec)
    A["ec_universe_size"] = NU
    del pec
    G.update(ec_list=uni, ec_index=ec_index, NU=NU, E=E)
    G["A_ctx"] = A
    log(f"context: {NG} panel genomes, EC universe {NU} ({time.time()-t0:.1f}s, "
        f"rss {rss_gb():.1f} GB)")


def load_cells(densities, folds, limit):
    t0 = time.time()
    usecols = ["amplicon", "md5", "region", "src_genome_id", "fold", "group_id", "density",
               "arm", "space", "outcome", "n_truth", "f1", "n_reps", "tier",
               "rank1_identity", "novelty_bin"]
    parts = []
    for ch in pd.read_csv(CELLS, dtype=str, usecols=usecols, low_memory=False,
                          chunksize=500_000, keep_default_na=False):
        ch = ch[(ch.arm == "prfba_exact_tie") & (ch.space == "ec") &
                (ch.outcome == "scored") & ch.density.isin(densities)]
        parts.append(ch)
    df = pd.concat(parts, ignore_index=True)
    out = {}
    for d in densities:
        s = df[df.density == d].copy()
        s["fold_i"] = s.fold.astype(int)
        s = s[s.fold_i.isin(folds)]
        s = s.sort_values(["fold_i", "amplicon"]).reset_index(drop=True)
        if limit:
            s = s.groupby("fold_i", group_keys=False).head(limit).reset_index(drop=True)
        out[d] = s
        log(f"cells {DTAG[d]}: {len(s)} scored exact_tie EC cells "
            f"({s.src_genome_id.nunique()} sources, {s.group_id.nunique()} groups)")
    log(f"cells loaded ({time.time()-t0:.1f}s)")
    return out


# ============================================================================ streaming
def stage0_lists(top, L):
    """select_references stage 0 + identity ranking + genome dedup: calibrate_pfunc.
    stage0_ranked (l.1585-1609) verbatim, plus each survivor's BV-BRC species key exactly
    as select_references' select_all_species_dedup builds it."""
    surv = []
    for h in top:
        ident = h["identity"]
        cov = round(min(h["aligned_len"] / L, 1.0), 4)
        ei = h.get("edlib_identity")
        agree = None if ei is None else round(abs(ident - ei), 4)
        if cov < STAGE0["cov_min"]:
            continue
        if agree is not None and agree > STAGE0["swed_max"]:
            continue
        if ident < STAGE0["family_floor"]:
            continue
        surv.append(h)
    surv.sort(key=lambda h: (-h["identity"], -h["align_score"]))
    seen, ranked = set(), []
    for h in surv:
        g = str(h["genome_id"])
        if g in seen:
            continue
        seen.add(g)
        lin = h.get("lineage") or {}
        sp = SRF._norm(lin.get("Species")) or (f"taxon:{h['taxon_id']}"
                                               if h.get("taxon_id") else g)
        ranked.append((g, float(h["identity"]), sp))
    return ranked


def stream_job(job):
    """One (density, fold): fold reference + EC cache check, reps from predictions.json,
    and the stage-0 neighbour lists of every wanted md5 from the hits file."""
    dtag, fold, md5s = job
    import ijson
    t0 = time.time()
    gidx, E, ec_index = G["gidx"], G["E"], G["ec_index"]
    ddir = os.path.join(D, "prfba", f"f{fold:02d}_{dtag}")
    want = set(md5s)
    A = collections.Counter()
    # (1) the fold reference = _cache_ec_subset keys; its EC sets == panel truth rows
    ecsub = jload(os.path.join(ddir, "_cache_ec_subset.json"))
    ref = []
    for g, v in ecsub.items():
        i = gidx.get(g)
        if i is None:
            A["ec_subset_genome_outside_panel"] += 1
            continue
        ref.append(i)
        a = sorted(ec_index.get(e, -1) for e in v)
        if a != np.flatnonzero(E[i]).tolist():
            A["ec_subset_mismatch_vs_panel_truth"] += 1
    ref = np.array(sorted(ref), np.int64)
    del ecsub
    # (2) reps (prFBA's order: reps[0] is the anchor)
    reps = {}
    with open(os.path.join(ddir, "predictions.json"), "rb") as fh:
        for m, rec in ijson.kvitems(fh, "pred.exact_tie", use_float=True):
            if m in want:
                reps[m] = (list(rec["reps"]), int(rec["n_reps"]),
                           np.array(sorted(ec_index.get(e, -1) for e in rec["ec"]), np.int32))
    t1 = time.time()
    # (3) hits: first 100 hits -> stage 0 -> ranked -> genome / species dedup
    names = ("g100", "s100", "g200")
    lists = {k: ([], [], [0]) for k in names}           # gi, id, ptr
    order, info = [], collections.defaultdict(list)
    path = os.path.join(ddir, "hits", "asv_top20_alignment_hits.json")
    with open(path, "rb") as fh:
        for m, rec in ijson.kvitems(fh, "", use_float=True):
            if m not in want or not isinstance(rec, dict):
                continue
            top = rec.get("top20") or []
            if [h.get("rank") for h in top] != list(range(1, len(top) + 1)):
                A["hits_rank_order_violation"] += 1
            L = rec.get("asv_len") or 1
            r100 = stage0_lists(top[:HITS_USED], L)
            r200 = stage0_lists(top, L)[:max(REP_K)]
            sp_seen, s100 = set(), []
            for g, x, sp in r100:
                if sp in sp_seen:
                    continue
                sp_seen.add(sp)
                s100.append((g, x, sp))
            for k, rl in (("g100", r100), ("s100", s100), ("g200", r200)):
                gi, idv, ptr = lists[k]
                for g, x, _sp in rl:
                    i = gidx.get(g)
                    if i is None:
                        A[f"hit_genome_outside_panel_{k}"] += 1
                        continue
                    gi.append(i)
                    idv.append(x)
                ptr.append(len(gi))
            Lq = max(1, int(rec.get("asv_len") or 1))
            qual = [h for h in top[:HITS_USED]
                    if (h.get("aligned_len") or 0) / Lq >= STAGE0["cov_min"]]
            last_q = float(qual[-1]["identity"]) if qual else None
            order.append(m)
            info["n_hits"].append(len(top))
            info["saturated100"].append(bool(len(top) >= HITS_USED and qual and
                                             (last_q or 0.0) >= STAGE0["family_floor"]))
            info["species_collapsed"].append(len(r100) - len(s100))
    out = {"dtag": dtag, "fold": fold, "ref": ref, "reps": reps, "order": order,
           "info": {k: np.asarray(v) for k, v in info.items()},
           "lists": {k: (np.asarray(v[0], np.int32), np.asarray(v[1], np.float64),
                         np.asarray(v[2], np.int64)) for k, v in lists.items()},
           "asserts": dict(A), "seconds_predictions": round(t1 - t0, 1),
           "seconds": round(time.time() - t0, 1), "peak_rss_gb": round(rss_gb(), 2)}
    return out


# ============================================================================ dataset
def build_dataset(dkey, cells, parts):
    """Merge cells, reps and neighbour lists; build the clade priors; run every join and
    leakage check. Returns the dataset dict (stored in G['DS'][dkey])."""
    t0 = time.time()
    gidx, E, NU = G["gidx"], G["E"], G["NU"]
    A = collections.Counter()
    parts = sorted(parts, key=lambda p: p["fold"])
    folds = [p["fold"] for p in parts]
    ref = {p["fold"]: p["ref"] for p in parts}
    in_ref = {}
    for f, r in ref.items():
        m = np.zeros(G["NG"], bool)
        m[r] = True
        in_ref[f] = m
    # ---- global record table: (fold, md5) -> record id; CSR neighbour lists (parts in
    # fold order; records of a part in its stream order)
    rec_of = {}
    rinfo = collections.defaultdict(list)
    for p in parts:
        base = len(rec_of)
        for j, m in enumerate(p["order"]):
            rec_of[(p["fold"], m)] = base + j
        for k, v in p["info"].items():
            rinfo[k].append(v)
        for k, v in p["asserts"].items():
            A[f"stream_{k}"] += v
    lists = {}
    for k in ("g100", "s100", "g200"):
        gis, ids_, ptrs, off = [], [], [np.zeros(1, np.int64)], 0
        for p in parts:
            gi, idv, ptr = p["lists"][k]
            gis.append(gi)
            ids_.append(idv)
            ptrs.append(ptr[1:] + off)
            off += len(gi)
        lists[k] = (np.concatenate(gis).astype(np.int64), np.concatenate(ids_),
                    np.concatenate(ptrs))
    for k in lists:
        gi, idv, ptr = lists[k]
        if len(ptr) != len(rec_of) + 1 or ptr[-1] != len(gi):
            raise RuntimeError(f"CSR {k}: ptr {len(ptr)} vs records {len(rec_of)}")
    rinfo = {k: np.concatenate(v) for k, v in rinfo.items()}
    nrec = len(rec_of)
    same_sp = np.zeros(nrec, bool)
    g_ptr, s_ptr = lists["g100"][2], lists["s100"][2]
    for r in range(nrec):
        a = lists["g100"][0][g_ptr[r]:g_ptr[r + 1]]
        b = lists["s100"][0][s_ptr[r]:s_ptr[r + 1]]
        same_sp[r] = (len(a) == len(b)) and bool(np.all(a == b))
    reps_of = {}
    for p in parts:
        for m, v in p["reps"].items():
            reps_of[(p["fold"], m)] = v

    # ---- per-cell arrays
    n_cells = len(cells)
    fold = cells.fold_i.to_numpy(np.int64)
    src = np.array([gidx.get(s, -1) for s in cells.src_genome_id], np.int64)
    rec = np.full(n_cells, -1, np.int64)
    reps_flat, reps_ptr = [], [0]
    keep = np.ones(n_cells, bool)
    for i, row in enumerate(cells.itertuples(index=False)):
        key = (int(row.fold_i), row.md5)
        rv = reps_of.get(key)
        r_id = rec_of.get(key)
        if rv is None:
            A["md5_missing_in_predictions"] += 1
            keep[i] = False
        if r_id is None:
            A["hits_record_missing"] += 1
            keep[i] = False
        if src[i] < 0:
            A["source_outside_panel"] += 1
            keep[i] = False
        reps = rv[0] if rv is not None else []
        ri = [gidx.get(r, -1) for r in reps]
        if not ri or -1 in ri:
            A["reps_empty_or_outside_panel"] += 1
            keep[i] = False
            ri = []
        if rv is not None and (len(reps) != int(row.n_reps) or len(reps) != rv[1]):
            A["n_reps_mismatch"] += 1
        rec[i] = r_id if r_id is not None else -1
        reps_flat.extend(ri)
        reps_ptr.append(len(reps_flat))
    if not keep.all():
        raise RuntimeError(f"{dkey}: {int((~keep).sum())} cells failed the joins: {dict(A)}")
    reps_flat = np.asarray(reps_flat, np.int64)
    reps_ptr = np.asarray(reps_ptr, np.int64)
    n = np.diff(reps_ptr)
    anchor = reps_flat[reps_ptr[:-1]]
    nb = np.minimum(n, 4) - 1                                   # 0:n1 1:n2 2:n3 3:n4+

    # ---- joins against cells.csv.gz and predictions.json (as calibrate_pfunc layer A)
    T = np.zeros(n_cells)
    Ns = np.zeros(n_cells)
    POSs = np.zeros(n_cells)
    F1U = np.zeros(n_cells)
    for i, row in enumerate(cells.itertuples(index=False)):
        r = reps_flat[reps_ptr[i]:reps_ptr[i + 1]]
        k = E[r].sum(0, dtype=np.int16)
        y = E[src[i]].astype(bool)
        u = k > 0
        T[i] = y.sum()
        Ns[i] = u.sum()
        POSs[i] = (u & y).sum()
        F1U[i] = 2 * POSs[i] / (Ns[i] + T[i]) if Ns[i] + T[i] else 0.0
        if not np.array_equal(np.flatnonzero(u), reps_of[(int(row.fold_i), row.md5)][2]):
            A["union_ec_ne_pred_ec"] += 1
        if abs(F1U[i] - float(row.f1)) > 1e-3:
            A["f1_mismatch_vs_cells"] += 1
        if int(T[i]) != int(row.n_truth):
            A["n_truth_mismatch_vs_cells"] += 1
        if row.group_id != G["group_name"][src[i]]:
            A["group_id_ne_fold_plan_genome_group"] += 1
        A["reps_outside_fold_reference"] += int((~in_ref[int(row.fold_i)][r]).sum())

    # ---- clades: per fold, per rank, label counts over the fold reference
    rc, rph = G["rank_code"], G["rank_ph"]
    cnt = {(f, rk): np.bincount(rc[rk][ref[f]], minlength=len(G["rank_labels"][rk]))
           for f in folds for rk in BACKOFF}
    clade_id = {}
    crow = np.full((len(CLADES), n_cells), -1, np.int64)
    crank = np.full((len(CLADES), n_cells), -1, np.int64)
    for i in range(n_cells):
        f = int(fold[i])
        a = anchor[i]
        for ci, start in enumerate(CLADES):
            for ri in range(BACKOFF.index(start), len(BACKOFF)):
                rk = BACKOFF[ri]
                code = int(rc[rk][a])
                if rph[rk][code] or cnt[(f, rk)][code] < MIN_CLADE:
                    continue
                key = (f, rk, code)
                if key not in clade_id:
                    clade_id[key] = len(clade_id)
                crow[ci, i] = clade_id[key]
                crank[ci, i] = ri
                break
    nclade = len(clade_id)
    PI = np.zeros((nclade, NU), np.float64)
    cm_flat, cm_ptr = [], [0]
    cmeta = [None] * nclade
    for (f, rk, code), j in sorted(clade_id.items(), key=lambda x: x[1]):
        mem = ref[f][rc[rk][ref[f]] == code]
        PI[j] = E[mem].mean(0)
        cm_flat.extend(mem.tolist())
        cm_ptr.append(len(cm_flat))
        cmeta[j] = (f, rk, G["rank_labels"][rk][code], len(mem))
    cm_flat = np.asarray(cm_flat, np.int64)
    cm_ptr = np.asarray(cm_ptr, np.int64)
    for ci, start in enumerate(CLADES):
        lab = collections.Counter(int(x) for x in crank[ci])
        A[f"clade_start_{start}_cells_at_rank"] = {
            ("none" if k < 0 else BACKOFF[k]): v for k, v in sorted(lab.items())}
    # why no admissible clade (genus start): the anchor's order label is a placeholder, or
    # a real order with < MIN_CLADE genomes in the fold reference
    why = collections.Counter()
    for i in np.flatnonzero(crow[0] < 0):
        code = int(rc["order"][anchor[i]])
        why["order_placeholder" if rph["order"][code] else
            f"order_lt{MIN_CLADE}_reference_genomes"] += 1
    A["clade_none_reason_at_order"] = dict(why)

    # ---- leakage: neighbours (all three lists) and clade members vs the source
    codes = {nm: G[f"{nm}_code"] for nm in ("genus", "species", "taxon", "group")}
    diag = collections.Counter()
    for i in range(n_cells):
        s = src[i]
        f = int(fold[i])
        a = anchor[i]
        A["leak_anchor_genus_eq_source_genus"] += int(codes["genus"][a] == codes["genus"][s])
        r = rec[i]
        for k in ("g100", "s100", "g200"):
            gi, _idv, ptr = lists[k]
            g = gi[ptr[r]:ptr[r + 1]]
            A[f"leak_source_in_neighbours_{k}"] += int((g == s).any())
            for nm, arr in codes.items():
                A[f"leak_neighbour_shares_{nm}_{k}"] += int((arr[g] == arr[s]).any())
            A[f"neighbours_outside_fold_reference_{k}"] += int((~in_ref[f][g]).sum())
        for ci, start in enumerate(CLADES):
            j = crow[ci, i]
            if j < 0:
                continue
            mem = cm_flat[cm_ptr[j]:cm_ptr[j + 1]]
            if cmeta[j][0] != f:
                A["clade_fold_ne_cell_fold"] += 1
            A[f"leak_source_in_clade_{start}"] += int((mem == s).any())
            for nm, arr in codes.items():
                A[f"leak_clade_member_shares_{nm}_{start}"] += int((arr[mem] == arr[s]).any())
            A[f"clade_members_outside_fold_reference_{start}"] += int((~in_ref[f][mem]).sum())
        # diagnostics: 100-hit truncation vs the 200-hit exact-tie reps
        gi, idv, ptr = lists["g100"]
        g = gi[ptr[r]:ptr[r + 1]]
        x = idv[ptr[r]:ptr[r + 1]]
        reps = set(reps_flat[reps_ptr[i]:reps_ptr[i + 1]].tolist())
        diag["cells"] += 1
        if len(g) == 0:
            diag["no_survivors_in_first_100_hits"] += 1
            continue
        tie = set(g[x == x[0]].tolist())
        diag["tie_set_100_ne_reps"] += int(tie != reps)
        diag["rep_missing_from_100hit_survivors"] += int(not reps <= set(g.tolist()))
        for K in K1:
            diag[f"fewer_than_K{K}_survivors_100"] += int(len(g) < K)
        gi2, idv2, ptr2 = lists["g200"]
        g2 = gi2[ptr2[r]:ptr2[r + 1]]
        x2 = idv2[ptr2[r]:ptr2[r + 1]]
        if len(g2):
            diag["tie_set_200_ne_reps"] += int(set(g2[x2 == x2[0]].tolist()) != reps)
        diag["survivors_100_sum"] += len(g)
    nsurv = np.diff(lists["g100"][2])
    diag["records"] = nrec
    diag["records_species_dedup_differs_from_genome_dedup"] = int((~same_sp).sum())
    diag["records_saturated_at_100_hits"] = int(rinfo["saturated100"].sum())
    diag["records_with_200_hits"] = int((rinfo["n_hits"] >= 200).sum())
    ds = dict(key=dkey, n_cells=n_cells, folds=folds, cells=cells, fold=fold, src=src,
              rec=rec, reps_flat=reps_flat, reps_ptr=reps_ptr, n=n, nb=nb, anchor=anchor,
              lists=lists, same_sp=same_sp, rinfo=rinfo, ref=ref,
              crow=crow, crank=crank, PI=PI, cm_flat=cm_flat, cm_ptr=cm_ptr, cmeta=cmeta,
              T=T, Ns=Ns, POSs=POSs, F1U=F1U, asserts=dict(A), diag=dict(diag),
              survivors_100_quantiles={str(q): float(np.quantile(nsurv, q))
                                       for q in (0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0)})
    log(f"dataset {dkey}: {n_cells} cells, {nrec} records, {nclade} clades "
        f"({time.time()-t0:.1f}s, rss {rss_gb():.1f} GB); diag {dict(diag)}")
    return ds


# ============================================================================ estimators
def csr_slice(lst, r):
    gi, idv, ptr = lst
    return gi[ptr[r]:ptr[r + 1]], idv[ptr[r]:ptr[r + 1]]


def compute_Q(ds, i):
    """Every variant's q over the universe for cell i: (NV, NU)."""
    E, NU, V, ix = G["E"], G["NU"], G["V"], G["VIX"]
    NV = len(V)
    s = ds["src"][i]
    y = E[s].astype(bool)
    r = ds["reps_flat"][ds["reps_ptr"][i]:ds["reps_ptr"][i + 1]]
    n = len(r)
    k = E[r[0]].astype(np.float64) if n == 1 else E[r].sum(0, dtype=np.float64)
    Q = np.empty((NV, NU))
    Q[0] = k / n
    rid = ds["rec"][i]
    fb = 0
    AK = {}
    nt, nk = len(TAU1), len(K1)
    for di, dname in enumerate(DEDUPS):
        rows = [ix["E1"][(dname, tau, K)] for tau in TAU1 for K in K1]
        if di == 1 and ds["same_sp"][rid]:
            Q[rows] = Q[G["E1_rows"][0]]
            for t3 in TAU3:
                AK[(dname, t3)] = AK[(DEDUPS[0], t3)]
            continue
        g, ids = csr_slice(ds["lists"]["g100" if di == 0 else "s100"], rid)
        if len(g) == 0:                 # no stage-0 survivor in the first 100 hits
            g, ids = r, np.ones(n)
            fb = 1
        X = E[g].astype(np.float64)
        dif = ids[0] - ids
        W = np.zeros((nt * nk, len(g)))
        for ti, tau in enumerate(TAU1):
            w = np.exp(-dif / tau)
            for ki, K in enumerate(K1):
                W[ti * nk + ki, :K] = w[:K]
        Am = W @ X
        den = W.sum(1)
        Q[rows] = Am / den[:, None]
        for t3 in TAU3:
            j = TAU1.index(t3) * nk + K1.index(K3)
            AK[(dname, t3)] = (Am[j], den[j])
    # 200-hit reproduction rows (calibrate_pfunc.topk_job: tau 0.01, genome dedup)
    g, ids = csr_slice(ds["lists"]["g200"], rid)
    if len(g):
        X = E[g].astype(np.float64)
        w = np.exp(-(ids[0] - ids) / REP_TAU)
        for K in REP_K:
            kk = min(K, len(g))
            Q[ix["R200"][K]] = (w[:kk] @ X[:kk]) / w[:kk].sum()
    else:
        for K in REP_K:
            Q[ix["R200"][K]] = 0.0
    # clade priors
    for ci, c in enumerate(CLADES):
        row = ds["crow"][ci, i]
        if row >= 0:
            pi, has = ds["PI"][row], 1.0
        else:
            pi, has = np.zeros(NU), 0.0
        for a in ALPHA2:
            Q[ix["E2"][(c, a)]] = (k + a * pi) / (n + a * has)
        for dname in DEDUPS:
            for t3 in TAU3:
                Am, den = AK[(dname, t3)]
                for a in ALPHA3:
                    Q[ix["E3"][(dname, t3, c, a)]] = (Am + a * pi) / (den + a * has)
    Q[Q >= 1 - 1e-9] = 1.0
    Q[Q <= 1e-12] = 0.0
    return Q, y, k > 0, n, fb


def basic_stats(Q, y, u, T):
    """Per-row additive statistics (pass A): SQ_u, LL_u, SQ_s, Sq_s, Sqq_s, Sqy_s, TP,
    PRED, F1, NSUP."""
    yf = y.astype(np.float64)
    Dm = Q - yf
    SQu = np.einsum("ij,ij->i", Dm, Dm)
    Qc = np.clip(Q, EPS, 1 - EPS)
    LLu = -np.log(np.where(y, Qc, 1 - Qc)).sum(1)
    Qs = Q[:, u]
    ys = y[u]
    Ds = Qs - ys
    SQs = np.einsum("ij,ij->i", Ds, Ds)
    Sqs = Qs.sum(1)
    Sqqs = np.einsum("ij,ij->i", Qs, Qs)
    Sqys = Qs[:, ys].sum(1)
    P = Q >= 0.5
    TP = P[:, y].sum(1)
    PRED = P.sum(1)
    F1 = np.where(PRED + T > 0, 2 * TP / np.maximum(PRED + T, 1), 0.0)
    NSUP = (Q > 0).sum(1)
    return np.stack([SQu, LLu, SQs, Sqs, Sqqs, Sqys, TP, PRED, F1, NSUP], 1)


NSA = 10
SA = dict(SQu=0, LLu=1, SQs=2, Sqs=3, Sqqs=4, Sqys=5, TP=6, PRED=7, F1=8, NSUP=9)


def fine_bins(Q):
    b = np.minimum((Q * NBF).astype(np.int64), NBF - 1) + 1
    b[Q <= 0] = 0
    b[Q >= 1] = NBF + 1
    return b


def passA_job(job):
    dkey, lo, hi = job
    t0 = time.time()
    ds = G["DS"][dkey]
    NV = len(G["V"])
    NB = NBF + 2
    st = np.zeros((hi - lo, NV, NSA))
    H = np.zeros((3, NV * NB))
    roff = (np.arange(NV) * NB)[:, None]
    buf_k, buf_y, buf_q = [], [], []
    nfb = 0

    def flush():
        if not buf_k:
            return
        kk = np.concatenate(buf_k)
        yy = np.concatenate(buf_y)
        qq = np.concatenate(buf_q)
        H[0] += np.bincount(kk, minlength=NV * NB)
        H[1] += np.bincount(kk[yy], minlength=NV * NB)
        H[2] += np.bincount(kk, weights=qq, minlength=NV * NB)
        buf_k.clear()
        buf_y.clear()
        buf_q.clear()

    for j, i in enumerate(range(lo, hi)):
        Q, y, u, n, fb = compute_Q(ds, i)
        nfb += fb
        st[j] = basic_stats(Q, y, u, ds["T"][i])
        buf_k.append((fine_bins(Q) + roff).ravel())
        buf_y.append(np.broadcast_to(y, Q.shape).ravel())
        buf_q.append(Q.ravel())
        if len(buf_k) >= 24:
            flush()
    flush()
    return {"dkey": dkey, "lo": lo, "hi": hi, "fold": int(ds["fold"][lo]), "stats": st,
            "H": H.reshape(3, NV, NB), "fallback_cells": nfb,
            "seconds": round(time.time() - t0, 1), "peak_rss_gb": round(rss_gb(), 2)}


# ---- pass B: the reported composites -------------------------------------------------
NSB = 49
SB_ = dict(SQu=0, LLu=1, Squ=2, Sqqu=3, Sqyu=4, SQs=5, Sqs=6, Sqqs=7, Sqys=8, TP=9,
           PRED=10, F1=11, NSUP=12)


def ece_bins(q):
    """calibrate's decile convention (topk_job), with exact 0 and exact 1 as own bins:
    0 = q == 0; 1..10 = (k/10, (k+1)/10]; 11 = q == 1."""
    b = np.clip(np.ceil(q * 10 - 1e-12).astype(np.int64) - 1, 0, 9) + 1
    b[q <= 0] = 0
    b[q >= 1] = 11
    return b


def full_stats(q, y, u, T):
    yf = y.astype(np.float64)
    d = q - yf
    qc = np.clip(q, EPS, 1 - EPS)
    qs = q[u]
    ys = y[u]
    ds_ = qs - ys
    p = q >= 0.5
    TP = float((p & y).sum())
    PRED = float(p.sum())
    out = np.zeros(NSB)
    out[0] = d @ d
    out[1] = -np.log(np.where(y, qc, 1 - qc)).sum()
    out[2] = q.sum()
    out[3] = q @ q
    out[4] = q[y].sum()
    out[5] = ds_ @ ds_
    out[6] = qs.sum()
    out[7] = qs @ qs
    out[8] = qs[ys].sum()
    out[9] = TP
    out[10] = PRED
    out[11] = 2 * TP / (PRED + T) if PRED + T > 0 else 0.0
    out[12] = (q > 0).sum()
    b = ece_bins(q)
    out[13:25] = np.bincount(b, minlength=12)
    out[25:37] = np.bincount(b[y], minlength=12)
    out[37:49] = np.bincount(b, weights=q, minlength=12)
    return out


def apply_map(q, mp_):
    if mp_ is None:
        return q
    r = np.interp(q, mp_[0], mp_[1])
    r[r >= 1 - 1e-9] = 1.0
    r[r <= 1e-12] = 0.0
    return r


def passB_job(job):
    dkey, lo, hi = job
    t0 = time.time()
    ds = G["DS"][dkey]
    comps = G["COMPS"][dkey]
    NC = len(comps)
    st = np.zeros((hi - lo, NC, NSB))
    qbuf = [[] for _ in range(NC)]
    ybuf, cbuf = [], []
    for j, i in enumerate(range(lo, hi)):
        Q, y, u, n, _fb = compute_Q(ds, i)
        f = int(ds["fold"][i])
        T = ds["T"][i]
        for c, comp in enumerate(comps):
            q = apply_map(Q[comp["rows"][f]], comp["maps"][f] if comp["maps"] else None)
            st[j, c] = full_stats(q, y, u, T)
            qbuf[c].append(q[u].astype(np.float32))
        ybuf.append(y[u])
        cbuf.append(np.full(int(u.sum()), i, np.int32))
    return {"dkey": dkey, "lo": lo, "hi": hi, "stats": st,
            "uq": np.stack([np.concatenate(x) for x in qbuf], 0),
            "uy": np.concatenate(ybuf), "ucell": np.concatenate(cbuf),
            "seconds": round(time.time() - t0, 1), "peak_rss_gb": round(rss_gb(), 2)}


# ============================================================================ selection
def fold_sums(dsA, ds):
    """Per (variant, fold) sums of the pass-A log-loss and the cell counts."""
    folds = ds["folds"]
    LL = np.zeros((len(G["V"]), len(folds)))
    NCf = np.zeros(len(folds))
    for fi, f in enumerate(folds):
        m = ds["fold"] == f
        LL[:, fi] = dsA["stats"][m][:, :, SA["LLu"]].sum(0)
        NCf[fi] = m.sum()
    return LL, NCf


def nested_select(dsA, ds):
    V, ix = G["V"], G["VIX"]
    NU = G["NU"]
    LL, NCf = fold_sums(dsA, ds)
    folds = ds["folds"]
    sel = {}
    for fam in FAMILIES:
        rows = sorted(ix[fam].values())
        per = {}
        for fi, f in enumerate(folds):
            tr = [x for x in range(len(folds)) if x != fi]
            ntr = NCf[tr].sum()
            ll = LL[rows][:, tr].sum(1) / (ntr * NU)
            best = ll.min()
            cand = [rows[j] for j in range(len(rows)) if ll[j] <= best + 1e-12 * abs(best)]
            choice = min(cand, key=lambda v: tie_key(V[v]))
            order = np.argsort(ll)
            ll_out = (LL[choice, fi] / (NCf[fi] * NU))
            per[f] = {"row": int(choice), "variant": V[choice]["name"],
                      "train_universe_logloss": float(ll[rows.index(choice)]),
                      "outer_universe_logloss": float(ll_out),
                      "exact_ties_at_minimum": [V[c]["name"] for c in cand],
                      "runners_up": [{"variant": V[rows[j]]["name"],
                                      "train_universe_logloss": float(ll[j])}
                                     for j in order[:6]]}
        cnt = collections.Counter(per[f]["row"] for f in folds)
        top = max(cnt.values())
        modal = min([v for v, c in cnt.items() if c == top], key=lambda v: tie_key(V[v]))
        # per-hyperparameter modes, reported only (the configuration mode is carried)
        hp_modes = {}
        for hp in ("K", "tau", "dedup", "clade", "alpha"):
            vals = collections.Counter(V[per[f]["row"]][hp] for f in folds)
            if any(v is not None for v in vals):
                hp_modes[hp] = vals.most_common()
        grids = {"E1": {"K": K1, "tau": TAU1}, "E2": {"alpha": ALPHA2},
                 "E3": {"tau": TAU3, "alpha": ALPHA3}}[fam]
        edge = {}
        for hp, grid in grids.items():
            vals = {V[per[f]["row"]][hp] for f in folds}
            edge[hp] = {"grid": list(grid), "chosen": sorted(vals),
                        "carried_forward_at": ("max" if V[modal][hp] == max(grid) else
                                               "min" if V[modal][hp] == min(grid) else
                                               "interior")}
        sel[fam] = {"per_fold": per,
                    "chosen_counts": {V[v]["name"]: c for v, c in cnt.most_common()},
                    "carried_forward_row": int(modal),
                    "carried_forward": {k: V[modal][k] for k in
                                        ("name", "K", "tau", "dedup", "clade", "alpha")},
                    "per_hyperparameter_modes": hp_modes,
                    "grid_edges": edge}
    return sel, LL, NCf


def fit_iso(Hsum):
    """PAVA on binned universe points: Hsum = (3, NB) [points, positives, sum q]."""
    from sklearn.isotonic import IsotonicRegression
    M, P, S = Hsum
    k = M > 0
    x = S[k] / M[k]
    yv = P[k] / M[k]
    ir = IsotonicRegression(y_min=0.0, y_max=1.0, increasing=True, out_of_bounds="clip")
    ir.fit(x, yv, sample_weight=M[k])
    return (np.asarray(ir.X_thresholds_, float), np.asarray(ir.y_thresholds_, float))


def make_composites(dkey, ds, dsA, sel, fixed=None):
    """d100: E0, the nested E1/E2/E3 composites and E4 of each (OOF isotonic: the map for
    outer fold f is fitted on the other folds' points of the variant selected for f).
    fixed: {label: row} -> one fixed variant for every fold (density sensitivity), plus
    E4 of the labels listed in fixed['_e4']."""
    folds = ds["folds"]
    Hf = dsA["H"]                                   # {fold: (3, NV, NB)}
    comps = [{"name": "E0", "base": "E0", "rows": {f: 0 for f in folds}, "maps": None}]
    if fixed is None:
        for fam in FAMILIES:
            comps.append({"name": fam, "base": fam,
                          "rows": {f: sel[fam]["per_fold"][f]["row"] for f in folds},
                          "maps": None})
    else:
        for lab, row in fixed.items():
            if lab.startswith("_"):
                continue
            comps.append({"name": lab, "base": lab, "rows": {f: row for f in folds},
                          "maps": None})
    e4_of = ([c["name"] for c in comps] if fixed is None else
             ["E0"] + list(fixed.get("_e4", [])))
    maps_log = {}
    for c in list(comps):
        if c["name"] not in e4_of:
            continue
        maps = {}
        for f in folds:
            v = c["rows"][f]
            Hs = sum(Hf[f2][:, v, :] for f2 in folds if f2 != f)
            maps[f] = fit_iso(Hs)
        comps.append({"name": f"E4({c['name']})", "base": c["name"], "rows": c["rows"],
                      "maps": maps})
        maps_log[f"E4({c['name']})"] = {str(f): {"x": maps[f][0].tolist(),
                                                 "y": maps[f][1].tolist()} for f in folds}
    return comps, maps_log


# ============================================================================ analysis
def strata_masks(ds):
    nb = ds["nb"]
    return {"all": np.ones(ds["n_cells"], bool), "n1": nb == 0, "n2": nb == 1,
            "n3": nb == 2, "n4+": nb == 3, "multi": nb >= 1}


def auroc_from_hist(Hp, Hn):
    """AUROC from per-bin positive/negative counts (last axis = bins, ascending q);
    ties within a bin count 1/2."""
    above = np.cumsum(Hp[..., ::-1], -1)[..., ::-1] - Hp
    num = (Hn * (above + 0.5 * Hp)).sum(-1)
    with np.errstate(all="ignore"):
        return num / (Hp.sum(-1) * Hn.sum(-1))


def auroc_exact(q, y):
    if len(q) == 0 or y.all() or not y.any():
        return float("nan")
    uq, inv = np.unique(q, return_inverse=True)
    Hp = np.bincount(inv, weights=y.astype(np.float64), minlength=len(uq))
    Hn = np.bincount(inv, weights=(~y).astype(np.float64), minlength=len(uq))
    return float(auroc_from_hist(Hp, Hn))


def derive_comp(S, L, names, NU):
    g = lambda k: S[..., L[k]]  # noqa: E731
    o = {}
    NC, T, Ns, POS, F1U = g("NCELL"), g("T"), g("N_s"), g("POS_s"), g("F1U")
    Nv = NC * NU
    with np.errstate(all="ignore"):
        uma = F1U / NC
        umi = 2 * POS / (Ns + T)
        o["union_rule"] = {"macro_F1": uma, "micro_F1": umi, "support_mean": Ns / NC,
                           "precision_on_union": POS / Ns}
        for c in names:
            b = g(f"{c}[]")
            v = {}
            v["universe_brier"] = b[..., 0] / Nv
            v["universe_logloss_eps0.01"] = b[..., 1] / Nv
            v["universe_ece_decile"] = np.abs(b[..., 37:49] - b[..., 25:37]).sum(-1) / Nv
            sxx = b[..., 3] - b[..., 2] ** 2 / Nv
            sxy = b[..., 4] - b[..., 2] * T / Nv
            v["universe_slope_ols"] = sxy / sxx
            v["union_brier"] = b[..., 5] / Ns
            sxx = b[..., 7] - b[..., 6] ** 2 / Ns
            sxy = b[..., 8] - b[..., 6] * POS / Ns
            v["union_slope_wls"] = sxy / sxx
            v["union_intercept_wls"] = POS / Ns - v["union_slope_wls"] * b[..., 6] / Ns
            v["union_mean_q"] = b[..., 6] / Ns
            v["macro_F1_q_ge_0.5"] = b[..., 11] / NC
            v["micro_F1_q_ge_0.5"] = 2 * b[..., 9] / (b[..., 10] + T)
            v["micro_P_q_ge_0.5"] = b[..., 9] / b[..., 10]
            v["micro_R_q_ge_0.5"] = b[..., 9] / T
            v["d_macro_F1_vs_union"] = v["macro_F1_q_ge_0.5"] - uma
            v["d_micro_F1_vs_union"] = v["micro_F1_q_ge_0.5"] - umi
            v["support_mean"] = b[..., 12] / NC
            v["pred_mean_q_ge_0.5"] = b[..., 10] / NC
            o[c] = v
    return o


CONTRAST_METRICS = ("universe_brier", "universe_logloss_eps0.01", "universe_ece_decile",
                    "union_brier", "union_slope_wls", "union_auroc", "macro_F1_q_ge_0.5",
                    "micro_F1_q_ge_0.5", "d_macro_F1_vs_union", "d_micro_F1_vs_union",
                    "universe_slope_ols")
# extra paired contrasts (in addition to every composite vs E0 and every E4(X) vs X)
PAIRS = (("E3", "E1"), ("E3", "E2"), ("E1", "E2"), ("E4(E3)", "E1"),
         ("E3*", "E1*"), ("E3*", "E2*"), ("E1*", "E2*"))


def analyze_job(job):
    """One (dataset, stratum): per-source sums, bootstrap schemes, every composite's
    metrics, AUROC (exact estimate, binned replicates) and paired contrasts vs E0."""
    dkey, stratum = job
    t0 = time.time()
    ds = G["DS"][dkey]
    PB = G["PB"][dkey]
    comps = G["COMPS"][dkey]
    names = [c["name"] for c in comps]
    NU = G["NU"]
    B, seed = G["B"], G["seed"]
    mask = strata_masks(ds)[stratum]
    idx = np.flatnonzero(mask)
    if len(idx) == 0:
        return {"dkey": dkey, "stratum": stratum, "res": None}
    srcs, sl = np.unique(ds["src"][idx], return_inverse=True)
    Gs = len(srcs)
    schemes = CP.make_schemes(srcs, B, seed, {"group": G["group_code"], "src": None,
                                              "genus": G["genus_code"]}, order=SCHEMES)
    snames = list(schemes)

    def ps(v):
        return np.bincount(sl, weights=np.asarray(v, float)[idx], minlength=Gs)
    X = CP.Cols(Gs)
    X.add("NCELL", np.bincount(sl, minlength=Gs))
    X.add("T", ps(ds["T"]))
    X.add("N_s", ps(ds["Ns"]))
    X.add("POS_s", ps(ds["POSs"]))
    X.add("F1U", ps(ds["F1U"]))
    st = PB["stats"]
    for c, nm in enumerate(names):
        blk = st[idx, c, :]
        X.add(f"{nm}[]", np.vstack([np.bincount(sl, weights=blk[:, j], minlength=Gs)
                                    for j in range(NSB)]).T)
    XM = X.matrix()
    S = XM.sum(0)
    SB = collections.OrderedDict(
        (s_, (W_ @ CP.agg_rows(XM, inv_, K_)) if W_ is not None else None)
        for s_, (inv_, W_, K_) in schemes.items())
    pt = derive_comp(S, X.L, names, NU)
    bts = [derive_comp(SB[s_], X.L, names, NU) if SB[s_] is not None else None
           for s_ in snames]
    # ---- AUROC on the union support
    c2l = np.full(ds["n_cells"], -1, np.int64)
    c2l[idx] = sl
    pts = np.flatnonzero(mask[PB["ucell"]])
    yy = PB["uy"][pts]
    lsrc = c2l[PB["ucell"][pts]]
    for c, nm in enumerate(names):
        q = PB["uq"][c, pts]
        pt[nm]["union_auroc"] = auroc_exact(q, yy)
        bins = np.where(q >= 1.0, NBA - 1, np.minimum((q * (NBA - 1)).astype(np.int64),
                                                       NBA - 2))
        Hp_all = np.bincount(bins[yy], minlength=NBA).astype(float)
        Hn_all = np.bincount(bins[~yy], minlength=NBA).astype(float)
        pt[nm]["union_auroc_binned"] = float(auroc_from_hist(Hp_all, Hn_all))
        for j, s_ in enumerate(snames):
            inv_, W_, K_ = schemes[s_]
            if W_ is None:
                continue
            cl = lsrc if inv_ is None else inv_[lsrc]
            key = cl * NBA + bins
            Hp = np.bincount(key[yy], minlength=K_ * NBA).reshape(K_, NBA).astype(float)
            Hn = np.bincount(key[~yy], minlength=K_ * NBA).reshape(K_, NBA).astype(float)
            a = auroc_from_hist(W_ @ Hp, W_ @ Hn)
            bts[j][nm]["union_auroc"] = a
            bts[j][nm]["union_auroc_binned"] = a
    # ---- paired contrasts vs E0 and between estimators (same draws)
    pairs = [(a, b) for a, b in PAIRS if a in names and b in names]
    pairs += [(f"E4({b})", b) for b in names
              if f"E4({b})" in names and b != "E0" and (f"E4({b})", b) not in pairs]
    for tree in [pt] + [b for b in bts if b is not None]:
        tree["contrasts_vs_E0"] = {nm: {m: tree[nm][m] - tree["E0"][m]
                                        for m in CONTRAST_METRICS}
                                   for nm in names if nm != "E0"}
        tree["contrasts_pairwise"] = {f"{a} - {b}": {m: tree[a][m] - tree[b][m]
                                                     for m in CONTRAST_METRICS}
                                      for a, b in pairs}
    res = CP.pack_schemes(pt, bts, snames)
    res["n_cells"] = int(len(idx))
    res["n_src"] = int(Gs)
    res["n_union_points"] = int(len(pts))
    res["n_universe_points"] = int(len(idx) * NU)
    res["ci_unit"] = snames[0]
    res["clusters"] = {s_: int(K_) for s_, (_i, _W, K_) in schemes.items()}
    res["seconds"] = round(time.time() - t0, 1)
    return {"dkey": dkey, "stratum": stratum, "res": res,
            "peak_rss_gb": round(rss_gb(), 2)}


def derive_grid(S, L, rows, NU):
    g = lambda k: S[..., L[k]]  # noqa: E731
    NC, T, Ns, POS, F1U = g("NCELL"), g("T"), g("N_s"), g("POS_s"), g("F1U")
    Nv = NC * NU
    o = {}
    with np.errstate(all="ignore"):
        e0 = g("v0[]")
        for v in rows:
            b = g(f"v{v}[]")
            d = {"universe_brier": b[..., 0] / Nv,
                 "universe_logloss_eps0.01": b[..., 1] / Nv,
                 "union_brier": b[..., 2] / Ns,
                 "union_slope_wls": (b[..., 5] - b[..., 3] * POS / Ns) /
                                    (b[..., 4] - b[..., 3] ** 2 / Ns),
                 "macro_F1_q_ge_0.5": b[..., 8] / NC,
                 "micro_F1_q_ge_0.5": 2 * b[..., 6] / (b[..., 7] + T),
                 "d_macro_F1_vs_union": b[..., 8] / NC - F1U / NC,
                 "d_micro_F1_vs_union": 2 * b[..., 6] / (b[..., 7] + T) - 2 * POS / (Ns + T),
                 "support_mean": b[..., 9] / NC,
                 "d_universe_brier_vs_E0": (b[..., 0] - e0[..., 0]) / Nv,
                 "d_union_brier_vs_E0": (b[..., 2] - e0[..., 2]) / Ns}
            o[G["V"][v]["name"]] = d
    return o


def grid_job(job):
    """Pass-A grid (every variant, pooled over folds, NOT nested): point estimates for
    every stratum; bootstrap CIs (3 schemes) for 'all'."""
    dkey, stratum, want_ci = job
    ds = G["DS"][dkey]
    st = G["PA"][dkey]["stats"]
    NU = G["NU"]
    mask = strata_masks(ds)[stratum]
    idx = np.flatnonzero(mask)
    srcs, sl = np.unique(ds["src"][idx], return_inverse=True)
    Gs = len(srcs)
    rows = list(range(len(G["V"])))

    def ps(v):
        return np.bincount(sl, weights=np.asarray(v, float)[idx], minlength=Gs)
    X = CP.Cols(Gs)
    X.add("NCELL", np.bincount(sl, minlength=Gs))
    X.add("T", ps(ds["T"]))
    X.add("N_s", ps(ds["Ns"]))
    X.add("POS_s", ps(ds["POSs"]))
    X.add("F1U", ps(ds["F1U"]))
    for v in rows:
        blk = st[idx, v, :]
        X.add(f"v{v}[]", np.vstack([np.bincount(sl, weights=blk[:, j], minlength=Gs)
                                    for j in range(NSA)]).T)
    XM = X.matrix()
    S = XM.sum(0)
    pt = derive_grid(S, X.L, rows, NU)
    if want_ci and G["B"]:
        schemes = CP.make_schemes(srcs, G["B"], G["seed"],
                                  {"group": G["group_code"], "src": None,
                                   "genus": G["genus_code"]}, order=SCHEMES)
        snames = list(schemes)
        bts = [derive_grid(W_ @ CP.agg_rows(XM, inv_, K_), X.L, rows, NU)
               if W_ is not None else None for s_, (inv_, W_, K_) in schemes.items()]
        res = CP.pack_schemes(pt, bts, snames)
    else:
        res = CP.pack_schemes(pt, [], [])
    return {"dkey": dkey, "stratum": stratum, "res": res, "n_cells": int(len(idx))}


# ============================================================================ outputs
def tojson(o):
    return CP.tojson(o)


def est(d, *ks):
    for k in ks:
        if d is None:
            return None
        d = d.get(k)
    return d


def fmt(d, k=4, sens=False):
    if not d or d.get("est") is None:
        return "n/a"
    s = f"{d['est']:.{k}f}"
    if d.get("lo") is not None:
        s += f" [{d['lo']:.{k}f}, {d['hi']:.{k}f}]"
    if sens:
        for nm in ("ci_genus", "ci_src"):
            v = d.get(nm)
            if v and v[0] is not None:
                s += f" ({nm[3:]} [{v[0]:.{k}f}, {v[1]:.{k}f}])"
    return s


def below0(d):
    return bool(d and d.get("hi") is not None and d["hi"] < 0)


def build_stage_o(ds, V, row, iso, fam_label, pop, cell_of_key):
    """stage_o_weights.json for one configuration (row of V; iso = all-fold map or None)."""
    v = V[row]
    fam = v["family"]
    gids = G["gids"]
    amps = collections.OrderedDict()
    for job in pop:
        i = cell_of_key[job["key"]]
        r = ds["reps_flat"][ds["reps_ptr"][i]:ds["reps_ptr"][i + 1]]
        reps = [gids[x] for x in r]
        if fam in ("E1", "E3"):
            g, ids = csr_slice(ds["lists"]["g100" if v["dedup"] == "genome" else "s100"],
                               ds["rec"][i])
            if len(g) == 0:              # no survivor in the first 100 hits: compute_Q's fallback
                g, ids = r, np.ones(len(r))
            g, ids = g[:v["K"]], ids[:v["K"]]
            w = np.exp(-(ids[0] - ids) / v["tau"])
            nb = [[gids[x], float(f"{wx:.12g}")] for x, wx in zip(g, w)]
        else:
            nb = [[x, 1.0] for x in reps]
        clade, alpha = None, 0.0
        if fam in ("E2", "E3"):
            ci = CLADES.index(v["clade"])
            j = ds["crow"][ci, i]
            if j >= 0:
                mem = ds["cm_flat"][ds["cm_ptr"][j]:ds["cm_ptr"][j + 1]]
                clade = {"rank": ds["cmeta"][j][1], "members": [gids[x] for x in mem]}
                alpha = float(v["alpha"])
        amps[job["key"]] = {"md5": job["md5"], "fold": int(job["fold"]), "neighbours": nb,
                            "clade": clade, "alpha": alpha, "n_reps": int(len(reps)),
                            "reps": reps}
    rule = None
    if fam in ("E2", "E3"):
        start = BACKOFF.index(v["clade"])
        rule = (" -> ".join(BACKOFF[start:]) + f"; start at {v['clade']}, back off one rank "
                f"while the anchor's clade (anchor = reps[0]) has < {MIN_CLADE} genomes in "
                f"the fold reference (_cache_ec_subset keys) or its panel label is a "
                f"placeholder ('{PLACEHOLDER}*'); none admissible -> clade null, alpha 0")
    meta = {"estimator": fam_label,
            "K": v["K"] if fam in ("E1", "E3") else None,
            "tau": v["tau"] if fam in ("E1", "E3") else None,
            "dedup": v["dedup"] if fam in ("E1", "E3") else None,
            "clade_rank_rule": rule,
            "alpha": float(v["alpha"]) if fam in ("E2", "E3") else None,
            "isotonic": ({"x": iso[0].tolist(), "y": iso[1].tolist()} if iso is not None
                         else None),
            "script_sha256": G["script_sha256"],
            "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z")}
    return {"meta": meta, "amplicons": amps}


def check_stage_o(obj, ds, row, cell_of_key, pop):
    """Recompute q from the file the way a consumer would and compare with compute_Q."""
    V = G["V"]
    gidx, E = G["gidx"], G["E"]
    bad = 0
    nchk = 0
    iso = obj["meta"]["isotonic"]
    for job in pop[:60]:
        a = obj["amplicons"][job["key"]]
        i = cell_of_key[job["key"]]
        Q, *_ = compute_Q(ds, i)
        num = np.zeros(G["NU"])
        den = 0.0
        for g, w in a["neighbours"]:
            num += w * E[gidx[g]]
            den += w
        if a["clade"] is not None:
            mem = [gidx[g] for g in a["clade"]["members"]]
            num += a["alpha"] * E[mem].mean(0)
            den += a["alpha"]
        q = num / den
        q[q >= 1 - 1e-9] = 1.0
        q[q <= 1e-12] = 0.0
        ref = Q[row]
        if iso is not None:
            q = apply_map(q, (np.asarray(iso["x"]), np.asarray(iso["y"])))
            ref = apply_map(ref, (np.asarray(iso["x"]), np.asarray(iso["y"])))
        nchk += 1
        bad += int(np.abs(q - ref).max() > 1e-6)
    return {"checked": nchk, "mismatch_gt_1e-6": bad}


# ============================================================================ main
def run_pool(fn, jobs, workers, label):
    out = []
    ctx = mp.get_context("fork")
    t0 = time.time()
    with ctx.Pool(workers, maxtasksperchild=8) as pool:
        for r in pool.imap_unordered(fn, jobs):
            out.append(r)
    log(f"{label}: {len(jobs)} jobs on {workers} workers ({time.time()-t0:.1f}s)")
    return out


def chunks_for(ds, size):
    jobs = []
    for f in ds["folds"]:
        cid = np.flatnonzero(ds["fold"] == f)
        lo, hi = int(cid[0]), int(cid[-1]) + 1
        assert hi - lo == len(cid)
        for a in range(lo, hi, size):
            jobs.append((ds["key"], a, min(a + size, hi)))
    return jobs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--densities", default="1.0,0.5,0.25")
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--limit-cells", type=int, default=0, help="dev: cells per fold")
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1729)
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--chunk", type=int, default=500)
    ap.add_argument("--outdir", default=os.path.join(D, "pfunc_estimators"))
    ap.add_argument("--cache-dir", default=None, help="scratch dir for the stream cache")
    ap.add_argument("--reuse-stream", action="store_true")
    ap.add_argument("--force-winner", default=None,
                    help="DEV ONLY: override the winner (tests the Stage O code paths)")
    args = ap.parse_args()
    T0 = time.time()
    densities = [d for d in args.densities.split(",") if d]
    folds = [int(f) for f in args.folds.split(",") if f]
    full = (folds == [0, 1, 2, 3, 4] and not args.limit_cells)
    if "1.0" not in densities:
        raise SystemExit("density 1.0 is the primary population")
    if len(folds) < 2:
        raise SystemExit("nested selection needs at least two folds")
    os.makedirs(args.outdir, exist_ok=True)
    G["B"], G["seed"] = args.bootstrap, args.seed
    G["script_sha256"] = sha256_file(os.path.abspath(__file__))
    V, VIX = build_variants()
    G["V"], G["VIX"] = V, VIX
    G["E1_rows"] = {di: [VIX["E1"][(d, t, K)] for t in TAU1 for K in K1]
                    for di, d in enumerate(DEDUPS)}
    timings = collections.OrderedDict()
    log(f"{len(V)} variants; densities {densities}; folds {folds}; limit "
        f"{args.limit_cells}; B {args.bootstrap}; workers {args.workers}; full={full}")

    load_context()
    cells = load_cells(densities, folds, args.limit_cells)
    timings["context_and_cells"] = round(time.time() - T0, 1)

    # ------------------------------------------------------------------ stream (once)
    cache = (os.path.join(args.cache_dir, f"stream_{'_'.join(DTAG[d] for d in densities)}"
                                          f"_f{''.join(map(str, folds))}_l{args.limit_cells}.pkl")
             if args.cache_dir else None)
    t1 = time.time()
    if args.reuse_stream and cache and os.path.exists(cache):
        with open(cache, "rb") as fh:
            parts = pickle.load(fh)
        log(f"stream cache reused: {cache}")
    else:
        jobs = [(DTAG[d], f, sorted(set(cells[d][cells[d].fold_i == f].md5)))
                for d in densities for f in folds]
        jobs.sort(key=lambda j: -len(j[2]))
        parts = run_pool(stream_job, jobs, min(args.workers, len(jobs)), "stream")
        if cache:
            os.makedirs(args.cache_dir, exist_ok=True)
            with open(cache, "wb") as fh:
                pickle.dump(parts, fh, protocol=5)
    stream_stats = [{"job": f"f{p['fold']:02d}_{p['dtag']}", "seconds": p["seconds"],
                     "peak_rss_gb": p["peak_rss_gb"]} for p in parts]
    timings["stream"] = round(time.time() - t1, 1)

    # ------------------------------------------------------------------ datasets
    t1 = time.time()
    G["DS"] = {}
    for d in densities:
        dk = DTAG[d]
        G["DS"][dk] = build_dataset(dk, cells[d], [p for p in parts if p["dtag"] == dk])
    del parts
    timings["datasets"] = round(time.time() - t1, 1)

    # ------------------------------------------------------------------ Stage O sample
    pop = jload(POPULATION)["jobs"]
    ds100 = G["DS"]["d100"]
    cell_of_key = {a: i for i, a in enumerate(ds100["cells"].amplicon)}
    popA = collections.Counter()
    for job in pop:
        i = cell_of_key.get(job["key"])
        if i is None:
            popA["key_not_in_scored_cells"] += 1
            continue
        r = [G["gids"][x] for x in ds100["reps_flat"][ds100["reps_ptr"][i]:
                                                      ds100["reps_ptr"][i + 1]]]
        popA["md5_mismatch"] += int(ds100["cells"].md5[i] != job["md5"])
        popA["fold_mismatch"] += int(int(ds100["fold"][i]) != int(job["fold"]))
        popA["reps_mismatch"] += int(r != list(job["reps"]))
        popA["src_mismatch"] += int(ds100["cells"].src_genome_id[i] != job["src_genome_id"])
        popA["group_mismatch"] += int(ds100["cells"].group_id[i] != job["group_id"])
    popA["jobs"] = len(pop)

    # ------------------------------------------------------------------ pass A
    t1 = time.time()
    G["PA"] = {}
    for dk, ds in G["DS"].items():
        res = run_pool(passA_job, chunks_for(ds, args.chunk), args.workers, f"pass A {dk}")
        st = np.zeros((ds["n_cells"], len(V), NSA))
        H = {f: np.zeros((3, len(V), NBF + 2)) for f in ds["folds"]}
        fb = 0
        for r in res:
            st[r["lo"]:r["hi"]] = r["stats"]
            H[r["fold"]] += r["H"]
            fb += r["fallback_cells"]
        G["PA"][dk] = {"stats": st, "H": H, "fallback_cells": fb,
                       "peak_rss_gb": max(r["peak_rss_gb"] for r in res)}
        del res
    timings["passA"] = round(time.time() - t1, 1)

    # ------------------------------------------------------------------ reproduction
    NU = G["NU"]
    PA = G["PA"]["d100"]
    stA = PA["stats"]
    NC100 = ds100["n_cells"]
    repro = collections.OrderedDict()
    e0 = stA[:, 0, :]
    repro["E0_universe_brier"] = float(e0[:, SA["SQu"]].sum() / (NC100 * NU))
    repro["E0_universe_logloss_eps0.01"] = float(e0[:, SA["LLu"]].sum() / (NC100 * NU))
    repro["E0_union_macro_F1"] = float(ds100["F1U"].mean())
    repro["E0_union_micro_F1"] = float(2 * ds100["POSs"].sum() /
                                       (ds100["Ns"].sum() + ds100["T"].sum()))
    repro["E0_union_support_brier"] = float(e0[:, SA["SQs"]].sum() / ds100["Ns"].sum())
    repro["E0_p_ge_0.5_macro_F1"] = float(e0[:, SA["F1"]].mean())
    for K in REP_K:
        b = stA[:, VIX["R200"][K], :]
        repro[f"E1_200hit_K{K}_tau0.01_genome_universe_brier"] = float(
            b[:, SA["SQu"]].sum() / (NC100 * NU))
        repro[f"E1_200hit_K{K}_tau0.01_genome_q_ge_0.5_macro_F1"] = float(b[:, SA["F1"]].mean())
    row20 = VIX["E1"][("genome", 0.01, 20)]
    repro["E1_100hit_K20_tau0.01_genome_universe_brier"] = float(
        stA[:, row20, SA["SQu"]].sum() / (NC100 * NU))
    repro["E1_100hit_K20_tau0.01_genome_q_ge_0.5_macro_F1"] = float(
        stA[:, row20, SA["F1"]].mean())
    calib = None
    targets = {}
    if os.path.exists(CALIB_SCORES):
        cs = jload(CALIB_SCORES)
        t = cs["exploratory_topk"]["all|all_records"]
        et = t["exact_tie_k_over_n"]
        targets = {"E0_universe_brier": et["brier_universe"]["est"],
                   "E0_universe_logloss_eps0.01": et["logloss_universe_eps0.01"]["est"],
                   "E0_union_macro_F1": et["union.macro_F1"]["est"],
                   "E0_union_micro_F1": et["union.micro_F1"]["est"],
                   "E0_union_support_brier": et["brier_union_support"]["est"],
                   "E0_p_ge_0.5_macro_F1": et["p_ge_0.5.macro_F1"]["est"]}
        for K in REP_K:
            targets[f"E1_200hit_K{K}_tau0.01_genome_universe_brier"] = \
                t[f"k={K}"]["brier_universe"]["est"]
            targets[f"E1_200hit_K{K}_tau0.01_genome_q_ge_0.5_macro_F1"] = \
                t[f"k={K}"]["q_ge_0.5.macro_F1"]["est"]
        pc = cs["primary_controlled"]
        for d in densities:
            if d != "1.0":
                targets[f"E0_universe_brier_{DTAG[d]}"] = \
                    pc[f"exact_tie|{DTAG[d]}"]["ec"]["all|all=all"]["universe"]["brier_p"]["est"]
        calib = {"path": CALIB_SCORES, "n_cells_exploratory": t["n_cells"]}
        del cs
    for d in densities:
        if d != "1.0":
            dk = DTAG[d]
            s_ = G["PA"][dk]["stats"]
            repro[f"E0_universe_brier_{dk}"] = float(s_[:, 0, SA["SQu"]].sum() /
                                                     (G["DS"][dk]["n_cells"] * NU))
    log("reproduction: " + ", ".join(f"{k}={v:.8g} (target {targets.get(k)})"
                                     for k, v in repro.items()))

    # ------------------------------------------------------------------ nested selection
    sel, LLfold, NCf = nested_select(PA, ds100)
    for fam in FAMILIES:
        log(f"selection {fam}: " + "; ".join(f"f{f}: {sel[fam]['per_fold'][f]['variant']}"
                                             for f in ds100["folds"]) +
            f" | carried forward {sel[fam]['carried_forward']['name']}")

    # ------------------------------------------------------------------ pass B d100
    t1 = time.time()
    G["COMPS"], G["PB"], MAPS = {}, {}, {}
    comps, maps_log = make_composites("d100", ds100, PA, sel)
    G["COMPS"]["d100"] = comps
    MAPS["d100"] = maps_log

    def run_passB(dk):
        ds = G["DS"][dk]
        res = run_pool(passB_job, chunks_for(ds, args.chunk), args.workers, f"pass B {dk}")
        res.sort(key=lambda r: r["lo"])
        st = np.zeros((ds["n_cells"], len(G["COMPS"][dk]), NSB))
        for r in res:
            st[r["lo"]:r["hi"]] = r["stats"]
        G["PB"][dk] = {"stats": st, "uq": np.concatenate([r["uq"] for r in res], 1),
                       "uy": np.concatenate([r["uy"] for r in res]),
                       "ucell": np.concatenate([r["ucell"] for r in res]),
                       "peak_rss_gb": max(r["peak_rss_gb"] for r in res)}
        del res

    run_passB("d100")
    # consistency: pass-B E0/composite stats equal pass-A rows
    pbchk = collections.Counter()
    stB = G["PB"]["d100"]["stats"]
    for c, comp in enumerate(comps):
        if comp["maps"] is not None:
            continue
        rows_cell = np.array([comp["rows"][int(f)] for f in ds100["fold"]])
        a = stA[np.arange(NC100), rows_cell][:, [SA["SQu"], SA["LLu"], SA["SQs"], SA["F1"]]]
        b = stB[:, c][:, [SB_["SQu"], SB_["LLu"], SB_["SQs"], SB_["F1"]]]
        pbchk[comp["name"]] = int((np.abs(a - b) > 1e-6 * np.maximum(1, np.abs(a))).sum())
    timings["passB_d100"] = round(time.time() - t1, 1)

    # ------------------------------------------------------------------ analysis d100
    t1 = time.time()
    ana = run_pool(analyze_job, [("d100", s) for s in STRATA],
                   min(args.workers, len(STRATA)), "analysis d100")
    comp_res = {"d100": {a["stratum"]: a["res"] for a in ana}}
    ana_peak = max(a.get("peak_rss_gb", 0) for a in ana)
    gjobs = [("d100", s, s == "all") for s in STRATA]
    for d in densities:
        if d != "1.0":
            gjobs.append((DTAG[d], "all", False))
    gr = run_pool(grid_job, gjobs, min(args.workers, len(gjobs)), "grid")
    grid = collections.defaultdict(dict)
    for g_ in gr:
        grid[g_["dkey"]][g_["stratum"]] = g_["res"]
    timings["analysis_d100"] = round(time.time() - t1, 1)

    # ------------------------------------------------------------------ decision, winner
    allr = comp_res["d100"]["all"]
    decision = collections.OrderedDict()
    cands = []
    for c in comps:
        nm = c["name"]
        if nm == "E0":
            continue
        cv = allr["contrasts_vs_E0"][nm]
        ub, sb = cv["universe_brier"], cv["union_brier"]
        beats = below0(ub) and below0(sb)
        decision[nm] = {"d_universe_brier": ub, "d_union_brier": sb, "beats_k_over_n": beats,
                        "universe_logloss": allr[nm]["universe_logloss_eps0.01"]["est"],
                        "universe_brier": allr[nm]["universe_brier"]["est"],
                        "union_brier": allr[nm]["union_brier"]["est"]}
        if beats:
            cands.append(nm)
    rank = {k: sorted(cands, key=lambda n: decision[n][k]) for k in
            ("universe_logloss", "universe_brier", "union_brier")}
    winner = rank["universe_logloss"][0] if cands else "E0"
    rule_dependent = bool(cands) and len({r[0] for r in rank.values()}) > 1
    if args.force_winner:
        if full:
            raise SystemExit("--force-winner is a development flag; not allowed on the full run")
        log(f"DEV: winner forced {winner} -> {args.force_winner}")
        winner = args.force_winner
    win_comp = next(c for c in comps if c["name"] == winner)
    base = win_comp["base"]
    base_fam = base if base in FAMILIES else "E0"
    if base_fam == "E0":
        cf_row = 0
    else:
        cf_row = sel[base_fam]["carried_forward_row"]
    iso_all = None
    if winner.startswith("E4"):
        Hs = sum(PA["H"][f][:, cf_row, :] for f in ds100["folds"])
        iso_all = fit_iso(Hs)
    log(f"decision: beats k/n: {cands}; winner {winner} (base {base_fam}, carried-forward "
        f"{V[cf_row]['name']}); ranking {rank}; rule-dependent {rule_dependent}")

    # ------------------------------------------------------------------ densities
    t1 = time.time()
    dens_res = {}
    for d in densities:
        if d == "1.0":
            continue
        dk = DTAG[d]
        ds = G["DS"][dk]
        fixed = collections.OrderedDict()
        for fam in FAMILIES:
            fixed[f"{fam}*"] = sel[fam]["carried_forward_row"]
        fixed["_e4"] = ([f"{base_fam}*"] if (winner.startswith("E4") and base_fam != "E0")
                        else [])
        cmp_, mlog = make_composites(dk, ds, G["PA"][dk], sel, fixed=fixed)
        G["COMPS"][dk] = cmp_
        MAPS[dk] = mlog
        run_passB(dk)
        ana = run_pool(analyze_job, [(dk, s) for s in STRATA],
                       min(args.workers, len(STRATA)), f"analysis {dk}")
        dens_res[dk] = {a["stratum"]: a["res"] for a in ana}
        ana_peak = max(ana_peak, max(a.get("peak_rss_gb", 0) for a in ana))
    comp_res.update(dens_res)
    timings["densities"] = round(time.time() - t1, 1)

    # ------------------------------------------------------------------ Stage O weights
    stage_o = {}
    pop_here = [j for j in pop if j["key"] in cell_of_key]       # all 360 on the full run
    if pop_here:
        lab = "E4" if winner.startswith("E4") else base_fam
        stage_o["stage_o_weights.json"] = build_stage_o(ds100, V, cf_row, iso_all, lab,
                                                        pop_here, cell_of_key)
        if winner != "E1":
            stage_o["stage_o_weights_E1.json"] = build_stage_o(
                ds100, V, sel["E1"]["carried_forward_row"], None, "E1", pop_here, cell_of_key)
    so_chk = {}
    for fn, obj in stage_o.items():
        row = cf_row if fn == "stage_o_weights.json" else sel["E1"]["carried_forward_row"]
        so_chk[fn] = check_stage_o(obj, ds100, row, cell_of_key, pop_here)
        so_chk[fn]["amplicons"] = len(obj["amplicons"])
    log(f"stage O weights: {[(fn, len(o['amplicons'])) for fn, o in stage_o.items()]}; "
        f"consumer check {so_chk}")

    # ------------------------------------------------------------------ assertions
    req = collections.OrderedDict()
    tol = {"E0_universe_brier": 1e-6, "E0_union_macro_F1": 1e-6}
    if full:
        for k in ("E0_universe_brier", "E0_union_macro_F1"):
            req[f"a_{k}_reproduces_calibration_within_{tol[k]:g}"] = bool(
                k in targets and abs(repro[k] - targets[k]) <= tol[k])
        k20 = "E1_200hit_K20_tau0.01_genome_universe_brier"
        req["b_E1_K20_tau0.01_genome_200hit_universe_brier_reproduces_0.06273832"] = bool(
            k20 in targets and abs(repro[k20] - targets[k20]) <= 1e-8)
        for k in ("E0_universe_logloss_eps0.01", "E0_union_micro_F1", "E0_union_support_brier",
                  "E0_p_ge_0.5_macro_F1"):
            req[f"c_{k}_reproduces_calibration_within_1e-6"] = bool(
                k in targets and abs(repro[k] - targets[k]) <= 1e-6)
        for d in densities:
            if d != "1.0":
                k = f"E0_universe_brier_{DTAG[d]}"
                req[f"d_{k}_reproduces_calibration_within_1e-6"] = bool(
                    k in targets and abs(repro[k] - targets[k]) <= 1e-6)
    diag_repro = {k: {"value": v, "target": targets.get(k),
                      "abs_diff": (abs(v - targets[k]) if k in targets else None)}
                  for k, v in repro.items()}
    for dk, ds in G["DS"].items():
        a = ds["asserts"]

        def s_(prefix):
            return sum(v for k, v in a.items() if k.startswith(prefix) and
                       isinstance(v, (int, np.integer)))
        req[f"e_{dk}_source_never_a_neighbour"] = s_("leak_source_in_neighbours")
        for nm in ("genus", "species", "taxon", "group"):
            req[f"f_{dk}_no_neighbour_shares_{nm}_with_source"] = s_(f"leak_neighbour_shares_{nm}")
        req[f"g_{dk}_every_neighbour_in_fold_reference"] = s_("neighbours_outside_fold_reference")
        req[f"h_{dk}_source_never_a_clade_member"] = s_("leak_source_in_clade")
        for nm in ("genus", "species", "taxon", "group"):
            req[f"i_{dk}_no_clade_member_shares_{nm}_with_source"] = s_(
                f"leak_clade_member_shares_{nm}")
        req[f"j_{dk}_every_clade_member_in_fold_reference"] = s_("clade_members_outside_fold_reference")
        req[f"j2_{dk}_clade_drawn_from_the_cells_fold"] = a.get("clade_fold_ne_cell_fold", 0)
        req[f"k_{dk}_anchor_genus_differs_from_source_genus"] = a.get(
            "leak_anchor_genus_eq_source_genus", 0)
        req[f"l_{dk}_reps_in_fold_reference"] = a.get("reps_outside_fold_reference", 0)
        req[f"m_{dk}_len_reps_equals_cells_and_pred_n_reps"] = a.get("n_reps_mismatch", 0)
        req[f"n_{dk}_union_of_reps_EC_equals_pred_ec"] = a.get("union_ec_ne_pred_ec", 0)
        req[f"o_{dk}_recomputed_union_f1_equals_cells_f1_within_1e-3"] = a.get(
            "f1_mismatch_vs_cells", 0)
        req[f"p_{dk}_recomputed_n_truth_equals_cells"] = a.get("n_truth_mismatch_vs_cells", 0)
        req[f"q_{dk}_cells_group_id_equals_fold_plan_genome_group"] = a.get(
            "group_id_ne_fold_plan_genome_group", 0)
        req[f"r_{dk}_ec_subset_cache_equals_panel_truth"] = a.get(
            "stream_ec_subset_mismatch_vs_panel_truth", 0) + a.get(
            "stream_ec_subset_genome_outside_panel", 0)
        req[f"s_{dk}_hits_in_align_score_rank_order"] = a.get("stream_hits_rank_order_violation", 0)
        req[f"t_{dk}_hit_genomes_in_panel"] = sum(v for k, v in a.items()
                                                  if k.startswith("stream_hit_genome_outside"))
        req[f"u_{dk}_200hit_tie_set_equals_exact_tie_reps"] = ds["diag"].get("tie_set_200_ne_reps", 0)
    req["v_panel_genomes_without_fold_plan_group"] = G["A_ctx"]["panel_genomes_without_fold_plan_group"]
    req["w_genera_spanning_more_than_one_group"] = G["A_ctx"]["panel_taxonomy"]["genera_spanning_more_than_one_group"]
    req["x_population_jobs_match_scored_cells_and_reps"] = (
        sum(v for k, v in popA.items() if k != "jobs") if full else 0)
    if full:
        req["x2_stage_o_weights_cover_all_360_amplicons"] = int(
            any(len(o["amplicons"]) != len(pop) for o in stage_o.values()) or not stage_o)
    req["y_passB_rows_equal_passA_rows"] = sum(pbchk.values())
    for fn, c in so_chk.items():
        req[f"z_{fn}_consumer_recomputation_matches"] = c["mismatch_gt_1e-6"]
    all_ok = all((v is True) or (v == 0 and v is not False) for v in req.values())
    assertions = {"all_required_zero_or_true": bool(all_ok), "full_population": full,
                  "required": req, "reproduction": diag_repro,
                  "per_dataset": {dk: ds["asserts"] for dk, ds in G["DS"].items()},
                  "diagnostics": {dk: {**ds["diag"],
                                       "survivors_100_quantiles": ds["survivors_100_quantiles"],
                                       "passA_cells_without_100hit_survivors_fallback_k_over_n":
                                           G["PA"][dk]["fallback_cells"]}
                                  for dk, ds in G["DS"].items()},
                  "population_check": dict(popA), "passB_vs_passA": dict(pbchk),
                  "stage_o_consumer_check": so_chk, "context": G["A_ctx"],
                  "notes": {
                      "a-d": "reproductions of data/calibration/scores.json (exploratory_topk "
                             "all|all_records exact_tie_k_over_n and k=20; primary_controlled "
                             "exact_tie|d050/d025 universe brier_p); b uses the 200-hit records "
                             "as calibrate did, the 100-hit value is in 'reproduction'",
                      "e-k": "LEAKAGE, every scored cell: neighbours = every genome of the "
                             "100-hit genome- and species-dedup lists and the 200-hit "
                             "reproduction list; clade members = both clade starts (genus, "
                             "family) after back-off; codes from panel_eligible (genus, species, "
                             "taxon_id) and fold_plan genome_group; fold reference = "
                             "_cache_ec_subset.json keys of the cell's fold run",
                      "l-r": "joins, as calibrate_pfunc layer A",
                      "u": "stage 0 + identity ranking + genome dedup of the 200-hit record "
                           "reproduces the exact-tie reps as the identity-1 tie set",
                      "x": "population.json (Stage O sample) keys are scored d100 cells with "
                           "the same md5, fold, reps, source and group",
                      "y": "pass B recomputes the composites' rows: they equal pass A's",
                      "z": "q recomputed from stage_o_weights.json with the consumer formula "
                           "equals the estimator's q (60 amplicons)"}}

    # ------------------------------------------------------------------ write
    sel_out = {"rule": "per outer fold, each family's variant minimises the pooled universe "
                       "log-loss (eps 0.01) of the other four folds; exact ties: smaller K, "
                       "larger tau, genome before species, genus before family, smaller alpha. "
                       "Carried forward: the configuration chosen in most outer folds (same "
                       "tie order); per-hyperparameter modes are reported for reference.",
               "families": sel,
               "E4_isotonic_maps_out_of_fold": MAPS,
               "winner": {"name": winner, "base_family": base_fam,
                          "carried_forward_variant": V[cf_row]["name"],
                          "isotonic_all_folds": ({"x": iso_all[0].tolist(),
                                                  "y": iso_all[1].tolist()}
                                                 if iso_all is not None else None),
                          "beats_k_over_n": cands, "ranking": rank,
                          "rule_dependent": rule_dependent,
                          "rule": "among composites whose nested-CV universe-Brier and "
                                  "union-support-Brier contrasts vs E0 both have group CIs "
                                  "below 0 (the design's decision rule), the lowest pooled "
                                  "nested-CV universe log-loss (the tuning criterion)"},
               "train_universe_logloss_by_variant_and_outer_fold": {
                   V[v]["name"]: {str(f): float((LLfold[v].sum() - LLfold[v, fi]) /
                                                ((NCf.sum() - NCf[fi]) * NU))
                                  for fi, f in enumerate(ds100["folds"])}
                   for v in range(len(V)) if V[v]["family"] in FAMILIES}}
    meta = {"script": os.path.abspath(__file__), "script_sha256": G["script_sha256"],
            "design": DESIGN, "design_sha256": sha256_file(DESIGN),
            "git_commit": CP.git_commit(), "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
            "args": vars(args), "full_population": full,
            "grids": {"E1": {"K": K1, "tau": TAU1, "dedup": DEDUPS},
                      "E2": {"clade_start": CLADES, "alpha": ALPHA2,
                             "backoff": BACKOFF, "min_clade": MIN_CLADE},
                      "E3": {"K": K3, "tau": TAU3, "dedup": DEDUPS, "clade_start": CLADES,
                             "alpha": ALPHA3}},
            "hits_used": HITS_USED, "eps": EPS, "isotonic_bins": NBF, "auroc_bins": NBA,
            "bootstrap": {"replicates": args.bootstrap, "seed": args.seed,
                          "draw": "calibrate_pfunc.boot_weights (score_folds.cluster_bootstrap"
                                  " draw verbatim), percentile 95%",
                          "unit": "PRIMARY holdout group -> lo/hi; sensitivities source "
                                  "genome (ci_src) and genus (ci_genus)",
                          "sesoi_delta_f1": G["SESOI"]},
            "timings_seconds": timings, "stream_jobs": stream_stats}
    scores = {"meta": meta, "reproduction": diag_repro,
              "composites": comp_res, "decision_d100_all": decision,
              "winner": sel_out["winner"],
              "composite_definitions": {dk: [{"name": c["name"], "base": c["base"],
                                              "rows_by_fold": {str(f): V[r]["name"]
                                                               for f, r in c["rows"].items()},
                                              "isotonic": c["maps"] is not None}
                                             for c in cl] for dk, cl in G["COMPS"].items()},
              "grid_pooled_not_nested": grid}
    timings["total_before_write"] = round(time.time() - T0, 1)
    meta["peak_rss_gb"] = {"parent": round(rss_gb(), 2),
                           "largest_child": round(resource.getrusage(
                               resource.RUSAGE_CHILDREN).ru_maxrss / 1e6, 2)}
    with open(os.path.join(args.outdir, "scores.json"), "wb") as fh:
        fh.write(orjson.dumps(tojson(scores), option=orjson.OPT_INDENT_2))
    with open(os.path.join(args.outdir, "assertions.json"), "wb") as fh:
        fh.write(orjson.dumps(tojson(assertions), option=orjson.OPT_INDENT_2))
    with open(os.path.join(args.outdir, "selected.json"), "wb") as fh:
        fh.write(orjson.dumps(tojson(sel_out), option=orjson.OPT_INDENT_2))
    for fn, obj in stage_o.items():
        with open(os.path.join(args.outdir, fn), "wb") as fh:
            fh.write(orjson.dumps(tojson(obj), option=orjson.OPT_INDENT_2))
    write_readme(os.path.join(args.outdir, "README.md"), tojson(scores), tojson(assertions),
                 tojson(sel_out), {fn: tojson(o["meta"]) for fn, o in stage_o.items()})
    log(f"wrote {args.outdir}: scores.json, assertions.json, selected.json, README.md"
        + "".join(f", {fn}" for fn in stage_o) +
        f"; all required zero/true: {all_ok}; total {time.time()-T0:.0f}s; peak rss "
        f"parent {meta['peak_rss_gb']['parent']} GB, largest child "
        f"{meta['peak_rss_gb']['largest_child']} GB")
    for k, v in req.items():
        log(f"  assert {k}: {v}")
    if not all_ok:
        raise SystemExit("[pfunc] a required assertion failed (see assertions.json)")


# ============================================================================ README
def write_readme(path, scores, assertions, sel, so_meta):
    L = []
    meta = scores["meta"]
    comp = scores["composites"]
    d100 = comp.get("d100", {})
    names = [c["name"] for c in scores["composite_definitions"]["d100"]]
    W = L.append
    W("# Stage E: function-probability estimators vs the shipped k/n")
    W("")
    W(f"`scripts/pfunc_estimators.py` (sha256 `{meta['script_sha256'][:12]}`), "
      f"{meta['generated']}, git `{meta['git_commit'][:10]}`. Design: "
      f"`PFUNC_ESTIMATORS_AND_OBJECTIVE_TUNING.md` section \"Stage E\" (sha256 "
      f"`{meta['design_sha256'][:12]}`). **Post hoc** (declared after the calibration "
      f"study's exploratory top-k had been seen); every tuned choice is nested by outer fold.")
    W("")
    W(f"All required assertions zero/true: **{assertions['all_required_zero_or_true']}** "
      f"(full population: {assertions['full_population']}).")
    W("")
    # ---------------- key results
    W("## Key results (d100, nested CV, all 57,021 cells unless stated)")
    W("")
    W("95% CIs: cluster bootstrap over holdout groups (2,000 draws, seed 1729, score_folds' "
      "draw); genus / source-genome intervals in scores.json. Contrasts are paired (same draws).")
    W("")
    W("| estimator | universe Brier | d vs E0 | universe log-loss | d vs E0 | union Brier | d vs E0 "
      "| union WLS slope | union AUROC | macro dF1 vs union | beats k/n |")
    W("|---|---|---|---|---|---|---|---|---|---|---|")
    al = d100.get("all", {})
    dec = scores.get("decision_d100_all", {})
    for nm in names:
        r = al.get(nm, {})
        cv = al.get("contrasts_vs_E0", {}).get(nm, {})
        W(f"| {nm} | {fmt(r.get('universe_brier'))} | {fmt(cv.get('universe_brier'))} | "
          f"{fmt(r.get('universe_logloss_eps0.01'))} | {fmt(cv.get('universe_logloss_eps0.01'))} | "
          f"{fmt(r.get('union_brier'))} | {fmt(cv.get('union_brier'))} | "
          f"{fmt(r.get('union_slope_wls'), 3)} | {fmt(r.get('union_auroc'), 3)} | "
          f"{fmt(r.get('d_macro_F1_vs_union'))} | "
          f"{'yes' if dec.get(nm, {}).get('beats_k_over_n') else ('-' if nm == 'E0' else 'no')} |")
    W("")
    win = scores["winner"]
    W(f"**Winner: {win['name']}** (rule: {win['rule']}). Estimators that beat k/n by the "
      f"design's rule: {', '.join(win['beats_k_over_n']) or 'none'}. Ranking by universe "
      f"log-loss / universe Brier / union Brier: {win['ranking']}; the choice is "
      f"{'RULE-DEPENDENT' if win['rule_dependent'] else 'the same under all three orderings'}. "
      f"Carried-forward configuration: `{win['carried_forward_variant']}`"
      + (" with an all-fold isotonic map." if win.get("isotonic_all_folds") else "."))
    W("")
    wb = f"{win['name']} - {win['base_family']}"
    wv = al.get("contrasts_pairwise", {}).get(wb)
    if wv:
        W(f"Winner vs its own base ({wb}, paired): universe Brier {fmt(wv.get('universe_brier'))}; "
          f"universe log-loss {fmt(wv.get('universe_logloss_eps0.01'))}; union Brier "
          f"{fmt(wv.get('union_brier'))}; union WLS slope {fmt(wv.get('union_slope_wls'), 3)}; "
          f"union AUROC {fmt(wv.get('union_auroc'), 3)}; macro dF1 "
          f"{fmt(wv.get('d_macro_F1_vs_union'))}. A Brier CI that covers 0 here means the "
          f"recalibration step is not distinguishable from its base on that score; "
          f"`stage_o_weights.json` carries the base configuration in full, so dropping its "
          f"isotonic map gives the base estimator.")
        W("")
    # ---------------- sensitivity CIs for the headline contrasts
    W("Headline contrasts vs E0 with all three bootstrap units (group [primary]; genus; "
      "source genome):")
    W("")
    for nm in names:
        if nm == "E0":
            continue
        cv = al.get("contrasts_vs_E0", {}).get(nm, {})
        W(f"- {nm}: universe Brier {fmt(cv.get('universe_brier'), 4, True)}; union Brier "
          f"{fmt(cv.get('union_brier'), 4, True)}")
    W("")
    # ---------------- pairwise
    pw = al.get("contrasts_pairwise", {})
    if pw:
        W("### Paired contrasts between estimators (all cells; group CI)")
        W("")
        W("| contrast | universe Brier | universe log-loss | union Brier | union WLS slope | "
          "union AUROC | macro dF1 vs union |")
        W("|---|---|---|---|---|---|---|")
        for k_, v in pw.items():
            W(f"| {k_} | {fmt(v.get('universe_brier'))} | "
              f"{fmt(v.get('universe_logloss_eps0.01'))} | {fmt(v.get('union_brier'))} | "
              f"{fmt(v.get('union_slope_wls'), 3)} | {fmt(v.get('union_auroc'), 3)} | "
              f"{fmt(v.get('d_macro_F1_vs_union'))} |")
        W("")
        W("Same contrasts on n = 1 and n >= 2 (universe Brier; union Brier):")
        W("")
        W("| contrast | n1 universe | n1 union | multi universe | multi union |")
        W("|---|---|---|---|---|")
        for k_ in pw:
            r1 = est(d100.get("n1"), "contrasts_pairwise", k_) or {}
            rm = est(d100.get("multi"), "contrasts_pairwise", k_) or {}
            W(f"| {k_} | {fmt(r1.get('universe_brier'))} | {fmt(r1.get('union_brier'))} | "
              f"{fmt(rm.get('universe_brier'))} | {fmt(rm.get('union_brier'))} |")
        W("")
    # ---------------- by n
    W("### By selection size (n = number of exact-tie reps)")
    W("")
    W("The WLS slope is undefined for E0 at n = 1 (q = 1 for every union EC) and unstable for "
      "E4(E0) there (its q at n = 1 is the per-fold constant r(1)); AUROC of E0 at n = 1 is 0.5 "
      "by construction.")
    W("")
    for metric, lab, k in (("universe_brier", "universe Brier", 4),
                           ("union_brier", "union-support Brier", 4),
                           ("universe_logloss_eps0.01", "universe log-loss", 4),
                           ("universe_ece_decile", "universe ECE (deciles)", 4),
                           ("d_macro_F1_vs_union", "macro dF1 vs the union (q >= 0.5)", 4),
                           ("d_micro_F1_vs_union", "micro dF1 vs the union (q >= 0.5)", 4),
                           ("union_auroc", "union-support AUROC", 3),
                           ("union_slope_wls", "union-support WLS slope", 3)):
        W(f"**{lab}** (estimate [group CI]; contrast vs E0 below):")
        W("")
        W("| estimator | " + " | ".join(STRATA) + " |")
        W("|---|" + "---|" * len(STRATA))
        for nm in names:
            W(f"| {nm} | " + " | ".join(fmt(est(d100.get(s), nm, metric), k)
                                       for s in STRATA) + " |")
        for nm in names:
            if nm == "E0":
                continue
            W(f"| {nm} - E0 | " + " | ".join(
                fmt(est(d100.get(s), "contrasts_vs_E0", nm, metric), k) for s in STRATA) + " |")
        W("")
    W("Cells per stratum: " + ", ".join(f"{s} {d100.get(s, {}).get('n_cells', 0):,}"
                                        for s in STRATA) + ".")
    W("")
    W("Decision metrics vs the union, read against the SESOI of "
      f"{meta['bootstrap']['sesoi_delta_f1']} with score_folds' rule (macro dF1, all cells):")
    W("")
    for nm in names:
        if nm == "E0":
            continue
        d = est(d100.get("all"), nm, "d_macro_F1_vs_union")
        dm = est(d100.get("all"), nm, "d_micro_F1_vs_union")
        ss = meta['bootstrap']['sesoi_delta_f1']
        W(f"- {nm}: macro {fmt(d)} -> {CP.sesoi_verdict(d, ss)}; "
          f"micro {fmt(dm)} -> {CP.sesoi_verdict(dm, ss)}")
    W("")
    # ---------------- selection
    W("## Nested selection")
    W("")
    W(sel["rule"])
    W("")
    W("| family | " + " | ".join(f"fold {f}" for f in sorted(
        next(iter(sel["families"].values()))["per_fold"])) + " | carried forward |")
    W("|---|" + "---|" * (len(next(iter(sel["families"].values()))["per_fold"]) + 1))
    for fam, s in sel["families"].items():
        W(f"| {fam} | " + " | ".join(f"`{s['per_fold'][f]['variant']}`"
                                    for f in sorted(s["per_fold"])) +
          f" | `{s['carried_forward']['name']}` |")
    W("")
    W("Grid edges (carried-forward value relative to the declared grid; a value at the edge "
      "means the optimum may lie outside the grid, which this study does not explore):")
    W("")
    for fam, s in sel["families"].items():
        W(f"- {fam}: " + "; ".join(f"{hp} chosen {e['chosen']} of {e['grid']} -> carried "
                                   f"forward at the {e['carried_forward_at']}"
                                   for hp, e in s.get("grid_edges", {}).items()))
    W("")
    # ---------------- grid
    W("## Grid (pooled over folds, NOT nested; point estimates, all cells)")
    W("")
    g = scores.get("grid_pooled_not_nested", {}).get("d100", {}).get("all", {})
    rows = [(k, v) for k, v in g.items() if isinstance(v, dict) and "universe_brier" in v]
    rows.sort(key=lambda kv: kv[1]["universe_logloss_eps0.01"]["est"])
    W("| variant | universe Brier | universe log-loss | union Brier | union slope | macro dF1 vs union |")
    W("|---|---|---|---|---|---|")
    for k, v in rows[:25]:
        W(f"| `{k}` | {fmt(v['universe_brier'])} | {fmt(v['universe_logloss_eps0.01'])} | "
          f"{fmt(v['union_brier'])} | {fmt(v['union_slope_wls'], 3)} | "
          f"{fmt(v['d_macro_F1_vs_union'])} |")
    W("")
    W(f"(top 25 of {len(rows)} by universe log-loss; the full grid with all strata is in "
      f"scores.json `grid_pooled_not_nested`.)")
    W("")
    # ---------------- densities
    dks = [k for k in comp if k != "d100"]
    if dks:
        W("## Density sensitivity (carried-forward configurations, no re-tuning)")
        W("")
        W("| density | estimator | universe Brier d vs E0 | union Brier d vs E0 | "
          "universe log-loss d vs E0 | macro dF1 vs union |")
        W("|---|---|---|---|---|---|")
        for dk in dks:
            a = comp[dk].get("all", {})
            for nm in [c["name"] for c in scores["composite_definitions"][dk]]:
                if nm == "E0":
                    W(f"| {dk} | E0 | universe Brier {fmt(a.get('E0', {}).get('universe_brier'))} "
                      f"| union Brier {fmt(a.get('E0', {}).get('union_brier'))} | | |")
                    continue
                cv = a.get("contrasts_vs_E0", {}).get(nm, {})
                W(f"| {dk} | {nm} | {fmt(cv.get('universe_brier'))} | "
                  f"{fmt(cv.get('union_brier'))} | {fmt(cv.get('universe_logloss_eps0.01'))} | "
                  f"{fmt(a.get(nm, {}).get('d_macro_F1_vs_union'))} |")
            for k_, v in a.get("contrasts_pairwise", {}).items():
                W(f"| {dk} | {k_} | {fmt(v.get('universe_brier'))} | "
                  f"{fmt(v.get('union_brier'))} | {fmt(v.get('universe_logloss_eps0.01'))} | "
                  f"{fmt(v.get('d_macro_F1_vs_union'))} (difference) |")
        W("")
    # ---------------- reproduction + assertions
    W("## Reproduction and assertions")
    W("")
    W("| quantity | this script | calibration scores.json |")
    W("|---|---|---|")
    for k, v in scores["reproduction"].items():
        t = v.get("target")
        W(f"| {k} | {v['value']:.8g} | {('%.8g' % t) if t is not None else '-'} |")
    W("")
    bad = [k for k, v in assertions["required"].items() if not (v is True or v == 0)]
    W(f"Required assertions: {len(assertions['required'])}, failing: {bad or 'none'}. "
      "Leakage (every cell, every neighbour list and both clade starts): the source genome is "
      "never a neighbour or clade member; no neighbour or clade member shares its genus, "
      "species, taxon_id or holdout group; every one is in the fold reference; the anchor's "
      "genus differs from the source's.")
    W("")
    dg = assertions["diagnostics"].get("d100", {})
    W(f"Diagnostics (d100): {dg}")
    W("")
    pdl = assertions["per_dataset"].get("d100", {})
    W("Clade back-off (cells at each rank): genus start "
      f"{pdl.get('clade_start_genus_cells_at_rank')}; family start "
      f"{pdl.get('clade_start_family_cells_at_rank')}. Cells without an admissible clade up "
      f"to order (alpha term dropped) by reason: {pdl.get('clade_none_reason_at_order')}.")
    W("")
    # ---------------- stage O
    if so_meta:
        W("## Stage O weights")
        W("")
        for fn, m in so_meta.items():
            mm = {k: (v if k != "isotonic" else (None if v is None else
                                                 f"{len(v['x'])}-knot map"))
                  for k, v in m.items()}
            W(f"- `{fn}`: {mm}")
        W("")
    # ---------------- deviations
    W("## Deviations and interpretations")
    W("")
    for s in DEVIATIONS:
        W(f"- {s}")
    W("")
    W("## Caveats")
    W("")
    for s in CAVEATS:
        W(f"- {s}")
    W("")
    W(f"Runtime {meta['timings_seconds']}; peak RSS {meta.get('peak_rss_gb')}.")
    with open(path, "w") as fh:
        fh.write("\n".join(L) + "\n")


DEVIATIONS = [
    "Hits truncation: the first 100 hits of each record in align-score rank order (rank "
    "1..100, asserted), then stage 0, then identity/score ranking and dedup. K = 100 means "
    "every survivor of those 100 hits (fewer than 100 after stage 0 and dedup). The E1 "
    "reproduction assertion uses the 200-hit records exactly as calibrate_pfunc did; the "
    "100-hit value is reported beside it.",
    "Species dedup uses select_references' species key on the hits' BV-BRC lineage "
    "(_norm(Species) -> taxon -> genome). The benchmark panel is species-dereplicated "
    "(8,872 genomes = 8,872 species), so the species variants are expected to equal the "
    "genome variants (records where they differ are counted in the diagnostics); exact ties "
    "resolve to genome dedup.",
    "Anchor = reps[0] of predictions.json (prFBA's own anchor: highest identity-decayed "
    "representativeness, then genome id). The clade includes the anchor and any other rep in it.",
    "Clades come from panel_eligible lineage (not the hits' lineage). Labels starting "
    "'unclassified_' (panel placeholders, e.g. 'unclassified_family:Lactobacillaceae' as a "
    "genus) are treated as 'no clade at this rank' and back off, like a clade with < 3 "
    "reference genomes. The back-off chain stops at order; with no admissible clade the "
    "alpha term is dropped (E2 = k/n, E3 = E1 for that cell).",
    "E4 isotonic is fitted on binned universe points (exact 0, 4,000 uniform bins, exact 1; "
    "x = bin mean, weight = points), with sklearn PAVA, and applied to the outer fold's exact "
    "q by the fitted piecewise-linear map. The map for outer fold f uses only the other four "
    "folds and the variant selected for f.",
    "Stage O's single isotonic map (if an E4 estimator wins) is fitted on all five folds "
    "with the carried-forward configuration; it is therefore in-sample for the 360 benchmark "
    "amplicons (0.6% of cells); the nested-CV numbers above are the honest ones.",
    "Winner rule (the design names none): among estimators whose universe-Brier and "
    "union-Brier contrasts vs E0 both have group CIs below 0, the lowest nested-CV universe "
    "log-loss (the tuning criterion). Carried-forward tie-breaks beyond 'smaller K, larger "
    "tau': genome before species, genus before family, smaller alpha.",
    "AUROC point estimates are exact (ties 1/2); bootstrap replicates use 1,001-bin "
    "histograms (scores.json gives the binned point estimate too).",
    "Universe ECE uses calibrate's decile bins with exact 0 and exact 1 as their own bins "
    "(for unrecalibrated q this is calibrate's 'ece_decile_universe'). ECE is a sum of "
    "absolute bin errors and is biased upward under resampling, so for a well-calibrated "
    "(E4) estimator its percentile interval can lie above the point estimate; read it as a "
    "descriptive number, not a contrast.",
    "Cells whose anchor has no admissible clade up to order (placeholder order label or an "
    "order with < 3 fold-reference genomes; counts in the diagnostics) keep the alpha term "
    "out: E2 = k/n and E3 = E1 there. The schema's clade ranks stop at order, so the "
    "back-off was not extended to class.",
    "Density sensitivity: the carried-forward configurations at d050/d025 without re-tuning "
    "(E4 refitted out of fold within the density); the clade prior uses that density's fold "
    "reference.",
]

CAVEATS = [
    "Post hoc: declared after the calibration study's exploratory top-k analysis had been "
    "seen. Nested selection protects the reported numbers from the tuning, not from the "
    "choice of the grids, which was informed by that analysis.",
    "Transfer: the benchmark reference is ~6,800 species-dereplicated genomes per fold. "
    "Production searches ~315,000 genomes, where the first 100 hits of a well-sequenced "
    "species can be one species and identity weighting collapses onto the tie set; the "
    "species-dedup variant cannot speak to this here (the panel holds one genome per "
    "species), only the density sensitivity can, and only partially.",
    "Genus-blocked folds: every number is cross-genus transfer (the source's genus is never "
    "in the reference); the clade prior is the nearest available clade, never the source's.",
    "Support expansion: E1-E4 give non-zero q outside the union (and E4 maps q = 0 to the "
    "base rate). Universe metrics include it; Stage O does not use it (union support only).",
]


if __name__ == "__main__":
    main()
