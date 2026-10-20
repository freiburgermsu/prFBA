#!/usr/bin/env python
"""Calibration check of prFBA's probabilistic weighting p(function) = k/n.

QUESTION
  When prFBA merges a selected multi-genome set R (n = |R| representatives), every
  function x of the union gets the weight p(x) = k_x / n, where k_x of the n genomes
  carry x (prFBA.tex eq:frequency; bvbrc_to_kbase_genome.create_synthetic_genome).
  Is p a calibrated probability that x is present in the organism the amplicon came
  from?  Read-only on existing artefacts: no rerun, no GPU, no network.

POPULATIONS
  PRIMARY      the controlled genus-blocked shared-reference benchmark: score_folds'
               scored cells (cells.csv.gz space=ec, outcome=scored) of every requested
               (select mode, density, forced) run. The cell's own src_genome_id is the
               truth genome; each cell is scored on its own, as score_folds does.
  SECONDARY    the as-shipped region-validation selections (picrust2_comparison/data/
               prfba_selection.json), regimes include and exclude_0.987, EC unit via
               genome_function_sets.json on compare.py's shared EC vocabulary.
               include-self is a CEILING (the source genome votes); exclude_0.987 is a
               post-hoc exclusion on truncated top-20 records (SPEC D14 rejects it).
  EXPLORATORY  top-k identity-weighted prevalence from the controlled hits files.
               NOT pre-registered; reported as exploratory only.

UNITS (--spaces)
  ec        EC numbers in panel_ec_universe.json (2,235; the scored vocabulary).
  rxn_ec    ModelSEED reactions via function_maps.ec_rxn, exactly score_folds' rxn space
            (truth = projection of the truth ECs; universe = reactions of the universe
            ECs). Reaction weight w(R) = max over member ECs reaching R of p_EC.
  rxn_role  reactions via function_maps.pgfam_rxn (PGFam role route; prFBA-only, never
            in the fair comparison). w(R) = max over member PGFams reaching R of p_PGFam.
  pgfam     PGFam families (truth and references from the five d100 fold PGFam caches).
  product   PGFam family_product strings (pgfam_products.json); empty strings and the
            generic 'hypothetical protein' variants (optionally 'FIGnnnnnn:'-prefixed,
            'conserved'/'putative'/'uncharacterized'-qualified, 'PREDICTED:') are
            removed. A proxy for the RAST function strings the real merge counts.

DENOMINATORS
  union     the union's own functions (p > 0) per scored cell -- the headline.
  universe  the full fixed vocabulary with p = 0 as a real bin (ec, rxn_ec, rxn_role).

UNCERTAINTY
  Cluster bootstrap, score_folds.cluster_bootstrap's draw verbatim (numpy
  default_rng(seed).integers(0, K, size=(B, K)) over the stratum's sorted cluster keys).
  PRIMARY unit = the holdout group (cells.csv.gz group_id == fold_plan genome_group),
  the benchmark's pre-specified convention; source genome and genus are reported as
  sensitivity CIs (ci_src, ci_genus). Source genomes nest in genera, genera in groups,
  so every scheme reweights the same per-source additive sums. The same draws serve
  every unit of a stratum, so unit contrasts are paired. The as-shipped secondary has
  no holdout groups: its primary unit is the source genus (truth_lineage), source
  genome is the sensitivity.

OUTPUTS (--outdir)
  reliability.csv  scores.json  assertions.json  README.md

USAGE
  PY=~/Documents/py_venv/bin/python
  $PY calibrate_pfunc.py --modes exact_tie,legacy,graduated,exact_tie_guards \\
      --densities 1.0,0.5,0.25 --forced --spaces ec,rxn_ec,rxn_role,pgfam,product \\
      --bootstrap 2000 --seed 1729 --eps 0.01,0.001 --outdir data/calibration/ \\
      --cache-dir <scratch dir> --workers 4

Interpreter: /home/freiburger/Documents/py_venv/bin/python (the ONLY interpreter).
"""
from __future__ import annotations

import os

# single-threaded BLAS per process: the worker count is the CPU budget
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
           "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")

import argparse            # noqa: E402
import collections         # noqa: E402
import csv                 # noqa: E402
import gzip                # noqa: E402
import hashlib             # noqa: E402
import math                # noqa: E402
import multiprocessing as mp  # noqa: E402
import pickle              # noqa: E402
import re                  # noqa: E402
import resource            # noqa: E402
import sys                 # noqa: E402
import time                # noqa: E402
from math import gcd       # noqa: E402

import numpy as np         # noqa: E402
import orjson              # noqa: E402
import pandas as pd        # noqa: E402

# ----------------------------------------------------------------------------- paths
HERE = os.path.dirname(os.path.abspath(__file__))
SR = os.path.dirname(HERE)                                   # .../shared_reference
D = os.path.join(SR, "data")
RV = os.path.dirname(SR)                                     # .../region_validation
PCD = os.path.join(RV, "picrust2_comparison", "data")
CELLS = os.path.join(D, "cells.csv.gz")
FOLD_PLAN = os.path.join(D, "fold_plan.json")
PANEL = os.path.join(D, "panel_eligible.json")
EC_UNI = os.path.join(D, "panel_ec_universe.json")
PANEL_EC = os.path.join(D, "_cache_panel_ec.json")
MAPS = os.path.join(PCD, "function_maps.json")
PRODUCTS = os.path.join(PCD, "pgfam_products.json")
SELECTION = os.path.join(PCD, "prfba_selection.json")
FUNCSETS = os.path.join(PCD, "genome_function_sets.json")
CMP_CELLS = os.path.join(PCD, "comparison_cells.csv.gz")
AS_TRUTH = os.path.join(RV, "data", "benchmark_truth.json")
FANOUT = os.path.join(RV, "data", "amplicon_md5_fanout.json")
P2_DEFAULTS = "/home/freiburger/opt/picrust2/src/picrust2-2.6.3/picrust2/default_files"

DTAG = {"1.0": "d100", "0.5": "d050", "0.25": "d025"}
ALL_SPACES = ("ec", "rxn_ec", "rxn_role", "pgfam", "product")
UNIVERSE_SPACES = ("ec", "rxn_ec", "rxn_role")
RULES = ("union", "p_ge_0.5", "p_gt_0.5", "p_eq_1")
SCHEMES = ("group", "src", "genus")        # primary first; see UNCERTAINTY
CONTRAST_STRATA = ("multi|all=all", "all|all=all", "n2|all=all", "n3|all=all",
                   "n4+|all=all")
CONTRAST_METRICS = ("murphy_reliability", "calib_slope_wls", "bss_vs_const",
                    "ece_distinct_p", "brier_p")
EXACT_MIN_POINTS = 1000
EXACT_MAX_N = 4
TOPK_K = (5, 10, 20, 50, 100)
TOPK_TAU = 0.01
STAGE0 = dict(cov_min=0.90, swed_max=0.02, family_floor=0.865)     # select_references
K2 = 200
AS_IDENT_EDGES = [0.0, 0.865, 0.945, 0.987, 0.995, 1.0001]           # compare.py
AS_IDENT_LABELS = ["<0.865", "0.865-0.945", "0.945-0.987", "0.987-0.995", ">=0.995"]

G: dict = {}            # global context, populated in the parent BEFORE forking workers


def log(msg):
    print(f"[calib {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def jload(path):
    with open(path, "rb") as fh:
        return orjson.loads(fh.read())


def rss_gb():
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6


def norm_ec(e):
    """score_folds.norm_ec verbatim: strip 'EC:'; drop partial ECs (containing '-')."""
    e = str(e).strip()
    if e.upper().startswith("EC:"):
        e = e[3:].strip()
    return None if ("-" in e or not e) else e


_FIG = re.compile(r"^fig\d+\s*:\s*")
_HYP = re.compile(r"(?:predicted:\s*)?(?:(?:conserved|putative|uncharacteri[sz]ed)\s+)*"
                  r"hypothetical\s+protein\.?")


def is_generic_hypothetical(s):
    t = _FIG.sub("", s.strip().casefold()).strip()
    return bool(_HYP.fullmatch(t))


def csr(lists):
    lens = np.fromiter((len(x) for x in lists), np.int64, len(lists))
    ptr = np.zeros(len(lists), np.int64)
    if len(lists):
        ptr[1:] = np.cumsum(lens)[:-1]
    idx = np.fromiter((j for x in lists for j in x), np.int64, int(lens.sum()))
    return ptr, idx


def in_sorted(u, T):
    """Boolean mask: elements of sorted-unique u present in sorted-unique T."""
    if len(T) == 0 or len(u) == 0:
        return np.zeros(len(u), bool)
    i = np.searchsorted(T, u)
    i[i == len(T)] = 0
    return T[i] == u


def ident_bin_as(v):
    if v is None:
        return "unknown"
    for i in range(len(AS_IDENT_LABELS)):
        if AS_IDENT_EDGES[i] <= v < AS_IDENT_EDGES[i + 1]:
            return AS_IDENT_LABELS[i]
    return "unknown"


def pseudo_fold(src):
    return int(hashlib.md5(src.encode()).hexdigest(), 16) % 5


# ============================================================================ context
def build_runs(args):
    runs = []
    for m in args.modes:
        runs.append(dict(key=f"{m}|d100", mode=m, pred_key=m, arm=f"prfba_{m}",
                         density="1.0", forced=False, dir="d100"))
    for d in args.densities:
        if d == "1.0":
            continue
        for m in args.extra_density_modes:
            runs.append(dict(key=f"{m}|{DTAG[d]}", mode=m, pred_key=m, arm=f"prfba_{m}",
                             density=d, forced=False, dir=DTAG[d]))
    if args.forced:
        for m in args.forced_modes:
            runs.append(dict(key=f"{m}_forced|d100", mode=m, pred_key=f"{m}_forced",
                             arm=f"prfba_{m}_forced", density="1.0", forced=True,
                             dir="d100_forced"))
    return runs


def load_context(args, runs):
    """Everything genome-side, as dense integer structures shared by forked workers."""
    t0 = time.time()
    A = collections.OrderedDict()
    panel = jload(PANEL)["genomes"]
    gids = sorted(panel)
    gidx = {g: i for i, g in enumerate(gids)}
    NG = len(gids)
    G.update(gids=gids, gidx=gidx, NG=NG)

    # ---- taxonomy + holdout-group codes (leakage checks, cluster labels) -------------
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
        "groups": int(G["group_code"].max() + 1), "genera": int(G["genus_code"].max() + 1),
        "species": int(G["species_code"].max() + 1),
        "taxon_ids": int(G["taxon_code"].max() + 1),
        "genera_spanning_more_than_one_group": int(sum(
            1 for n in collections.Counter(
                gen for gen, _ in set(zip(G["genus_code"].tolist(),
                                          G["group_code"].tolist()))).values() if n > 1))}
    del fp, gg

    # ---- EC universe + truth (score_folds: panel_ec_universe after norm_ec) ----------
    uni, seen = [], set()
    for e in jload(EC_UNI)["ec"]:
        n = norm_ec(e)
        if n and n not in seen:
            seen.add(n)
            uni.append(n)
    ec_index = {e: i for i, e in enumerate(uni)}
    NU = len(uni)
    pec = jload(PANEL_EC)
    E = np.zeros((NG, NU), np.uint8)
    raw_ec = {}
    n_changed = 0
    for g, v in pec.items():
        i = gidx.get(g)
        if i is None:
            continue
        # truth ECs through norm_ec, as score_folds.load_functions does
        ids = []
        for e in v["ec"]:
            ne = norm_ec(e)
            n_changed += int(ne != e)
            if ne is not None and ne in ec_index:
                ids.append(ec_index[ne])
        E[i, ids] = 1
        raw_ec[i] = frozenset(v["ec"])
    A["panel_ec_entries_changed_by_norm_ec"] = n_changed
    A["panel_genomes"] = NG
    A["panel_genomes_without_panel_ec_record"] = sum(1 for g in gids if g not in pec)
    A["ec_universe_size"] = NU
    del pec
    G.update(ec_list=uni, ec_index=ec_index, NU=NU, E=E)
    log(f"EC: universe {NU}, panel {NG} genomes ({time.time()-t0:.1f}s)")

    maps = jload(MAPS)
    # ---- rxn_ec: score_folds.load_ec_rxn verbatim (sorted keys, norm_ec, overwrite) ----
    ec_rxn = {}
    for e, rx in sorted(maps["ec_rxn"].items()):
        n = norm_ec(e)
        if n is None or n not in ec_index:
            continue
        ec_rxn[ec_index[n]] = list(rx)
    rxn_vocab = sorted({r for v in ec_rxn.values() for r in v})
    r_index = {r: i for i, r in enumerate(rxn_vocab)}
    pre = [[] for _ in rxn_vocab]
    for e, rx in sorted(ec_rxn.items()):
        for r in set(rx):
            pre[r_index[r]].append(e)
    ec_ptr, ec_idx = csr(pre)
    Rg_ec = np.maximum.reduceat(E[:, ec_idx], ec_ptr, axis=1).astype(np.uint8)
    G.update(ec_ptr=ec_ptr, ec_idx=ec_idx, Rg_ec=Rg_ec, NR_ec=len(rxn_vocab))
    A["rxn_ec_universe_size"] = len(rxn_vocab)
    log(f"rxn_ec: {len(rxn_vocab)} reactions reachable from the EC universe")

    need_pg = any(s in args.spaces for s in ("rxn_role", "pgfam", "product")) or True
    # ---- PGFam: the five d100 fold caches; conflict check over EVERY fold directory ----
    pg_vocab = {}
    pg_arr = [None] * NG
    pg_hash = [None] * NG
    conf = collections.OrderedDict()
    dirs = []
    for f in range(5):
        dirs.append(f"f{f:02d}_d100")
    for f in range(5):
        for tag in ("d100_forced", "d050", "d025"):
            dirs.append(f"f{f:02d}_{tag}")
    n_out = 0
    for dname in dirs:
        path = os.path.join(D, "prfba", dname, "_cache_pgfam_subset.json")
        if not os.path.exists(path):
            conf[dname] = "missing"
            continue
        sub = jload(path)
        nconf = nnew = nchk = 0
        for g, v in sub.items():
            i = gidx.get(g)
            if i is None:
                n_out += 1
                continue
            h = (len(v), hash(frozenset(v)))
            if pg_arr[i] is None:
                if not dname.endswith("_d100"):
                    nnew += 1          # first seen outside d100 (should not happen)
                ids = np.fromiter((pg_vocab.setdefault(s, len(pg_vocab)) for s in v),
                                  np.int64, len(v))
                pg_arr[i] = np.unique(ids).astype(np.int32)
                pg_hash[i] = h
            else:
                nchk += 1
                if pg_hash[i] != h:
                    nconf += 1
        conf[dname] = {"genomes": len(sub), "compared": nchk, "conflicts": nconf,
                       "first_seen_outside_d100": nnew}
        del sub
    A["pgfam_cache_conflicts"] = conf
    A["pgfam_cache_conflicts_total"] = sum(v["conflicts"] for v in conf.values()
                                           if isinstance(v, dict))
    A["pgfam_cache_genomes_outside_panel"] = n_out
    A["pgfam_panel_genomes_without_cache"] = sum(1 for a in pg_arr if a is None)
    for i in range(NG):
        if pg_arr[i] is None:
            pg_arr[i] = np.zeros(0, np.int32)
    inv_pg = [None] * len(pg_vocab)
    for s, i in pg_vocab.items():
        inv_pg[i] = s
    NPG = len(inv_pg)
    G.update(pg_arr=pg_arr, NPG=NPG)
    A["pgfam_vocab_panel"] = NPG
    log(f"PGFam: {NPG} families over {NG} genomes; conflicts "
        f"{A['pgfam_cache_conflicts_total']} ({time.time()-t0:.1f}s, rss {rss_gb():.1f} GB)")

    # ---- assertion: pgfam_ec projection of the PGFam caches reproduces _cache_panel_ec --
    pgfam_ec = maps["pgfam_ec"]
    has = np.zeros(NPG, bool)
    pe = {}
    for s, v in pgfam_ec.items():
        i = pg_vocab.get(s)
        if i is not None:
            has[i] = True
            pe[i] = v
    n_raw_mis = n_uni_mis = n_chk = 0
    for i in range(NG):
        if i not in raw_ec:
            continue
        n_chk += 1
        a = pg_arr[i]
        proj = set()
        for j in a[has[a]].tolist():
            proj.update(pe[j])
        if proj != raw_ec[i]:
            n_raw_mis += 1
        pu = sorted({ec_index[e] for e in proj if e in ec_index})
        if pu != np.flatnonzero(E[i]).tolist():
            n_uni_mis += 1
    A["pgfam_ec_projection"] = {"genomes_checked": n_chk,
                                "raw_set_mismatches": n_raw_mis,
                                "universe_restricted_mismatches": n_uni_mis}
    del raw_ec, pe

    # ---- rxn_role: pgfam_rxn incidence restricted to panel PGFams ------------------------
    prx = maps["pgfam_rxn"]
    role_pg = sorted(pg_vocab[s] for s in prx if s in pg_vocab)
    role_dense = np.full(NPG, -1, np.int64)
    role_dense[role_pg] = np.arange(len(role_pg))
    role_rxn_vocab = sorted({r for p in role_pg for r in prx[inv_pg[p]]})
    rr_index = {r: i for i, r in enumerate(role_rxn_vocab)}
    pre = [[] for _ in role_rxn_vocab]
    for d, p in enumerate(role_pg):
        for r in set(prx[inv_pg[p]]):
            pre[rr_index[r]].append(d)
    role_ptr, role_idx = csr(pre)
    NPR = len(role_pg)
    Pr = np.zeros((NG, NPR), np.uint8)
    for i in range(NG):
        d = role_dense[pg_arr[i]]
        Pr[i, d[d >= 0]] = 1
    Rg_role = np.zeros((NG, len(role_rxn_vocab)), np.uint8)
    for s in range(0, NG, 500):
        Rg_role[s:s + 500] = np.maximum.reduceat(Pr[s:s + 500][:, role_idx], role_ptr, axis=1)
    G.update(Pr=Pr, role_ptr=role_ptr, role_idx=role_idx, Rg_role=Rg_role,
             NR_role=len(role_rxn_vocab))
    A["rxn_role"] = {"pgfams_with_reactions_in_panel": NPR,
                     "reaction_universe_size": len(role_rxn_vocab)}
    del maps, prx
    log(f"rxn_role: {NPR} role PGFams -> {len(role_rxn_vocab)} reactions "
        f"({time.time()-t0:.1f}s)")

    # ---- product strings -----------------------------------------------------------------
    if "product" in args.spaces:
        raw = jload(PRODUCTS)
        prod_of = np.full(NPG, -1, np.int64)
        pv = {}
        n_missing = n_empty = n_hyp = 0
        hyp_strings = collections.Counter()
        for p, s in enumerate(inv_pg):
            if s not in raw:
                n_missing += 1
                continue
            st = (((raw[s] or {}).get("family_product")) or "").strip()
            if not st:
                n_empty += 1
                continue
            if is_generic_hypothetical(st):
                n_hyp += 1
                hyp_strings[st] += 1
                continue
            prod_of[p] = pv.setdefault(st, len(pv))
        del raw
        pr_arr = []
        for i in range(NG):
            q = prod_of[pg_arr[i]]
            pr_arr.append(np.unique(q[q >= 0]).astype(np.int32))
        G.update(pr_arr=pr_arr)
        A["product_unit"] = {
            "panel_pgfams": NPG, "pgfams_without_product_record": n_missing,
            "pgfams_with_empty_product": n_empty,
            "pgfams_with_generic_hypothetical_product_dropped": n_hyp,
            "distinct_generic_hypothetical_strings_dropped": len(hyp_strings),
            "most_common_dropped": hyp_strings.most_common(8),
            "distinct_product_strings_kept": len(pv),
            "rule": "strip; drop empty; drop if, after removing a leading 'FIGnnnnnn:' and "
                    "casefolding, the string is exactly '[predicted:] [conserved|putative|"
                    "uncharacterized ]*hypothetical protein[.]'; qualified strings such as "
                    "'hypothetical protein; putative signal peptide' are kept"}
        log(f"product: {len(pv)} strings kept; {n_hyp} PGFams with generic hypothetical "
            f"dropped ({time.time()-t0:.1f}s, rss {rss_gb():.1f} GB)")

    # ---- twin md5s (run_prfba_fold: fanout genomes intersected with the panel) -----------
    fan = jload(FANOUT)
    pset = set(gids)
    twin = set()
    for m, cl in fan.items():
        if len({c.split("__", 1)[0] for c in cl} & pset) > 1:
            twin.add(m)
    G["twin"] = twin
    del fan
    G["A_ctx"] = A
    log(f"context ready: {len(twin)} twin md5s ({time.time()-t0:.1f}s, rss {rss_gb():.1f} GB)")


def load_cells(runs, want_nn1=True):
    t0 = time.time()
    arms = {r["arm"] for r in runs}
    if want_nn1:
        arms.add("base_nn1")
    dens = {r["density"] for r in runs}
    usecols = ["amplicon", "md5", "region", "src_genome_id", "fold", "group_id",
               "density", "arm", "space", "outcome", "n_truth", "f1", "n_reps",
               "novelty_bin", "rank1_rank", "rank1_identity"]
    parts = []
    for ch in pd.read_csv(CELLS, dtype=str, usecols=usecols, low_memory=False,
                          chunksize=500_000, keep_default_na=False):
        ch = ch[ch.arm.isin(arms) & ch.density.isin(dens) & (ch.outcome == "scored")]
        parts.append(ch)
    df = pd.concat(parts, ignore_index=True)
    cells_ec, cells_rxn = {}, {}
    for r in runs:
        s = df[(df.arm == r["arm"]) & (df.density == r["density"])]
        e = s[s.space == "ec"].sort_values(["fold", "amplicon"]).reset_index(drop=True)
        cells_ec[r["key"]] = {f: g.reset_index(drop=True) for f, g in e.groupby("fold")}
        x = s[s.space == "rxn"]
        cells_rxn[r["key"]] = dict(zip(x.amplicon, zip(x.f1.astype(float),
                                                       x.n_truth.astype(int))))
    nn1 = df[(df.arm == "base_nn1") & (df.density == "1.0") & (df.space == "ec")]
    G["nn1_f1"] = dict(zip(nn1.amplicon, nn1.f1.astype(float)))
    G["cells_ec"] = cells_ec
    G["cells_rxn"] = cells_rxn
    log(f"cells: {sum(sum(len(v) for v in d.values()) for d in cells_ec.values())} scored "
        f"EC cells over {len(runs)} runs ({time.time()-t0:.1f}s)")


# ============================================================================ layer A
def compress(k, T, n):
    """k: per-item count of carrying reps (0..n); T: truth membership (bool, same axis)."""
    nz = k > 0
    kc = k[nz].astype(np.int64)
    m = np.bincount(kc, minlength=n + 1)
    pos = np.bincount(kc[T[nz]], minlength=n + 1)
    return m, pos


def counts_sorted(arrs, T, n):
    if n == 1:
        a = arrs[0]
        m = np.zeros(2, np.int64)
        pos = np.zeros(2, np.int64)
        m[1] = len(a)
        pos[1] = int(in_sorted(a, T).sum())
        return m, pos
    u, cnt = np.unique(np.concatenate(arrs), return_counts=True)
    y = in_sorted(u, T)
    return np.bincount(cnt, minlength=n + 1), np.bincount(cnt[y], minlength=n + 1)


def layerA_job(job):
    dirtag, fold, run_keys, spaces, limit = job
    t0 = time.time()
    gidx, E, ec_index = G["gidx"], G["E"], G["ec_index"]
    ddir = os.path.join(D, "prfba", f"f{fold:02d}_{dirtag}")
    out = {"dir": dirtag, "fold": fold, "runs": {}, "dir_asserts": {}}

    # (e) reference-side EC cache agreement
    ecsub = jload(os.path.join(ddir, "_cache_ec_subset.json"))
    mis = outside = 0
    for g, v in ecsub.items():
        i = gidx.get(g)
        if i is None:
            outside += 1
            continue
        a = sorted(ec_index.get(e, -1) for e in v)
        if a != np.flatnonzero(E[i]).tolist():
            mis += 1
    out["dir_asserts"]["ec_subset_cache"] = {"genomes": len(ecsub), "mismatches": mis,
                                             "outside_panel": outside}
    ref_set = set(ecsub)
    del ecsub

    P = jload(os.path.join(ddir, "predictions.json"))
    pred, meta = P["pred"], P["meta"]
    del P
    first = next(iter(pred.values()))
    out["dir_asserts"]["twin_md5"] = {
        "recomputed": sum(1 for m in first if m in G["twin"]),
        "meta_n_twin_md5": meta.get("n_twin_md5")}
    out["dir_asserts"]["meta_assertions"] = {
        k: meta.get("assertions", {}).get(k) for k in
        ("species_dedup_collisions", "k2_saturated_records", "reps_outside_reference",
         "n_reference_genomes_missing_function_set",
         "n_reference_genomes_missing_pgfam_cache")}

    for rk in run_keys:
        run = G["runs_by_key"][rk]
        Pm = pred[run["pred_key"]]
        sub = G["cells_ec"][rk].get(str(fold))
        if sub is None:
            continue
        if limit:
            sub = sub.iloc[:limit]
        rxn_rows = G["cells_rxn"][rk]
        Acnt = collections.Counter()
        meta_cols = collections.defaultdict(list)
        reps_flat, reps_len = [], []
        acc = {s: {"T": [], "lc": [], "c": [], "m": [], "pos": []} for s in spaces}
        diag = {s: collections.Counter() for s in spaces}
        li = 0
        for row in sub.itertuples(index=False):
            Acnt["cells_in"] += 1
            rec = Pm.get(row.md5)
            if rec is None:
                Acnt["md5_missing_in_predictions"] += 1
                continue
            reps = rec["reps"]
            n = len(reps)
            ri = [gidx.get(r, -1) for r in reps]
            si = gidx.get(row.src_genome_id, -1)
            if n == 0 or -1 in ri or si < 0:
                Acnt["reps_or_src_outside_panel_or_empty"] += 1
                continue
            Acnt["reps_outside_fold_reference"] += sum(1 for r in reps if r not in ref_set)
            # ---- leakage (cells): source never a rep; no rep shares its genus,
            # species, taxon_id or holdout group; cells' group_id == fold_plan's
            ra = np.asarray(ri, np.int64)
            Acnt["leak_source_in_reps"] += int(si in ri)
            for nm in ("genus", "species", "taxon", "group"):
                arr = G[f"{nm}_code"]
                Acnt[f"leak_rep_shares_{nm}_with_source"] += int((arr[ra] == arr[si]).any())
            if row.group_id != G["group_name"][si]:
                Acnt["group_id_ne_fold_plan_genome_group"] += 1
            if n != int(row.n_reps) or n != int(rec["n_reps"]):
                Acnt["n_reps_mismatch"] += 1
            # ---- EC (dense) ----
            k = E[ri[0]] if n == 1 else E[ri].sum(0, dtype=np.int16)
            T = E[si].astype(bool)
            union_ids = np.flatnonzero(k).tolist()
            if union_ids != sorted(ec_index.get(e, -1) for e in rec["ec"]):
                Acnt["union_ec_ne_pred_ec"] += 1
            m, pos = compress(k, T, n)
            tp, ps, ts = int(pos.sum()), int(m.sum()), int(T.sum())
            f1 = 2 * tp / (ps + ts) if (ps + ts) else 0.0
            if abs(f1 - float(row.f1)) > 1e-3:
                Acnt["f1_mismatch"] += 1
            if ts != int(row.n_truth):
                Acnt["n_truth_mismatch"] += 1
            Acnt["cells"] += 1
            res = {"ec": (ts, m, pos)}
            # ---- reactions via ec_rxn (score_folds rxn space) ----
            if "rxn_ec" in spaces:
                kR = np.maximum.reduceat(k[G["ec_idx"]], G["ec_ptr"])
                TR = G["Rg_ec"][si].astype(bool)
                mR, pR = compress(kR, TR, n)
                res["rxn_ec"] = (int(TR.sum()), mR, pR)
                rr = rxn_rows.get(row.amplicon)
                if rr is None:
                    Acnt["rxn_row_missing"] += 1
                else:
                    Acnt["rxn_cells_checked"] += 1
                    tpr, psr, tsr = int(pR.sum()), int(mR.sum()), int(TR.sum())
                    f1r = 2 * tpr / (psr + tsr) if (psr + tsr) else 0.0
                    if abs(f1r - rr[0]) > 1e-3:
                        Acnt["rxn_f1_mismatch"] += 1
                    if tsr != rr[1]:
                        Acnt["rxn_n_truth_mismatch"] += 1
                if n > 1:
                    anyc = G["Rg_ec"][ri].sum(0, dtype=np.int16)
                    u = kR > 0
                    diag["rxn_ec"]["points"] += int(u.sum())
                    diag["rxn_ec"]["max_rule_ne_genome_prevalence"] += int((anyc[u] != kR[u]).sum())
            # ---- reactions via pgfam_rxn (role route) ----
            if "rxn_role" in spaces:
                Pr = G["Pr"]
                kP = Pr[ri[0]] if n == 1 else Pr[ri].sum(0, dtype=np.int16)
                kR = np.maximum.reduceat(kP[G["role_idx"]], G["role_ptr"])
                TR = G["Rg_role"][si].astype(bool)
                mR, pR = compress(kR, TR, n)
                res["rxn_role"] = (int(TR.sum()), mR, pR)
                if n > 1:
                    anyc = G["Rg_role"][ri].sum(0, dtype=np.int16)
                    u = kR > 0
                    diag["rxn_role"]["points"] += int(u.sum())
                    diag["rxn_role"]["max_rule_ne_genome_prevalence"] += int((anyc[u] != kR[u]).sum())
            if "pgfam" in spaces:
                pa = G["pg_arr"]
                Tp = pa[si]
                m2, p2 = counts_sorted([pa[r] for r in ri], Tp, n)
                res["pgfam"] = (len(Tp), m2, p2)
            if "product" in spaces:
                pa = G["pr_arr"]
                Tp = pa[si]
                m2, p2 = counts_sorted([pa[r] for r in ri], Tp, n)
                res["product"] = (len(Tp), m2, p2)
            for s in spaces:
                ts_, m_, p_ = res[s]
                a = acc[s]
                a["T"].append(ts_)
                cs = np.flatnonzero(m_)
                a["lc"].append(np.full(len(cs), li, np.int32))
                a["c"].append(cs.astype(np.int16))
                a["m"].append(m_[cs].astype(np.int64))
                a["pos"].append(p_[cs].astype(np.int64))
            r1 = rec.get("rank1_genome")
            meta_cols["amplicon"].append(row.amplicon)
            meta_cols["md5"].append(row.md5)
            meta_cols["src"].append(si)
            meta_cols["fold"].append(int(row.fold))
            meta_cols["group_id"].append(row.group_id)
            meta_cols["region"].append(row.region)
            meta_cols["novelty_bin"].append(row.novelty_bin)
            meta_cols["rank1_rank"].append(row.rank1_rank)
            meta_cols["rank1_identity"].append(float(row.rank1_identity)
                                               if row.rank1_identity else float("nan"))
            meta_cols["rank1_gi"].append(gidx.get(r1, -1) if r1 else -1)
            meta_cols["n"].append(n)
            meta_cols["f1_cells"].append(float(row.f1))
            meta_cols["twin"].append(row.md5 in G["twin"])
            reps_flat.extend(ri)
            reps_len.append(n)
            li += 1
        sp_out = {}
        for s in spaces:
            a = acc[s]
            sp_out[s] = {
                "T": np.asarray(a["T"], np.int64),
                "lc": np.concatenate(a["lc"]) if a["lc"] else np.zeros(0, np.int32),
                "c": np.concatenate(a["c"]) if a["c"] else np.zeros(0, np.int16),
                "m": np.concatenate(a["m"]) if a["m"] else np.zeros(0, np.int64),
                "pos": np.concatenate(a["pos"]) if a["pos"] else np.zeros(0, np.int64),
                "diag": dict(diag[s])}
        out["runs"][rk] = {"meta": {k: np.asarray(v) for k, v in meta_cols.items()},
                           "reps_flat": np.asarray(reps_flat, np.int32),
                           "reps_len": np.asarray(reps_len, np.int32),
                           "spaces": sp_out, "asserts": dict(Acnt)}
    del pred
    out["seconds"] = round(time.time() - t0, 1)
    out["peak_rss_gb"] = round(rss_gb(), 2)
    return out


def merge_layerA(results, runs, spaces):
    """Concatenate the per-(dir, fold) partials of every run in fold order."""
    merged = {}
    dir_asserts = collections.OrderedDict()
    for r in sorted(results, key=lambda x: (x["dir"], x["fold"])):
        dir_asserts[f"f{r['fold']:02d}_{r['dir']}"] = r["dir_asserts"]
    for run in runs:
        rk = run["key"]
        parts = [r["runs"][rk] for r in sorted(results, key=lambda x: x["fold"])
                 if rk in r["runs"]]
        if not parts:
            continue
        meta = {k: np.concatenate([p["meta"][k] for p in parts]) for k in parts[0]["meta"]}
        off = np.cumsum([0] + [len(p["meta"]["n"]) for p in parts])
        sp = {}
        for s in spaces:
            sp[s] = {"T": np.concatenate([p["spaces"][s]["T"] for p in parts]),
                     "lc": np.concatenate([p["spaces"][s]["lc"].astype(np.int64) + off[j]
                                           for j, p in enumerate(parts)]),
                     "c": np.concatenate([p["spaces"][s]["c"] for p in parts]).astype(np.int64),
                     "m": np.concatenate([p["spaces"][s]["m"] for p in parts]),
                     "pos": np.concatenate([p["spaces"][s]["pos"] for p in parts])}
            dg = collections.Counter()
            for p in parts:
                dg.update(p["spaces"][s]["diag"])
            sp[s]["diag"] = dict(dg)
        asserts = collections.Counter()
        for p in parts:
            asserts.update(p["asserts"])
        merged[rk] = {"meta": meta,
                      "reps_flat": np.concatenate([p["reps_flat"] for p in parts]),
                      "reps_len": np.concatenate([p["reps_len"] for p in parts]),
                      "spaces": sp, "asserts": dict(asserts)}
    return merged, dir_asserts


# ============================================================================ stats
def fit_logistic(x, M, Y, iters=100):
    """Grouped logistic regression y ~ a + b*x by Newton-Raphson with step halving on
    the log-likelihood, on x centred and scaled by its M-weighted mean / sd (the raw
    Hessian is near-singular when x is uncentred), vectorised over leading axes.
    x: (J,); M, Y: (..., J) trials / successes. Returns (a, b) on the original x scale;
    NaN where fewer than two x levels carry weight."""
    M = np.asarray(M, float)
    Y = np.asarray(Y, float)
    x = np.broadcast_to(np.asarray(x, float), M.shape)
    shape = M.shape[:-1]
    tot = M.sum(-1)
    safe = np.where(tot > 0, tot, 1.0)
    xm = (M * x).sum(-1) / safe
    xv = (M * (x - xm[..., None]) ** 2).sum(-1) / safe
    xs = np.sqrt(np.where(xv > 0, xv, 1.0))
    z = (x - xm[..., None]) / xs[..., None]
    ybar = np.clip(Y.sum(-1) / safe, 1e-9, 1 - 1e-9)
    a = np.log(ybar / (1 - ybar)) * np.ones(shape)
    b = np.zeros(shape)

    def loglik(a_, b_):
        eta = a_[..., None] + b_[..., None] * z
        return (Y * eta - M * np.logaddexp(0.0, eta)).sum(-1)

    with np.errstate(all="ignore"):
        cur = loglik(a, b)
        for _ in range(iters):
            eta = a[..., None] + b[..., None] * z
            mu = 1.0 / (1.0 + np.exp(-eta))
            w = M * mu * (1 - mu)
            r = Y - M * mu
            g0 = r.sum(-1)
            g1 = (r * z).sum(-1)
            h00 = w.sum(-1)
            h01 = (w * z).sum(-1)
            h11 = (w * z * z).sum(-1)
            det = h00 * h11 - h01 ** 2
            da = np.nan_to_num((h11 * g0 - h01 * g1) / det)
            db = np.nan_to_num((-h01 * g0 + h00 * g1) / det)
            step = np.ones(shape)
            acc = np.zeros(shape, bool)
            na, nb, nl = a.copy(), b.copy(), cur.copy()
            for _h in range(50):
                ta = a + step * da
                tb = b + step * db
                tl = loglik(ta, tb)
                ok = (tl >= cur - 1e-12 * np.abs(cur)) & np.isfinite(tl) & ~acc
                na = np.where(ok, ta, na)
                nb = np.where(ok, tb, nb)
                nl = np.where(ok, tl, nl)
                acc |= ok
                if acc.all():
                    break
                step = np.where(acc, step, step / 2)
            moved = np.where(acc, np.maximum(np.abs(step * da), np.abs(step * db)), 0.0)
            a, b, cur = na, nb, nl
            if np.all(moved < 1e-10):
                break
    nlev = (M > 0).sum(-1)
    b_o = b / xs
    a_o = a - b * xm / xs
    bad = (nlev < 2) | ~np.isfinite(a_o) | ~np.isfinite(b_o)
    return np.where(bad, np.nan, a_o), np.where(bad, np.nan, b_o)


def logit_clip(p, eps):
    q = np.clip(np.asarray(p, float), eps, 1 - eps)
    return np.log(q / (1 - q))


def isotonic_fit_predict(xtr, ytr, wtr, xte):
    from sklearn.isotonic import IsotonicRegression
    if len(xtr) == 1:
        return np.full(len(xte), float(ytr[0]))
    ir = IsotonicRegression(y_min=0.0, y_max=1.0, increasing=True, out_of_bounds="clip")
    ir.fit(xtr, ytr, sample_weight=wtr)
    return ir.predict(xte)


def boot_weights(K, B, seed):
    """score_folds.cluster_bootstrap's draw verbatim -- rng = default_rng(seed);
    idx = rng.integers(0, K, size=(B, K)) over the stratum's K sorted cluster keys --
    returned as the (B, K) multiplicity matrix (W @ per-cluster sums == sums[idx].sum(1))."""
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, K, size=(B, K))
    flat = (idx + (np.arange(B, dtype=np.int64) * K)[:, None]).ravel()
    return np.bincount(flat, minlength=B * K).reshape(B, K).astype(np.float64)


def make_schemes(src_codes, B, seed, label_arrays, order=SCHEMES):
    """Cluster schemes for one stratum. src_codes: the stratum's sorted unique source
    indices (analyze_stratum's own order). label_arrays: {scheme: per-source label array
    indexed by source index, or None for the source itself}. Returns an ordered dict
    {scheme: (inv or None, W or None, K)}: inv maps each source to its sorted cluster;
    W is None when B == 0 or fewer than 2 clusters (score_folds: no interval)."""
    out = collections.OrderedDict()
    for s in order:
        lab = label_arrays.get(s)
        if lab is None:
            inv, K = None, len(src_codes)
        else:
            _, inv = np.unique(np.asarray(lab)[src_codes], return_inverse=True)
            K = int(inv.max()) + 1 if len(inv) else 0
        W = boot_weights(K, B, seed) if (B and K >= 2) else None
        out[s] = (inv, W, K)
    return out


def agg_rows(XM, inv, K):
    """Sum the rows of XM (sources) into their K clusters."""
    if inv is None:
        return XM
    order = np.argsort(inv, kind="stable")
    si = inv[order]
    starts = np.flatnonzero(np.r_[True, si[1:] != si[:-1]])
    out = np.add.reduceat(XM[order], starts, axis=0)
    assert len(starts) == K
    return out


class Cols:
    """Column builder for the per-source additive feature matrix X (Gs x F)."""

    def __init__(self, Gs):
        self.Gs = Gs
        self.blocks = []
        self.L = {}
        self.F = 0

    def add(self, name, v):
        v = np.asarray(v, float).reshape(self.Gs, -1)
        self.L[name] = slice(self.F, self.F + v.shape[1]) if v.shape[1] > 1 or \
            name.endswith("[]") else self.F
        self.blocks.append(v)
        self.F += v.shape[1]

    def matrix(self):
        return np.hstack(self.blocks) if self.blocks else np.zeros((self.Gs, 0))


def ci_of(est, reps, q=(2.5, 97.5)):
    reps = np.asarray(reps, float)
    ok = np.isfinite(reps)
    if est is None or not np.isfinite(est) or ok.sum() < 10:
        return {"est": None if est is None or not np.isfinite(est) else float(est),
                "lo": None, "hi": None}
    lo, hi = np.percentile(reps[ok], q)
    return {"est": float(est), "lo": float(lo), "hi": float(hi)}


def precompute_space(sp, n_cell, NU, eps_list):
    """Per-row and per-cell quantities of one run x space (stratum-independent)."""
    lc, c, m, pos = sp["lc"], sp["c"], sp["m"].astype(float), sp["pos"].astype(float)
    n_row = n_cell[lc]
    p = c / n_row
    Nc = len(n_cell)
    T = sp["T"].astype(float)
    pc = {"p": p, "n_row": n_row}
    pc["N"] = np.bincount(lc, m, Nc)
    pc["POS"] = np.bincount(lc, pos, Nc)
    pc["SQ"] = np.bincount(lc, pos * (1 - p) ** 2 + (m - pos) * p ** 2, Nc)
    pc["FN"] = T - pc["POS"]
    pc["T"] = T
    for e in eps_list:
        q = np.clip(p, e, 1 - e)
        pc[f"LL{e}"] = np.bincount(lc, -(pos * np.log(q) + (m - pos) * np.log(1 - q)), Nc)
    c0 = {"union": np.ones(Nc, np.int64), "p_ge_0.5": (n_cell + 1) // 2,
          "p_gt_0.5": n_cell // 2 + 1, "p_eq_1": n_cell.astype(np.int64)}
    for r in RULES:
        flag = c >= c0[r][lc]
        tp = np.bincount(lc, pos * flag, Nc)
        pr = np.bincount(lc, m * flag, Nc)
        pc[f"{r}.TP"] = tp
        pc[f"{r}.PRED"] = pr
        with np.errstate(all="ignore"):
            pc[f"{r}.F1"] = np.where(pr + T > 0, 2 * tp / (pr + T), 0.0)
            pc[f"{r}.P"] = np.where(pr > 0, tp / np.where(pr > 0, pr, 1), 0.0)
            pc[f"{r}.R"] = np.where(T > 0, tp / np.where(T > 0, T, 1), 0.0)
        pc[f"{r}.EMPTY"] = (pr == 0).astype(float)
    # bins: decile index 0..9 for (k/10,(k+1)/10], 10 for p == 1 (integer arithmetic)
    pc["dbin"] = np.where(c == n_row, 10, (10 * c + n_row - 1) // n_row - 1)
    pc["ekey"] = np.where(n_row <= EXACT_MAX_N, n_row * 100 + c, -1)
    g = np.gcd(c, n_row)
    pc["dpkey"] = (n_row // g) * 1000 + (c // g)
    return pc


DEC_LABELS = [f"({k/10:.1f},{(k+1)/10:.1f}]" for k in range(9)] + ["(0.9,1.0)", "p=1"]
DEC_LO = [k / 10 for k in range(10)] + [1.0]
DEC_HI = [(k + 1) / 10 for k in range(10)] + [1.0]


def analyze_stratum(sp, pc, cmeta, cmask, NU, eps_list, schemes, pair_id, want_rel=True,
                    label=None, keep_reps=False):
    """All metrics + reliability rows for one (run, space, stratum).
    schemes: make_schemes() output; the first scheme gives lo/hi, the others ci_<name>."""
    lc, m_all, pos_all = sp["lc"], sp["m"].astype(float), sp["pos"].astype(float)
    ci = np.flatnonzero(cmask)
    if len(ci) == 0:
        return None, []
    src_all = cmeta["src"]
    fold_all = cmeta["fold"]
    srcs, sl = np.unique(src_all[ci], return_inverse=True)
    Gs = len(srcs)
    names = list(schemes)
    for s_, (inv_, _W, K_) in schemes.items():
        if (inv_ is None and K_ != Gs) or (inv_ is not None and len(inv_) != Gs):
            raise RuntimeError(f"scheme {s_}: cluster map does not match the stratum's sources")
    Nc = len(cmask)
    c2l = np.full(Nc, -1, np.int64)
    c2l[ci] = sl
    rows = np.flatnonzero(cmask[lc])
    if len(rows) == 0:
        return None, []
    rsl = c2l[lc[rows]]
    m = m_all[rows]
    pos = pos_all[rows]
    p = pc["p"][rows]
    rfold = fold_all[lc[rows]]
    cfold = fold_all[ci]

    def per_src(v):
        return np.bincount(sl, weights=v[ci], minlength=Gs)

    def per_src_bin(binidx, K, w):
        return np.bincount(rsl * K + binidx, weights=w, minlength=Gs * K).reshape(Gs, K)

    X = Cols(Gs)
    X.add("NCELL", np.bincount(sl, minlength=Gs))
    for k in ("N", "POS", "SQ", "FN", "T"):
        X.add(k, per_src(pc[k]))
    for e in eps_list:
        X.add(f"LL{e}", per_src(pc[f"LL{e}"]))
    for r in RULES:
        for f in ("TP", "PRED", "F1", "P", "R", "EMPTY"):
            X.add(f"{r}.{f}", per_src(pc[f"{r}.{f}"]))
    # p = 1 from all n >= 2 tied references vs p = 1 from a single reference (paired
    # within a stratum that holds both, e.g. 'all')
    is1 = pc["p"][rows] >= 1.0
    nr = pc["n_row"][rows]
    for nm, mk_ in (("P1M", is1 & (nr >= 2)), ("P1S", is1 & (nr == 1))):
        X.add(f"{nm}_M", np.bincount(rsl, m * mk_, Gs))
        X.add(f"{nm}_POS", np.bincount(rsl, pos * mk_, Gs))
    # ---- primary partition: exact (n, c) for n <= 4 with >= 1000 points, deciles otherwise
    ek = pc["ekey"][rows]
    ekeys, einv = np.unique(np.where(ek >= 0, ek, -1), return_inverse=True)
    emass = np.bincount(einv, weights=m)
    elig = np.array([(k >= 0) and (emass[j] >= EXACT_MIN_POINTS) for j, k in enumerate(ekeys)])
    ex_keys = [int(k) for j, k in enumerate(ekeys) if elig[j]]
    ex_pos = {k: i for i, k in enumerate(ex_keys)}
    Ne = len(ex_keys)
    is_ex = elig[einv]
    dbin = pc["dbin"][rows]
    pbin = np.where(is_ex, np.array([ex_pos.get(int(k), 0) for k in ekeys])[einv], Ne + dbin)
    KP = Ne + 11
    X.add("RB_M[]", per_src_bin(pbin, KP, m))
    X.add("RB_POS[]", per_src_bin(pbin, KP, pos))
    X.add("RB_SP[]", per_src_bin(pbin, KP, p * m))
    X.add("DB_M[]", per_src_bin(dbin, 11, m))
    X.add("DB_POS[]", per_src_bin(dbin, 11, pos))
    X.add("DB_SP[]", per_src_bin(dbin, 11, p * m))
    dpk = pc["dpkey"][rows]
    dkeys, dinv = np.unique(dpk, return_inverse=True)
    J = len(dkeys)
    pj = np.array([(k % 1000) / (k // 1000) for k in dkeys], float)
    X.add("DP_M[]", per_src_bin(dinv, J, m))
    X.add("DP_POS[]", per_src_bin(dinv, J, pos))

    # ---- out-of-fold recalibration (4 folds train, 5th test), within this stratum ----
    oof_ok = len(np.unique(cfold)) >= 2
    if oof_ok:
        r_iso = np.full(len(rows), np.nan)
        r_pl = np.full(len(rows), np.nan)
        r_c = np.full(len(rows), np.nan)
        if NU:
            r_iso_v = np.full(len(rows), np.nan)
            r_pl_v = np.full(len(rows), np.nan)
            r0_iso = np.full(len(ci), np.nan)
            r0_pl = np.full(len(ci), np.nan)
        M0_cell = (NU - pc["N"][ci]) if NU else None
        FN_cell = pc["FN"][ci]
        for f in np.unique(cfold):
            tr = rfold != f
            te = rfold == f
            if not tr.any():
                continue
            Mj = np.bincount(dinv[tr], weights=m[tr], minlength=J)
            Pj = np.bincount(dinv[tr], weights=pos[tr], minlength=J)
            keep = Mj > 0
            xs, ys, ws = pj[keep], Pj[keep] / Mj[keep], Mj[keep]
            r_iso[te] = isotonic_fit_predict(xs, ys, ws, p[te])
            a, b = fit_logistic(logit_clip(xs, eps_list[0]), ws, Pj[keep])
            if not np.isfinite(a):
                yb = min(max(Pj.sum() / Mj.sum(), 1e-6), 1 - 1e-6)
                a, b = math.log(yb / (1 - yb)), 0.0
            r_pl[te] = 1 / (1 + np.exp(-(a + b * logit_clip(p[te], eps_list[0]))))
            r_c[te] = Pj.sum() / Mj.sum()
            if NU:
                ctr = cfold != f
                cte = cfold == f
                M0 = M0_cell[ctr].sum()
                P0 = FN_cell[ctr].sum()
                xs_v = np.concatenate([[0.0], xs])
                ys_v = np.concatenate([[P0 / M0], ys])
                ws_v = np.concatenate([[M0], ws])
                r_iso_v[te] = isotonic_fit_predict(xs_v, ys_v, ws_v, p[te])
                r0_iso[cte] = isotonic_fit_predict(xs_v, ys_v, ws_v, np.zeros(1))[0]
                av, bv = fit_logistic(logit_clip(xs_v, eps_list[0]), ws_v,
                                      np.concatenate([[P0], Pj[keep]]))
                r_pl_v[te] = 1 / (1 + np.exp(-(av + bv * logit_clip(p[te], eps_list[0]))))
                r0_pl[cte] = 1 / (1 + np.exp(-(av + bv * logit_clip(0.0, eps_list[0]))))
        if np.isnan(r_iso).any():
            oof_ok = False
    if oof_ok:
        def sqrow(r):
            return pos * (1 - r) ** 2 + (m - pos) * r ** 2
        X.add("SQ_iso_u", np.bincount(rsl, sqrow(r_iso), Gs))
        X.add("SQ_platt_u", np.bincount(rsl, sqrow(r_pl), Gs))
        X.add("SQ_coof_u", np.bincount(rsl, sqrow(r_c), Gs))
        if NU:
            def sq0(r0):
                return FN_cell * (1 - r0) ** 2 + (M0_cell - FN_cell) * r0 ** 2
            X.add("SQ_iso_v", np.bincount(rsl, sqrow(r_iso_v), Gs) +
                  np.bincount(sl, sq0(r0_iso), Gs))
            X.add("SQ_platt_v", np.bincount(rsl, sqrow(r_pl_v), Gs) +
                  np.bincount(sl, sq0(r0_pl), Gs))
            X.add("SQ_coof_v", np.bincount(rsl, sqrow(r_c), Gs) +
                  np.bincount(sl, FN_cell, Gs))

    XM = X.matrix()
    S = XM.sum(0)
    SB = collections.OrderedDict(
        (s_, (W_ @ agg_rows(XM, inv_, K_)) if W_ is not None else None)
        for s_, (inv_, W_, K_) in schemes.items())
    L = X.L

    def derive(S):
        g = lambda k: S[..., L[k]]  # noqa: E731
        o = {}
        N, POS, SQ, FN, NC, T = g("N"), g("POS"), g("SQ"), g("FN"), g("NCELL"), g("T")
        with np.errstate(all="ignore"):
            c = POS / N
            u = {}
            u["c_pooled_precision"] = c
            u["brier_p"] = SQ / N
            u["brier_union"] = (N - POS) / N
            u["brier_const_c"] = c * (1 - c)
            u["d_brier_p_minus_union"] = u["brier_p"] - u["brier_union"]
            u["d_brier_p_minus_const"] = u["brier_p"] - u["brier_const_c"]
            u["bss_vs_const"] = 1 - u["brier_p"] / u["brier_const_c"]
            for e in eps_list:
                u[f"logloss_eps{e}"] = g(f"LL{e}") / N
            RBM, RBP, RBS = g("RB_M[]"), g("RB_POS[]"), g("RB_SP[]")
            DBM, DBP, DBS = g("DB_M[]"), g("DB_POS[]"), g("DB_SP[]")
            DPM, DPP = g("DP_M[]"), g("DP_POS[]")
            u["ece_primary_bins"] = np.abs(RBS - RBP).sum(-1) / N
            u["ece_decile_bins"] = np.abs(DBS - DBP).sum(-1) / N
            u["ece_distinct_p"] = np.abs(pj * DPM - DPP).sum(-1) / N
            oj = np.where(DPM > 0, DPP / np.where(DPM > 0, DPM, 1), 0.0)
            rel = (DPM * (pj - oj) ** 2).sum(-1) / N
            res = (DPM * (oj - c[..., None]) ** 2).sum(-1) / N
            unc = c * (1 - c)
            u["murphy_reliability"] = rel
            u["murphy_resolution"] = res
            u["murphy_uncertainty"] = unc
            u["murphy_identity_residual"] = u["brier_p"] - (rel - res + unc)
            # probability-scale calibration line: weighted least squares of observed
            # presence on p over the distinct-p bins, weights = points (identical to the
            # point-level OLS slope of y on p); needs no clipping constant
            pbar = (DPM * pj).sum(-1) / N
            sxx = (DPM * pj * pj).sum(-1) - N * pbar ** 2
            sxy = (DPP * pj).sum(-1) - pbar * POS
            u["calib_slope_wls"] = sxy / sxx
            u["calib_intercept_wls"] = c - u["calib_slope_wls"] * pbar
            p1m, p1s = g("P1M_POS") / g("P1M_M"), g("P1S_POS") / g("P1S_M")
            u["obs_p1_multi"] = p1m
            u["obs_p1_single"] = p1s
            u["obs_p1_multi_minus_single"] = p1m - p1s
            for e in eps_list:
                a, b = fit_logistic(logit_clip(pj, e), DPM, DPP)
                u[f"calib_intercept_eps{e}"] = a
                u[f"calib_slope_eps{e}"] = b
            if "SQ_iso_u" in L:
                for nm, key in (("isotonic", "SQ_iso_u"), ("platt", "SQ_platt_u"),
                                ("const_c_oof", "SQ_coof_u")):
                    bb = g(key) / N
                    u[f"oof_brier_{nm}"] = bb
                    u[f"oof_gap_p_minus_{nm}"] = u["brier_p"] - bb
                u["bss_vs_const_oof"] = 1 - u["brier_p"] / u["oof_brier_const_c_oof"]
            o["union"] = u
            if NU:
                v = {}
                Nv = NC * NU
                M0 = Nv - N
                v["brier_p"] = (SQ + FN) / Nv
                v["brier_union"] = (N - POS + FN) / Nv
                v["brier_const_c"] = (POS * (1 - c) ** 2 + (N - POS) * c ** 2 + FN) / Nv
                v["d_brier_p_minus_union"] = v["brier_p"] - v["brier_union"]
                v["d_brier_p_minus_const"] = v["brier_p"] - v["brier_const_c"]
                v["bss_vs_const"] = 1 - v["brier_p"] / v["brier_const_c"]
                for e in eps_list:
                    v[f"logloss_eps{e}"] = (g(f"LL{e}") - FN * math.log(e) -
                                            (M0 - FN) * math.log(1 - e)) / Nv
                v["ece_primary_bins"] = (np.abs(RBS - RBP).sum(-1) + FN) / Nv
                v["ece_distinct_p"] = (np.abs(pj * DPM - DPP).sum(-1) + FN) / Nv
                ob = (POS + FN) / Nv
                o0 = FN / M0
                relv = ((DPM * (pj - oj) ** 2).sum(-1) + M0 * o0 ** 2) / Nv
                resv = ((DPM * (oj - ob[..., None]) ** 2).sum(-1) + M0 * (o0 - ob) ** 2) / Nv
                v["murphy_reliability"] = relv
                v["murphy_resolution"] = resv
                v["murphy_uncertainty"] = ob * (1 - ob)
                v["murphy_identity_residual"] = v["brier_p"] - (relv - resv + ob * (1 - ob))
                v["obs_freq_p0"] = o0
                pbv = (DPM * pj).sum(-1) / Nv            # the p = 0 bin adds mass, not p
                sxxv = (DPM * pj * pj).sum(-1) - Nv * pbv ** 2
                sxyv = (DPP * pj).sum(-1) - pbv * (POS + FN)
                v["calib_slope_wls"] = sxyv / sxxv
                v["calib_intercept_wls"] = ob - v["calib_slope_wls"] * pbv
                for e in eps_list:
                    xv = np.concatenate([[logit_clip(0.0, e)], logit_clip(pj, e)])
                    Mv = np.concatenate([M0[..., None], DPM], axis=-1)
                    Yv = np.concatenate([FN[..., None], DPP], axis=-1)
                    a, b = fit_logistic(xv, Mv, Yv)
                    v[f"calib_intercept_eps{e}"] = a
                    v[f"calib_slope_eps{e}"] = b
                if "SQ_iso_v" in L:
                    for nm, key in (("isotonic", "SQ_iso_v"), ("platt", "SQ_platt_v"),
                                    ("const_c_oof", "SQ_coof_v")):
                        bb = g(key) / Nv
                        v[f"oof_brier_{nm}"] = bb
                        v[f"oof_gap_p_minus_{nm}"] = v["brier_p"] - bb
                    v["bss_vs_const_oof"] = 1 - v["brier_p"] / v["oof_brier_const_c_oof"]
                o["universe"] = v
            bi = {}
            f1u_mi = 2 * g("union.TP") / (g("union.PRED") + T)
            f1u_ma = g("union.F1") / NC
            for r in RULES:
                TP, PR = g(f"{r}.TP"), g(f"{r}.PRED")
                d = {"micro_P": TP / PR, "micro_R": TP / T, "micro_F1": 2 * TP / (PR + T),
                     "macro_P": g(f"{r}.P") / NC, "macro_R": g(f"{r}.R") / NC,
                     "macro_F1": g(f"{r}.F1") / NC,
                     "frac_cells_empty_prediction": g(f"{r}.EMPTY") / NC}
                d["d_micro_F1_vs_union"] = d["micro_F1"] - f1u_mi
                d["d_macro_F1_vs_union"] = d["macro_F1"] - f1u_ma
                bi[r] = d
            o["binary"] = bi
        return o

    pt = derive(S)
    bts = [derive(SB[s_]) if SB[s_] is not None else None for s_ in names]
    res = pack_schemes(pt, bts, names)
    if keep_reps:
        res["_reps"] = {"est": {k: float(pt["union"][k]) for k in CONTRAST_METRICS},
                        "reps": {s_: ({k: np.asarray(b["union"][k], float)
                                       for k in CONTRAST_METRICS} if b is not None else None)
                                 for s_, b in zip(names, bts)}}
    # never-bootstrapped bookkeeping
    res["ci_unit"] = names[0] if names else None
    res["clusters"] = {}
    for s_, (inv_, W_, K_) in schemes.items():
        nsrc_pts = XM[:, L["N"]]
        nsrc_cells = XM[:, L["NCELL"]]
        if inv_ is not None:
            nsrc_pts = np.bincount(inv_, nsrc_pts, K_)
            nsrc_cells = np.bincount(inv_, nsrc_cells, K_)
        res["clusters"][s_] = {
            "n": int(K_), "ci": W_ is not None,
            "largest_share_points": float(nsrc_pts.max() / nsrc_pts.sum()) if K_ else None,
            "largest_share_cells": float(nsrc_cells.max() / nsrc_cells.sum()) if K_ else None}
    pairs = pair_id[ci]
    res["n_cells"] = int(len(ci))
    res["n_src"] = int(Gs)
    res["n_distinct_pairs"] = int(len(np.unique(pairs)))
    res["n_points_union"] = int(S[L["N"]])
    if NU:
        res["n_points_universe"] = int(S[L["NCELL"]] * NU)
    res["n_folds"] = int(len(np.unique(cfold)))
    res["n_twin_md5_cells"] = int(cmeta["twin"][ci].sum())
    res["oof_recalibration"] = bool(oof_ok)
    res["n_distinct_p"] = int(J)

    rel_rows = []
    if want_rel:
        cell_of_row = lc[rows]
        pair_row = pair_id[cell_of_row]
        src_row = src_all[cell_of_row]

        def bin_counts(binidx, K):
            out = []
            for arr in (cell_of_row, src_row, pair_row):
                key = np.unique(binidx.astype(np.int64) * (int(arr.max()) + 1) + arr)
                out.append(np.bincount(key // (int(arr.max()) + 1), minlength=K))
            return out

        def emit(scheme, labels, lo, hi, Mk, Pk, Sk, binidx, K, min_pts=0):
            nc, ns, npairs = bin_counts(binidx, K)
            Mb = S[L[Mk]]
            Pb = S[L[Pk]]
            Sb = S[L[Sk]] if Sk else None
            obs = {}
            for s_ in names:
                if SB[s_] is not None:
                    with np.errstate(all="ignore"):
                        obs[s_] = SB[s_][:, L[Pk]] / SB[s_][:, L[Mk]]
            for j in range(K):
                if Mb[j] <= 0 or Mb[j] < min_pts:
                    continue
                rel_rows.append({**label, "denominator": "union", "bin_scheme": scheme,
                                 "bin": labels[j], "p_lo": lo[j], "p_hi": hi[j],
                                 "n_points": int(Mb[j]), "n_present": int(Pb[j]),
                                 "n_cells": int(nc[j]), "n_src": int(ns[j]),
                                 "n_distinct_pairs": int(npairs[j]),
                                 "mean_p": (Sb[j] / Mb[j]) if Sb is not None else pj[j],
                                 "obs_freq": Pb[j] / Mb[j],
                                 **ci_cols(Pb[j] / Mb[j],
                                           {s_: o[:, j] for s_, o in obs.items()}, names)})

        ex_labels = [f"n={k // 100},p={k % 100}/{k // 100}" for k in ex_keys]
        emit("exact_or_decile", ex_labels + DEC_LABELS,
             [(k % 100) / (k // 100) for k in ex_keys] + DEC_LO,
             [(k % 100) / (k // 100) for k in ex_keys] + DEC_HI,
             "RB_M[]", "RB_POS[]", "RB_SP[]", pbin, KP)
        emit("decile", DEC_LABELS, DEC_LO, DEC_HI, "DB_M[]", "DB_POS[]", "DB_SP[]", dbin, 11)
        dlabels = [f"p={k % 1000}/{k // 1000}" for k in dkeys]
        emit("distinct_p", dlabels, list(pj), list(pj), "DP_M[]", "DP_POS[]", None, dinv, J,
             min_pts=EXACT_MIN_POINTS)
        if NU:
            Nv = S[L["NCELL"]] * NU
            M0 = Nv - S[L["N"]]
            FN = S[L["FN"]]
            ob0 = {}
            for s_ in names:
                if SB[s_] is not None:
                    B_ = SB[s_]
                    with np.errstate(all="ignore"):
                        ob0[s_] = B_[:, L["FN"]] / (B_[:, L["NCELL"]] * NU - B_[:, L["N"]])
            rel_rows.append({**label, "denominator": "universe", "bin_scheme": "p0",
                             "bin": "p=0 (outside union)", "p_lo": 0.0, "p_hi": 0.0,
                             "n_points": int(M0), "n_present": int(FN),
                             "n_cells": int(len(ci)), "n_src": int(Gs),
                             "n_distinct_pairs": int(len(np.unique(pairs))),
                             "mean_p": 0.0, "obs_freq": FN / M0,
                             **ci_cols(FN / M0, ob0, names)})
    return res, rel_rows


def pack_schemes(pt, bts, names):
    """{est, lo, hi (first scheme), ci_<other scheme>: [lo, hi]} for every leaf."""
    if isinstance(pt, dict):
        return {k: pack_schemes(pt[k], [None if b is None else b[k] for b in bts], names)
                for k in pt}
    a = float(pt)
    out = {"est": a if np.isfinite(a) else None}
    for j, (s_, b) in enumerate(zip(names, bts)):
        c = ci_of(a, b) if b is not None else {"lo": None, "hi": None}
        if j == 0:
            out["lo"], out["hi"] = c["lo"], c["hi"]
        else:
            out[f"ci_{s_}"] = [c["lo"], c["hi"]]
    return out


def ci_cols(est, reps_by_scheme, names):
    """reliability.csv CI columns: ci_lo/ci_hi for the primary unit, then per scheme."""
    out = {"ci_unit": names[0] if names else None}
    for j, s_ in enumerate(names):
        r = reps_by_scheme.get(s_)
        c = ci_of(est, r) if r is not None else {"lo": None, "hi": None}
        if j == 0:
            out["ci_lo"], out["ci_hi"] = c["lo"], c["hi"]
        else:
            out[f"ci_lo_{s_}"], out[f"ci_hi_{s_}"] = c["lo"], c["hi"]
    return out


def unit_contrasts(per_space, spaces):
    """Paired unit differences a - b from the shared bootstrap draws of one stratum."""
    out = {}
    sp = [s for s in spaces if per_space.get(s)]
    for i, a in enumerate(sp):
        for b in sp[i + 1:]:
            A_, B_ = per_space[a], per_space[b]
            names = list(A_["reps"])
            d = {}
            for k in CONTRAST_METRICS:
                est = A_["est"][k] - B_["est"][k]
                bts = [None if (A_["reps"][s] is None or B_["reps"][s] is None) else
                       A_["reps"][s][k] - B_["reps"][s][k] for s in names]
                d[k] = pack_schemes(est, bts, names)
            out[f"{a}-{b}"] = d
    return out


def strata_for(cmeta, kinds, cross_pops=("all", "single", "multi")):
    n = cmeta["n"]
    pops = [("all", n >= 1), ("single", n == 1), ("multi", n >= 2), ("n2", n == 2),
            ("n3", n == 3), ("n4+", n >= 4)]
    out = []
    for pname, pm in pops:
        out.append((pname, "all", "all", pm))
        if pname in cross_pops:
            for kind in kinds:
                vals = cmeta[kind]
                for v in sorted(set(vals.tolist()), key=str):
                    mk = pm & (vals == v)
                    if mk.any():
                        out.append((pname, kind, str(v), mk))
    return out


def pair_ids(cmeta, reps_flat, reps_len):
    off = np.concatenate([[0], np.cumsum(reps_len)])
    d = {}
    out = np.empty(len(reps_len), np.int64)
    src = cmeta["src"]
    for i in range(len(reps_len)):
        key = (int(src[i]), tuple(sorted(reps_flat[off[i]:off[i + 1]].tolist())))
        out[i] = d.setdefault(key, len(d))
    return out


def twin_free_first(cmeta):
    """Per-md5 population with twin md5s removed (C2 convention): one cell per md5."""
    keep = np.zeros(len(cmeta["n"]), bool)
    seen = set()
    for i, (m, t) in enumerate(zip(cmeta["md5"], cmeta["twin"])):
        if t or m in seen:
            continue
        seen.add(m)
        keep[i] = True
    return keep


def layerB_run(job):
    rk, spaces, B, seed, eps_list = job
    t0 = time.time()
    R = G["A_merged"][rk]
    run = G["runs_by_key"][rk]
    cmeta = R["meta"]
    n_cell = cmeta["n"].astype(np.int64)
    pid = pair_ids(cmeta, R["reps_flat"], R["reps_len"])
    tf = twin_free_first(cmeta)
    strata = strata_for(cmeta, ("novelty_bin", "rank1_rank", "region"))
    for pname, pm in (("all", n_cell >= 1), ("single", n_cell == 1), ("multi", n_cell >= 2)):
        strata.append((pname, "population", "twin_free_per_md5", pm & tf))
    pcs = {s: precompute_space(R["spaces"][s], n_cell, nu_of(s), eps_list) for s in spaces}
    scores = {s: {} for s in spaces}
    scores["_unit_contrasts"] = {}
    rel = []
    labs = {"group": G["group_code"], "src": None, "genus": G["genus_code"]}
    for (pname, kind, val, mk) in strata:
        skey = f"{pname}|{kind}={val}"
        schemes = make_schemes(np.unique(cmeta["src"][mk]), B, seed, labs)
        keep = skey in CONTRAST_STRATA
        reps = {}
        for s in spaces:
            label = {"analysis": "primary_controlled", "run": rk, "mode": run["mode"],
                     "density": run["density"], "forced": run["forced"], "space": s,
                     "pop": pname, "stratum_kind": kind, "stratum_value": val}
            res, rows = analyze_stratum(R["spaces"][s], pcs[s], cmeta, mk, nu_of(s),
                                        eps_list, schemes, pid, label=label,
                                        keep_reps=keep)
            if res is None:
                continue
            if keep:
                reps[s] = res.pop("_reps")
            scores[s][skey] = res
            rel.extend(rows)
        if keep and len(reps) >= 2:
            scores["_unit_contrasts"][skey] = unit_contrasts(reps, spaces)
        del schemes, reps
    for s in spaces:
        scores[s]["_diag"] = R["spaces"][s]["diag"]
    return {"run": rk, "scores": scores, "rel": rel, "seconds": round(time.time() - t0, 1),
            "peak_rss_gb": round(rss_gb(), 2)}


def nu_of(space):
    return {"ec": G.get("NU"), "rxn_ec": G.get("NR_ec"), "rxn_role": G.get("NR_role"),
            "as_ec": G.get("AS_NU")}.get(space)


# ============================================================================ secondary
def secondary_prepare(args):
    """As-shipped region validation, EC unit, compare.py's shared EC vocabulary."""
    import ijson
    t0 = time.time()
    A = collections.OrderedDict()
    sel = jload(SELECTION)
    truth = jload(AS_TRUTH)
    keys = sorted(k for k in sel if k != "_meta")
    need = set()
    for k in keys:
        need.add(truth[k]["src_genome_id"])
        for reg in ("include", "exclude_0.987"):
            need.update((sel[k].get(reg) or {}).get("reps") or [])
    maps = jload(MAPS)
    bvbrc_vocab = {e for v in maps["pgfam_ec"].values() for e in v}
    del maps
    ref_vocab = set()
    for dom in ("bacteria", "archaea"):
        with gzip.open(os.path.join(P2_DEFAULTS, dom, "ec.txt.gz"), "rt") as fh:
            ref_vocab.update(fh.readline().rstrip("\n").split("\t")[1:])
    shared = sorted(bvbrc_vocab & ref_vocab)
    sh_index = {e: i for i, e in enumerate(shared)}
    A["shared_ec_size"] = len(shared)
    A["shared_ec_size_expected_comparison_summary"] = 2132
    gl = sorted(need)
    gi = {g: i for i, g in enumerate(gl)}
    M = np.zeros((len(gl), len(shared)), np.uint8)
    found = np.zeros(len(gl), bool)
    with open(FUNCSETS, "rb") as fh:
        for gid, rec in ijson.kvitems(fh, "", use_float=True):
            i = gi.get(gid)
            if i is None or not isinstance(rec, dict):
                continue
            ids = [sh_index[e] for e in rec.get("ec", ()) if e in sh_index]
            M[i, ids] = 1
            found[i] = True
    A["genomes_needed"] = len(gl)
    A["genomes_missing_from_function_sets"] = int((~found).sum())
    # cluster labels: the source's genus (truth_lineage), else the source itself
    gen_of_src = {}
    for k in keys:
        t = truth[k]
        gen = ((t.get("truth_lineage") or {}).get("Genus") or "").strip()
        gen_of_src[t["src_genome_id"]] = gen or f"src:{t['src_genome_id']}"
    glab = [gen_of_src.get(g, f"ref:{g}") for g in gl]
    u = {v: i for i, v in enumerate(sorted(set(glab)))}
    sec_genus_code = np.array([u[v] for v in glab], np.int64)
    A["sources_without_truth_genus"] = sum(1 for v in gen_of_src.values()
                                           if v.startswith("src:"))
    log(f"secondary: function sets streamed for {len(gl)} genomes "
        f"({time.time()-t0:.1f}s, rss {rss_gb():.1f} GB)")
    cmp = pd.read_csv(CMP_CELLS, dtype=str, keep_default_na=False,
                      usecols=["amplicon", "method", "space", "outcome", "n_truth", "f1",
                               "n_reps"])
    cmp = cmp[(cmp.space == "ec") & cmp.method.isin(["prfba_include", "prfba_exclude"])]
    cmpd = {(r.method, r.amplicon): (r.outcome, r.f1, r.n_truth, r.n_reps)
            for r in cmp.itertuples(index=False)}
    out = {}
    for reg, meth in (("include", "prfba_include"), ("exclude_0.987", "prfba_exclude")):
        Acnt = collections.Counter()
        meta = collections.defaultdict(list)
        acc = {"T": [], "lc": [], "c": [], "m": [], "pos": []}
        reps_flat, reps_len = [], []
        li = 0
        for k in keys:
            t = truth[k]
            src = t["src_genome_id"]
            r = sel[k].get(reg) or {}
            reps = r.get("reps") or []
            n = len(reps)
            si = gi[src]
            T = M[si].astype(bool)
            cm = cmpd.get((meth, k))
            if n == 0 or not T.any():
                Acnt["not_answered_or_truth_empty"] += 1
                if cm is not None and cm[0] == "ok":
                    Acnt["cmp_ok_but_not_scored_here"] += 1
                continue
            ri = [gi[g] for g in reps]
            kk = M[ri[0]] if n == 1 else M[ri].sum(0, dtype=np.int16)
            if not (kk > 0).any():
                Acnt["pred_empty_in_shared_ec"] += 1
                if cm is not None and cm[0] == "ok":
                    Acnt["cmp_ok_but_not_scored_here"] += 1
                continue
            m_, p_ = compress(kk, T, n)
            tp, ps, ts = int(p_.sum()), int(m_.sum()), int(T.sum())
            f1 = 2 * tp / (ps + ts)
            Acnt["cells"] += 1
            if cm is None or cm[0] != "ok":
                Acnt["scored_here_but_not_cmp_ok"] += 1
            else:
                if abs(f1 - float(cm[1])) > 1e-3:
                    Acnt["f1_mismatch_vs_comparison_cells"] += 1
                if ts != int(cm[2]):
                    Acnt["n_truth_mismatch_vs_comparison_cells"] += 1
                if n != int(cm[3]):
                    Acnt["n_reps_mismatch_vs_comparison_cells"] += 1
            cs = np.flatnonzero(m_)
            acc["T"].append(ts)
            acc["lc"].append(np.full(len(cs), li, np.int64))
            acc["c"].append(cs.astype(np.int64))
            acc["m"].append(m_[cs])
            acc["pos"].append(p_[cs])
            meta["amplicon"].append(k)
            meta["md5"].append(k)
            meta["src"].append(si)
            meta["fold"].append(pseudo_fold(src))
            meta["region"].append(t["region"])
            meta["ident_bin"].append(ident_bin_as(sel[k].get("best_nonself_identity")))
            meta["src_in_reps"].append(str(src in reps))
            meta["n"].append(n)
            meta["twin"].append(False)
            reps_flat.extend(ri)
            reps_len.append(n)
            li += 1
        n_ok = sum(1 for (mm, _), v in cmpd.items() if mm == meth and v[0] == "ok")
        Acnt["comparison_cells_ok_rows"] = n_ok
        sp = {"T": np.asarray(acc["T"], np.int64), "lc": np.concatenate(acc["lc"]),
              "c": np.concatenate(acc["c"]), "m": np.concatenate(acc["m"]),
              "pos": np.concatenate(acc["pos"]), "diag": {}}
        out[reg] = {"meta": {k: np.asarray(v) for k, v in meta.items()},
                    "reps_flat": np.asarray(reps_flat, np.int32),
                    "reps_len": np.asarray(reps_len, np.int32),
                    "spaces": {"as_ec": sp}, "asserts": dict(Acnt)}
    A["seconds"] = round(time.time() - t0, 1)
    A["peak_rss_gb"] = round(rss_gb(), 2)
    return {"runs": out, "asserts": A, "AS_NU": len(shared), "genus_code": sec_genus_code}


def secondary_layerB(job):
    reg, B, seed, eps_list = job
    t0 = time.time()
    R = G["SEC"]["runs"][reg]
    cmeta = R["meta"]
    n_cell = cmeta["n"].astype(np.int64)
    pid = pair_ids(cmeta, R["reps_flat"], R["reps_len"])
    strata = strata_for(cmeta, ("region", "ident_bin", "src_in_reps"))
    NU = G["SEC"]["AS_NU"]
    pc = precompute_space(R["spaces"]["as_ec"], n_cell, NU, eps_list)
    scores, rel = {}, []
    labs = {"genus": G["SEC"]["genus_code"], "src": None}
    for (pname, kind, val, mk) in strata:
        schemes = make_schemes(np.unique(cmeta["src"][mk]), B, seed, labs,
                               order=("genus", "src"))
        label = {"analysis": "secondary_as_shipped", "run": f"as_shipped|{reg}",
                 "mode": reg, "density": "as_shipped", "forced": False, "space": "ec",
                 "pop": pname, "stratum_kind": kind, "stratum_value": val}
        res, rows = analyze_stratum(R["spaces"]["as_ec"], pc, cmeta, mk, NU, eps_list,
                                    schemes, pid, label=label)
        if res is not None:
            scores[f"{pname}|{kind}={val}"] = res
            rel.extend(rows)
    return {"regime": reg, "scores": scores, "rel": rel,
            "seconds": round(time.time() - t0, 1)}


# ============================================================================ top-k
def stage0_ranked(rec):
    """select_references stage 0 verbatim, then identity ranking with genome dedup."""
    top = rec.get("top20") or []
    L = rec.get("asv_len") or 1
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
        ranked.append((g, float(h["identity"])))
    Lq = max(1, int(rec.get("asv_len") or 1))
    qual = [h for h in top if (h.get("aligned_len") or 0) / Lq >= STAGE0["cov_min"]]
    last_q = float(qual[-1]["identity"]) if qual else None
    saturated = bool(len(top) >= K2 and qual and (last_q or 0.0) >= STAGE0["family_floor"])
    return ranked, {"n_hits": len(top), "n_survivors": len(ranked), "saturated": saturated,
                    "last_qualified_identity": last_q}


def topk_job(job):
    fold, limit = job
    import ijson
    t0 = time.time()
    R = G["A_merged"]["exact_tie|d100"]
    cm = R["meta"]
    sel = np.flatnonzero(cm["fold"] == fold)
    if limit:
        sel = sel[:limit]
    want = set(cm["md5"][sel].tolist())
    ref_set = set(jload(os.path.join(D, "prfba", f"f{fold:02d}_d100",
                                     "_cache_ec_subset.json")))
    recs = {}
    path = os.path.join(D, "prfba", f"f{fold:02d}_d100", "hits", "asv_top20_alignment_hits.json")
    with open(path, "rb") as fh:
        for md5, rec in ijson.kvitems(fh, "", use_float=True):
            if md5 in want and isinstance(rec, dict):
                ranked, info = stage0_ranked(rec)
                recs[md5] = (ranked[:max(TOPK_K)], info)
    t1 = time.time()
    E = G["E"]
    gidx = G["gidx"]
    NU = G["NU"]
    off = np.concatenate([[0], np.cumsum(R["reps_len"])])
    K = TOPK_K
    feats = collections.defaultdict(list)
    A = collections.Counter()
    for i in sel:
        md5 = cm["md5"][i]
        if md5 not in recs:
            A["record_missing"] += 1
            continue
        ranked, info = recs[md5]
        reps = R["reps_flat"][off[i]:off[i + 1]]
        if any(g not in gidx for g, _ in ranked):
            A["hit_genome_outside_panel"] += 1
            ranked = [(g, x) for g, x in ranked if g in gidx]
        if not ranked:
            A["no_stage0_survivors"] += 1
            continue
        gl = [gidx[g] for g, _ in ranked]
        ids = np.array([x for _, x in ranked], float)
        # ---- leakage (top-k neighbours = every genome any k <= max(TOPK_K) uses)
        si_ = int(cm["src"][i])
        ga = np.asarray(gl, np.int64)
        A["leak_source_in_topk_neighbours"] += int(si_ in gl)
        for nm in ("genus", "species", "taxon", "group"):
            arr = G[f"{nm}_code"]
            A[f"leak_topk_neighbour_shares_{nm}_with_source"] += int((arr[ga] == arr[si_]).any())
        A["topk_neighbours_outside_fold_reference"] += sum(1 for g, _ in ranked
                                                           if g not in ref_set)
        tie = {gl[j] for j in range(len(gl)) if ids[j] == ids[0]}
        if tie != set(reps.tolist()):
            A["tie_set_ne_exact_tie_reps"] += 1
        A["cells"] += 1
        T = E[cm["src"][i]].astype(bool)
        ts = int(T.sum())
        w = np.exp(-(ids[0] - ids) / TOPK_TAU)
        Mx = E[gl].astype(np.float64) * w[:, None]
        cum = np.cumsum(Mx, axis=0)
        wc = np.cumsum(w)
        feats["cell"].append(i)
        feats["saturated"].append(info["saturated"])
        feats["n_survivors"].append(info["n_survivors"])
        # nearest-reference comparators
        r1 = cm["rank1_gi"][i]
        for nm, g in (("rank1_score_copy", r1), ("rank1_identity_copy", gl[0])):
            if g is None or g < 0:
                feats[f"{nm}.F1"].append(np.nan)
                feats[f"{nm}.TP"].append(0)
                feats[f"{nm}.PRED"].append(0)
                A[f"{nm}_missing"] += 1
                continue
            S = E[g].astype(bool)
            tp = int((S & T).sum())
            ps = int(S.sum())
            feats[f"{nm}.F1"].append(2 * tp / (ps + ts) if ps + ts else 0.0)
            feats[f"{nm}.TP"].append(tp)
            feats[f"{nm}.PRED"].append(ps)
        for k in K:
            kk = min(k, len(gl))
            q = cum[kk - 1] / wc[kk - 1]
            q = np.where(q >= 1 - 1e-9, 1.0, q)
            supp = q > 1e-12
            qs = q[supp]
            ys = T[supp]
            N = int(supp.sum())
            POS = int(ys.sum())
            FN = ts - POS
            sq = float(((qs - ys) ** 2).sum())
            pre = f"k{k}."
            feats[pre + "N"].append(N)
            feats[pre + "POS"].append(POS)
            feats[pre + "SQ"].append(sq)
            feats[pre + "FN"].append(FN)
            for e in (0.01, 0.001):
                qc = np.clip(qs, e, 1 - e)
                ll = float(-(ys * np.log(qc) + (~ys) * np.log(1 - qc)).sum())
                feats[pre + f"LLu{e}"].append(ll)
                feats[pre + f"LLv{e}"].append(ll - FN * math.log(e) -
                                              (NU - N - FN) * math.log(1 - e))
            db = np.where(qs >= 1.0, 10, np.clip(np.ceil(qs * 10 - 1e-12).astype(int) - 1, 0, 9))
            feats[pre + "DB_M"].append(np.bincount(db, minlength=11))
            feats[pre + "DB_POS"].append(np.bincount(db[ys], minlength=11))
            feats[pre + "DB_SP"].append(np.bincount(db, weights=qs, minlength=11))
            for rule, thr in (("q_ge_0.5", 0.5), ("q_gt_0", 1e-12)):
                pm = q >= thr if rule == "q_ge_0.5" else supp
                tp = int((pm & T).sum())
                ps = int(pm.sum())
                feats[pre + f"{rule}.TP"].append(tp)
                feats[pre + f"{rule}.PRED"].append(ps)
                feats[pre + f"{rule}.F1"].append(2 * tp / (ps + ts) if ps + ts else 0.0)
            # HEURISTIC, not a bound: the 200-hit cutoff is by align_score, not
            # identity, so a hit beyond rank 200 can out-rank the k-th survivor
            feats[pre + "exposed"].append(bool(info["saturated"] and (
                info["n_survivors"] < k or (info["last_qualified_identity"] is not None and
                                            ids[kk - 1] <= info["last_qualified_identity"]))))
            feats[pre + "short"].append(bool(info["saturated"] and info["n_survivors"] < k))
    out = {k: np.asarray(v) for k, v in feats.items()}
    return {"fold": fold, "feats": out, "asserts": dict(A),
            "seconds_stream": round(t1 - t0, 1), "seconds": round(time.time() - t0, 1),
            "peak_rss_gb": round(rss_gb(), 2)}


def topk_analyze(parts, B, seed, eps_list):
    R = G["A_merged"]["exact_tie|d100"]
    cm = R["meta"]
    F = {k: np.concatenate([p["feats"][k] for p in parts]) for k in parts[0]["feats"]}
    cell = F["cell"].astype(np.int64)
    NU = G["NU"]
    # exact-tie p and binary rules on the same cells, from the primary EC arrays
    n_cell = cm["n"].astype(np.int64)
    pc = precompute_space(R["spaces"]["ec"], n_cell, NU, eps_list)
    src = cm["src"][cell]
    n = n_cell[cell]
    pops = [("all", np.ones(len(cell), bool)), ("exact_tie_single", n == 1),
            ("exact_tie_multi", n >= 2)]
    sens = [("all_records", np.ones(len(cell), bool)), ("non_saturated_records",
                                                         ~F["saturated"].astype(bool))]
    out = {}
    rel = []
    for pname, pm in pops:
        for sname, sm in sens:
            mk = pm & sm
            if not mk.any():
                continue
            idx = np.flatnonzero(mk)
            srcs, sl = np.unique(src[idx], return_inverse=True)
            Gs = len(srcs)
            schemes = make_schemes(srcs, B, seed, {"group": G["group_code"], "src": None,
                                                   "genus": G["genus_code"]})
            names = list(schemes)
            X = Cols(Gs)

            def ps_(v):
                return np.bincount(sl, weights=np.asarray(v, float)[idx], minlength=Gs)
            X.add("NCELL", np.bincount(sl, minlength=Gs))
            X.add("T", ps_(pc["T"][cell]))
            X.add("et.N", ps_(pc["N"][cell]))
            X.add("et.POS", ps_(pc["POS"][cell]))
            X.add("et.SQ", ps_(pc["SQ"][cell]))
            X.add("et.FN", ps_(pc["FN"][cell]))
            for e in eps_list:
                X.add(f"et.LL{e}", ps_(pc[f"LL{e}"][cell]))
            for r in ("union", "p_ge_0.5"):
                for f in ("TP", "PRED", "F1"):
                    X.add(f"et.{r}.{f}", ps_(pc[f"{r}.{f}"][cell]))
            for nm in ("rank1_score_copy", "rank1_identity_copy"):
                X.add(f"{nm}.F1", ps_(np.nan_to_num(F[f"{nm}.F1"])))
                X.add(f"{nm}.TP", ps_(F[f"{nm}.TP"]))
                X.add(f"{nm}.PRED", ps_(F[f"{nm}.PRED"]))
            for k in TOPK_K:
                pre = f"k{k}."
                for f in ("N", "POS", "SQ", "FN", "LLv0.01", "LLv0.001", "LLu0.01",
                          "q_ge_0.5.TP", "q_ge_0.5.PRED", "q_ge_0.5.F1", "q_gt_0.TP",
                          "q_gt_0.PRED", "q_gt_0.F1", "exposed", "short"):
                    X.add(pre + f, ps_(F[pre + f]))
                for f in ("DB_M", "DB_POS", "DB_SP"):
                    arr = F[pre + f][idx]
                    X.add(pre + f + "[]", np.vstack([np.bincount(sl, weights=arr[:, j],
                                                                 minlength=Gs)
                                                     for j in range(11)]).T)
            XM = X.matrix()
            S = XM.sum(0)
            SB = collections.OrderedDict(
                (s_, (W_ @ agg_rows(XM, inv_, K_)) if W_ is not None else None)
                for s_, (inv_, W_, K_) in schemes.items())
            L = X.L

            def derive(S):
                g = lambda k: S[..., L[k]]  # noqa: E731
                o = {}
                NC, T = g("NCELL"), g("T")
                Nv = NC * NU
                with np.errstate(all="ignore"):
                    et = {"brier_union_support": g("et.SQ") / g("et.N"),
                          "brier_universe": (g("et.SQ") + g("et.FN")) / Nv,
                          "logloss_universe_eps0.01": (g("et.LL0.01") - g("et.FN") * math.log(0.01)
                                                       - (Nv - g("et.N") - g("et.FN")) *
                                                       math.log(0.99)) / Nv,
                          "support_size_mean": g("et.N") / NC}
                    for r in ("union", "p_ge_0.5"):
                        et[f"{r}.macro_F1"] = g(f"et.{r}.F1") / NC
                        et[f"{r}.micro_F1"] = 2 * g(f"et.{r}.TP") / (g(f"et.{r}.PRED") + T)
                    o["exact_tie_k_over_n"] = et
                    for nm in ("rank1_score_copy", "rank1_identity_copy"):
                        o[nm] = {"macro_F1": g(f"{nm}.F1") / NC,
                                 "micro_F1": 2 * g(f"{nm}.TP") / (g(f"{nm}.PRED") + T)}
                    for k in TOPK_K:
                        pre = f"k{k}."
                        N, POS, SQ, FN = g(pre + "N"), g(pre + "POS"), g(pre + "SQ"), g(pre + "FN")
                        DBM, DBP, DBS = g(pre + "DB_M[]"), g(pre + "DB_POS[]"), g(pre + "DB_SP[]")
                        d = {"support_size_mean": N / NC,
                             "brier_support": SQ / N,
                             "brier_universe": (SQ + FN) / Nv,
                             "logloss_universe_eps0.01": g(pre + "LLv0.01") / Nv,
                             "logloss_universe_eps0.001": g(pre + "LLv0.001") / Nv,
                             "logloss_support_eps0.01": g(pre + "LLu0.01") / N,
                             "ece_decile_support": np.abs(DBS - DBP).sum(-1) / N,
                             "ece_decile_universe": (np.abs(DBS - DBP).sum(-1) + FN) / Nv,
                             "frac_cells_truncation_exposed_heuristic": g(pre + "exposed") / NC,
                             "frac_cells_saturated_and_short": g(pre + "short") / NC}
                        d["d_brier_universe_vs_exact_tie"] = d["brier_universe"] - et["brier_universe"]
                        for rule in ("q_ge_0.5", "q_gt_0"):
                            mi = 2 * g(pre + f"{rule}.TP") / (g(pre + f"{rule}.PRED") + T)
                            ma = g(pre + f"{rule}.F1") / NC
                            d[f"{rule}.micro_F1"] = mi
                            d[f"{rule}.macro_F1"] = ma
                            d[f"{rule}.micro_P"] = g(pre + f"{rule}.TP") / g(pre + f"{rule}.PRED")
                            d[f"{rule}.micro_R"] = g(pre + f"{rule}.TP") / T
                        d["d_macro_F1_q_ge_0.5_vs_exact_tie_union"] = d["q_ge_0.5.macro_F1"] - et["union.macro_F1"]
                        d["d_micro_F1_q_ge_0.5_vs_exact_tie_union"] = d["q_ge_0.5.micro_F1"] - et["union.micro_F1"]
                        d["d_macro_F1_q_ge_0.5_vs_rank1_score_copy"] = d["q_ge_0.5.macro_F1"] - o["rank1_score_copy"]["macro_F1"]
                        d["d_macro_F1_q_ge_0.5_vs_rank1_identity_copy"] = d["q_ge_0.5.macro_F1"] - o["rank1_identity_copy"]["macro_F1"]
                        d["d_macro_F1_q_ge_0.5_vs_exact_tie_p_ge_0.5"] = d["q_ge_0.5.macro_F1"] - et["p_ge_0.5.macro_F1"]
                        o[f"k={k}"] = d
                    o["d_macro_F1_exact_tie_union_vs_rank1_score_copy"] = et["union.macro_F1"] - o["rank1_score_copy"]["macro_F1"]
                return o

            pt = derive(S)
            bts = [derive(SB[s_]) if SB[s_] is not None else None for s_ in names]
            res = pack_schemes(pt, bts, names)
            res["n_cells"] = int(len(idx))
            res["n_src"] = int(Gs)
            res["ci_unit"] = names[0]
            res["clusters"] = {s_: int(K_) for s_, (_i, _W, K_) in schemes.items()}
            out[f"{pname}|{sname}"] = res
            for k in TOPK_K:
                pre = f"k{k}."
                Mb, Pb, Sb = S[L[pre + "DB_M[]"]], S[L[pre + "DB_POS[]"]], S[L[pre + "DB_SP[]"]]
                obs = {}
                for s_ in names:
                    if SB[s_] is not None:
                        with np.errstate(all="ignore"):
                            obs[s_] = SB[s_][:, L[pre + "DB_POS[]"]] / SB[s_][:, L[pre + "DB_M[]"]]
                for j in range(11):
                    if Mb[j] <= 0:
                        continue
                    c_ = ci_cols(Pb[j] / Mb[j], {s_: o[:, j] for s_, o in obs.items()}, names)
                    rel.append({"analysis": "exploratory_topk", "run": f"topk|k={k}",
                                "mode": "exact_tie_cells", "density": "1.0", "forced": False,
                                "space": "ec", "pop": pname, "stratum_kind": "records",
                                "stratum_value": sname, "denominator": "union",
                                "bin_scheme": "decile", "bin": DEC_LABELS[j],
                                "p_lo": DEC_LO[j], "p_hi": DEC_HI[j], "n_points": int(Mb[j]),
                                "n_present": int(Pb[j]), "n_cells": None, "n_src": None,
                                "n_distinct_pairs": None, "mean_p": Sb[j] / Mb[j],
                                "obs_freq": Pb[j] / Mb[j], **c_})
    return out, rel


# ============================================================================ README
def fmt_ci(d, k=3):
    """'est [lo, hi]' with the PRIMARY-unit CI (holdout group; genus for the secondary)."""
    if not d or d.get("est") is None:
        return "n/a"
    if d.get("lo") is None:
        return f"{d['est']:.{k}f}"
    return f"{d['est']:.{k}f} [{d['lo']:.{k}f}, {d['hi']:.{k}f}]"


def fmt_sens(d, k=3):
    """The sensitivity CIs of one estimate: '(genus [..]; source [..])'."""
    if not d:
        return ""
    parts = []
    for key, nm in (("ci_genus", "genus"), ("ci_src", "source")):
        v = d.get(key)
        if v and v[0] is not None:
            parts.append(f"{nm} [{v[0]:.{k}f}, {v[1]:.{k}f}]")
    return f" ({'; '.join(parts)})" if parts else ""


def fmt_all(d, k=3):
    return fmt_ci(d, k) + fmt_sens(d, k)


def _e(d):
    return None if not d else d.get("est")


def sgn(d):
    """+1 / -1 when the primary CI excludes 0 on that side, 0 when it covers 0,
    None when there is no CI."""
    if not d or d.get("lo") is None:
        return None
    return 1 if d["lo"] > 0 else (-1 if d["hi"] < 0 else 0)


def sens_flip(d):
    """Names of the sensitivity schemes whose CI disagrees with the primary on 'covers 0'."""
    s0 = sgn(d)
    out = []
    for key, nm in (("ci_genus", "genus"), ("ci_src", "source")):
        v = (d or {}).get(key)
        if s0 is None or not v or v[0] is None:
            continue
        s = 1 if v[0] > 0 else (-1 if v[1] < 0 else 0)
        if (s == 0) != (s0 == 0):
            out.append(nm)
    return out


def word0(d, pos="> 0", neg="< 0", zero="CI covers 0: not distinguishable from 0"):
    s = sgn(d)
    return {1: pos, -1: neg, 0: zero, None: "no CI"}[s]


def sesoi_verdict(d, sesoi):
    """score_folds' rule: exceeded iff |delta| >= SESOI and the CI excludes 0."""
    if not d or d.get("est") is None or d.get("lo") is None:
        return "indeterminate (no CI)"
    if abs(d["est"]) >= sesoi and (d["lo"] > 0 or d["hi"] < 0):
        lb = d["lo"] if d["est"] > 0 else -d["hi"]
        return ("exceeded" + ("" if lb >= sesoi else
                              f" (but the CI bound nearer 0, {lb:+.3f}, is inside the SESOI)"))
    if d["lo"] > 0 or d["hi"] < 0:
        return "not exceeded (non-zero, below the SESOI)"
    return "not exceeded (CI covers 0)"


def _rel(rel_df, run, space, pop, scheme, binlabel):
    if rel_df is None:
        return None
    x = rel_df[(rel_df.run == run) & (rel_df.space == space) & (rel_df["pop"] == pop) &
               (rel_df.stratum_kind == "all") & (rel_df.bin_scheme == scheme) &
               (rel_df.bin == binlabel)]
    return None if not len(x) else x.iloc[0]


def _relfmt(r, k=3):
    if r is None:
        return "n/a"
    ci = (f" [{r.ci_lo:.{k}f}, {r.ci_hi:.{k}f}]" if r.ci_lo == r.ci_lo and
          r.ci_lo is not None else "")
    return f"{r.obs_freq:.{k}f}{ci}"


def _maxrule_note(P):
    try:
        d1 = P["exact_tie|d100"]["rxn_ec"]["_diag"]
        d2 = P["exact_tie|d100"]["rxn_role"]["_diag"]
        return (" The max-over-genes reaction weight differs from genome-level prevalence on "
                f"{100 * d1['max_rule_ne_genome_prevalence'] / d1['points']:.2f}% (rxn_ec) and "
                f"{100 * d2['max_rule_ne_genome_prevalence'] / d2['points']:.2f}% (rxn_role) of "
                "multi-cell union reactions.")
    except (KeyError, TypeError, ZeroDivisionError):
        return ""


def write_readme(path, scores, assertions, args, rel_df):
    P = scores["primary_controlled"]
    sesoi = scores["meta"]["bootstrap"].get("sesoi_delta_f1") or 0.02
    L = []
    w = L.append

    def g(run, space, stratum):
        return P.get(run, {}).get(space, {}).get(stratum)

    def et(space, stratum="multi|all=all"):
        return g("exact_tie|d100", space, stratum)

    R = et("ec")
    if not R:
        return
    U = R["union"]
    allc = et("ec", "all|all=all")
    n_all = allc["n_cells"] if allc else None
    cl = R.get("clusters", {})
    w("# Calibration of prFBA's p = k/n weighting\n")
    _sha = hashlib.sha256(open(os.path.abspath(__file__), "rb").read()).hexdigest()[:12]
    w(f"`scripts/calibrate_pfunc.py` (sha256 {_sha}), {time.strftime('%Y-%m-%d')} (commit "
      f"{scores['meta'].get('git_commit', '?')[:10]}); read-only on existing benchmark "
      "artefacts. Question: prFBA gives each function x of a selected reference set of n "
      "genomes the weight p(x) = k_x/n. Is p a calibrated probability that x is present in "
      "the organism the amplicon came from?\n")
    w("## Key results\n")
    w(f"Population: controlled genus-blocked benchmark, exact_tie d100, EC, multi-genome "
      f"cells (n >= 2; {R['n_cells']:,} cells, {R['n_src']:,} source genomes, "
      f"{cl.get('group', {}).get('n', '?')} holdout groups), union denominator, unless "
      f"stated. **95% CIs are a cluster bootstrap over holdout groups** (score_folds' "
      f"pre-specified unit, {args.bootstrap} draws, seed {args.seed}); genus and "
      "source-genome intervals are sensitivities. The largest group holds "
      f"{100 * (cl.get('group', {}).get('largest_share_points') or 0):.0f}% of the "
      "multi-cell points, so the group interval mostly reflects how often that one cluster "
      "is drawn; the genus interval is the middle ground.\n")

    # 1. calibration: reliability term + probability-scale slope first
    sl = U["calib_slope_wls"]
    over = sl.get("hi") is not None and sl["hi"] < 1
    r12 = _rel(rel_df, "exact_tie|d100", "ec", "n2", "exact_or_decile", "n=2,p=1/2")
    r14 = _rel(rel_df, "exact_tie|d100", "ec", "n4+", "exact_or_decile", "n=4,p=1/4")
    p1m = _rel(rel_df, "exact_tie|d100", "ec", "multi", "distinct_p", "p=1/1")
    w(f"1. **{'p is overconfident' if over else 'Calibration'}.** The Murphy reliability "
      f"term (the part of p's Brier score that perfect recalibration would remove) is "
      f"{fmt_all(U['murphy_reliability'], 4)} on a Brier of {fmt_ci(U['brier_p'], 4)}; "
      f"out-of-fold isotonic recalibration recovers "
      f"{fmt_ci(U.get('oof_gap_p_minus_isotonic'), 4)} of it. The probability-scale "
      f"calibration slope (weighted least squares of observed presence on p over the "
      f"distinct-p bins, weights = points; 1 = calibrated, < 1 = too extreme) is "
      f"{fmt_all(sl, 2)}, intercept {fmt_ci(U['calib_intercept_wls'], 2)}. "
      + ("The overconfidence is at p >= 1/2" if (r12 is not None and p1m is not None and
                                                r12.obs_freq < 0.5 and p1m.obs_freq < 1)
         else "At p >= 1/2") +
      f" (n = 2, p = 1/2 is present {_relfmt(r12)}; p = 1 from n >= 2 references "
      f"{_relfmt(p1m)}); "
      + ("low p is slightly under-confident" if (r14 is not None and r14.obs_freq > 0.25)
         else "at low p") +
      f" (n = 4, p = 1/4: {_relfmt(r14)}). The logistic slope on logit p is set by where "
      f"p = 1 is clipped and is secondary: {fmt_ci(U['calib_slope_eps0.01'], 2)} at eps "
      f"0.01, {fmt_ci(U['calib_slope_eps0.001'], 2)} at eps 0.001.")

    # 2. ranking
    if rel_df is not None:
        ex = rel_df[(rel_df.run == "exact_tie|d100") & (rel_df.space == "ec") &
                    (rel_df.stratum_kind == "all") & (rel_df.bin_scheme == "exact_or_decile") &
                    (rel_df["pop"].isin(["n2", "n3", "n4+"])) & rel_df.bin.str.startswith("n=")]
        mono = all(np.all(np.diff(grp.sort_values("p_lo").obs_freq.values) > 0)
                   for _, grp in ex.groupby("pop"))
        sep = all(np.all(grp.sort_values("p_lo").ci_hi.values[:-1] <
                         grp.sort_values("p_lo").ci_lo.values[1:])
                  for _, grp in ex.groupby("pop"))
        n3 = ex[ex["pop"] == "n3"].sort_values("p_lo")
        cells = "; ".join(f"{b.split(',')[1]}: {_relfmt(r)}" for b, r in
                          zip(n3.bin, n3.itertuples(index=False)))
        w(f"2. **{'p ranks functions correctly' if mono else 'p does not rank monotonically'}.** "
          f"Within n = 2, 3 and 4 observed presence rises strictly with p"
          + (", and adjacent bins' group CIs do not overlap" if sep else
             " (point estimates; some adjacent CIs overlap)") + f"; n = 3: {cells}.")

    # 3. skill vs constant, by n bucket
    nb = {pop: et("ec", f"{pop}|all=all") for pop in ("n2", "n3", "n4+")}
    s2 = nb["n2"]["union"]["bss_vs_const"] if nb["n2"] else None
    sh2c = nb["n2"]["n_cells"] / R["n_cells"] if nb["n2"] else float("nan")
    sh2p = nb["n2"]["n_points_union"] / R["n_points_union"] if nb["n2"] else float("nan")
    if sgn(s2) == 0:
        n2txt = ("**for two-way ties (n = 2) p does not beat a constant confidence**: "
                 f"skill {fmt_all(s2)}, CI covers 0")
    else:
        n2txt = f"for two-way ties (n = 2) skill is {fmt_all(s2)} ({word0(s2)})"
    s3, s4 = (nb[k]["union"]["bss_vs_const"] if nb[k] else None for k in ("n3", "n4+"))
    head3 = ("Skill over a constant confidence comes from n >= 3 only"
             if (sgn(s2) == 0 and sgn(s3) == 1 and sgn(s4) == 1) else
             "Skill over a constant confidence, by n bucket")
    w(f"3. **{head3}.** Pooled multi "
      f"Brier skill vs the constant c = pooled precision is {fmt_all(U['bss_vs_const'])} "
      f"({word0(U['bss_vs_const'])}). By bucket, {n2txt}. n = 2 is "
      f"{100 * sh2c:.0f}% of multi cells and {100 * sh2p:.0f}% of their points. "
      f"n = 3: {fmt_ci(nb['n3']['union']['bss_vs_const'])}; n >= 4: "
      f"{fmt_ci(nb['n4+']['union']['bss_vs_const'])}.")

    # 4. p = 1 multi vs single, p = 0
    s1 = et("ec", "single|all=all")
    dm = allc["union"].get("obs_p1_multi_minus_single") if allc else None
    head4 = {1: "p = 1 is more often right when several references agree", 0: "Agreement among references "
             "does not resolvably help", -1: "Agreement among references hurts",
             None: "p = 1 multi vs single"}[sgn(dm)]
    w(f"4. **{head4}.** p = 1 from n >= 2 references is present "
      f"{_relfmt(p1m)} of the time, p = 1 from a single (cross-genus) reference "
      f"{fmt_ci(s1['union']['c_pooled_precision'])}; paired difference "
      f"{fmt_all(dm)} ({word0(dm)}). Functions outside the union (p = 0) are present "
      f"{fmt_ci(R['universe']['obs_freq_p0'])}.")

    # 5. binary rules vs SESOI (macro F1 = the benchmark's metric)
    b = {pop: (nb[pop] or {}).get("binary") for pop in nb}
    ba = allc["binary"] if allc else None
    if all(b.values()) and ba:
        def sh(pop):
            return f"{100 * nb[pop]['n_cells'] / n_all:.1f}%"
        m3 = b["n3"]["p_ge_0.5"]["d_macro_F1_vs_union"]
        m4 = b["n4+"]["p_ge_0.5"]["d_macro_F1_vs_union"]
        mall = ba["p_ge_0.5"]["d_macro_F1_vs_union"]
        inter = [sgn(b[k]["p_eq_1"]["d_macro_F1_vs_union"]) for k in ("n2", "n3", "n4+")]
        head5 = ("Thresholding p barely moves the benchmark"
                 if sesoi_verdict(mall, sesoi).startswith("not exceeded") else
                 "Thresholding p moves the benchmark")
        w(f"5. **{head5}.** Against the pre-registered "
          f"SESOI of {sesoi} F1 (fold_plan meta.frozen.sesoi_delta_f1; score_folds' rule: "
          f"exceeded iff |dF1| >= SESOI and the CI excludes 0), a majority vote (p >= 0.5) "
          f"vs the union changes macro F1 by {fmt_ci(m3)} at n = 3 ({sh('n3')} of the "
          f"{n_all:,} scored cells; micro {fmt_ci(b['n3']['p_ge_0.5']['d_micro_F1_vs_union'])}): "
          f"{sesoi_verdict(m3, sesoi)}; and by {fmt_all(m4)} at n >= 4 ({sh('n4+')} of "
          f"cells; micro {fmt_ci(b['n4+']['p_ge_0.5']['d_micro_F1_vs_union'])}): "
          f"{sesoi_verdict(m4, sesoi)}. At n = 2 it is the union. Over all scored cells "
          f"the rule changes the benchmark's macro F1 by {fmt_ci(mall, 4)}: "
          f"{sesoi_verdict(mall, sesoi)}. The intersection (p = 1; at n = 2 also p > 0.5) "
          + ("loses to the union in every n bucket" if all(v == -1 for v in inter) else
             "does not lose in every bucket") +
          f" (macro dF1 at n = 2 {fmt_ci(b['n2']['p_eq_1']['d_macro_F1_vs_union'])}, n >= 4 "
          f"{fmt_ci(b['n4+']['p_eq_1']['d_macro_F1_vs_union'])}).")

    # 6. units (paired contrasts)
    C = P.get("exact_tie|d100", {}).get("_unit_contrasts", {}).get("multi|all=all", {})

    def con(a, bb, k="murphy_reliability"):
        x = C.get(f"{a}-{bb}") or {}
        return x.get(k)
    trio = [("ec", "rxn_ec"), ("ec", "rxn_role"), ("rxn_ec", "rxn_role")]
    tri = [(a, bb, con(a, bb)) for a, bb in trio if con(a, bb)]
    indist = all(sgn(d) == 0 for _, _, d in tri)
    uparts = []
    for sp in args.spaces:
        x = et(sp)
        if x:
            uparts.append(f"{sp} {fmt_ci(x['union']['murphy_reliability'], 3)}")
    s_er, s_eo, s_ro = (sgn(con(a, bb)) for a, bb in trio)
    if indist:
        head6 = "EC and the two reaction layers are indistinguishable in calibration"
    elif s_er == 1 and s_eo == 1 and s_ro == 0:
        head6 = (f"EC is marginally worse calibrated than either reaction layer (by "
                 f"{_e(con('ec', 'rxn_ec')):.4f} and {_e(con('ec', 'rxn_role')):.4f} on a "
                 f"reliability term of {_e(U['murphy_reliability']):.3f}); the two reaction "
                 "layers are indistinguishable")
    else:
        head6 = "EC and the reaction layers differ in calibration"
    w(f"6. **{head6}** (reliability term {'; '.join(uparts)}). Paired differences from "
      "the shared draws: " + "; ".join(f"{a} - {bb} {fmt_ci(d, 4)}" for a, bb, d in tri)
      + (" (all cover 0)" if indist else "") + ". "
      + ("PGFam and the product-string proxy are worse calibrated"
         if (sgn(con("ec", "pgfam")) == -1 and sgn(con("ec", "product")) == -1) else
         "PGFam and product strings vs EC")
      + f" (ec - pgfam {fmt_ci(con('ec', 'pgfam'), 3)}, ec - product "
      f"{fmt_ci(con('ec', 'product'), 3)}). Skill vs c: " + "; ".join(
          f"{sp} {fmt_ci(et(sp)['union']['bss_vs_const'])}"
          f"{'' if sgn(et(sp)['union']['bss_vs_const']) else ' (covers 0)'}"
          for sp in args.spaces if et(sp)) + ".")

    # 7. novelty, forced
    nv = [(st, et("ec", f"multi|novelty_bin={st}")) for st in ("0.865-0.945", "0.945-0.987")]
    fo = g("exact_tie_forced|d100", "ec", "multi|novelty_bin=<0.865")
    if all(x for _, x in nv):
        txt = "; ".join(f"rank-1 identity {st} ({x['n_cells']:,} cells) "
                        f"{fmt_all(x['union']['bss_vs_const'])}, {word0(x['union']['bss_vs_const'])}"
                        for st, x in nv)
        sv = [sgn(x["union"]["bss_vs_const"]) for _, x in nv]
        head7 = ("Skill exists only when the nearest reference is close"
                 if sv == [0, 1] else "Skill by how close the nearest reference is")
        w(f"7. **{head7}**: {txt}." +
          (f" Forced answers below the family floor (exact_tie_forced, rank-1 identity "
           f"< 0.865, {fo['n_cells']:,} multi cells): {fmt_ci(fo['union']['bss_vs_const'])}"
           f"{', worse than a constant' if sgn(fo['union']['bss_vs_const']) == -1 else ''}."
           if fo else ""))

    # 8. other runs
    mods = []
    for rk in ("legacy|d100", "graduated|d100", "exact_tie_guards|d100", "exact_tie|d050",
               "exact_tie|d025", "exact_tie_forced|d100"):
        x = g(rk, "ec", "multi|all=all")
        if x:
            d = x["union"]["bss_vs_const"]
            fl = sens_flip(d)
            mods.append(f"{rk} {fmt_ci(d)}" + (" **covers 0**" if sgn(d) == 0 else "") +
                        (f" (excludes 0 under {', '.join(fl)})" if fl and sgn(d) == 0 else
                         (f" (covers 0 under {', '.join(fl)})" if fl else "")))
    if mods:
        w("8. **Other runs** (EC multi skill vs c): " + "; ".join(mods) + ".")

    # 9. secondary, 10. exploratory
    sec = scores.get("secondary_as_shipped", {})
    try:
        i_ = sec["include"]["multi|all=all"]["union"]
        e_ = sec["exclude_0.987"]["multi|all=all"]["union"]
        w(f"9. **As shipped (secondary; CIs over source genus)**: include-self multi skill "
          f"{fmt_ci(i_['bss_vs_const'])} is a ceiling (the source genome votes); "
          f"exclude_0.987 {fmt_ci(e_['bss_vs_const'])}, WLS slope "
          f"{fmt_ci(e_['calib_slope_wls'], 2)}, a post-hoc exclusion on truncated top-20 "
          f"records.")
    except (KeyError, TypeError):
        pass
    tk = scores.get("exploratory_topk", {})
    ta, tn = tk.get("all|all_records"), tk.get("all|non_saturated_records")
    if ta:
        d = ta.get("k=20", {})
        dn = (tn or {}).get("k=20", {})
        w(f"10. **Exploratory, not pre-registered**: top-20 identity-weighted prevalence "
          f"(tau 0.01) thresholded at 0.5 on all {ta['n_cells']:,} exact_tie cells changes "
          f"macro F1 vs the shipped union by "
          f"{fmt_all(d.get('d_macro_F1_q_ge_0.5_vs_exact_tie_union'), 4)} "
          f"(SESOI rule, were it pre-registered: "
          f"{sesoi_verdict(d.get('d_macro_F1_q_ge_0.5_vs_exact_tie_union'), sesoi)}); "
          f"universe Brier {fmt_ci(d.get('brier_universe'), 4)} vs "
          f"{fmt_ci(ta['exact_tie_k_over_n']['brier_universe'], 4)} for k/n. On non-saturated "
          f"records ({(tn or {}).get('n_cells', 0):,} cells: the 200th coverage-qualified hit is "
          f"already below the family floor; every record holds 200 hits, so this too assumes "
          f"identity tracks align_score beyond rank 200) the "
          f"macro F1 change is {fmt_ci(dn.get('d_macro_F1_q_ge_0.5_vs_exact_tie_union'), 4)}; "
          f"that is the only truncation evidence (the 'exposed' flag in scores.json is a "
          f"heuristic: the cutoff is by align_score, not identity). A hypothesis for a "
          f"pre-registered test, not a result.")
    w("")

    # ------------------------------------------------------------------ tables
    w("## Tables\n")
    w("exact_tie d100, multi-genome cells, union denominator. CIs over holdout groups. "
      "REL = Murphy reliability term; WLS = probability-scale slope; logit slope at eps "
      "0.01 / 0.001.\n")
    w("| unit | REL | WLS slope | skill vs c | OOF isotonic gap | logit slope .01 / .001 |")
    w("|---|---|---|---|---|---|")
    for s in args.spaces:
        r = et(s)
        if not r:
            continue
        u = r["union"]
        w(f"| {s} | {fmt_ci(u['murphy_reliability'], 4)} | {fmt_ci(u['calib_slope_wls'], 2)} | "
          f"{fmt_ci(u['bss_vs_const'], 3)} | {fmt_ci(u.get('oof_gap_p_minus_isotonic'), 4)} | "
          f"{_e(u['calib_slope_eps0.01']):.2f} / {_e(u['calib_slope_eps0.001']):.2f} |")
    w("\nBy n bucket (EC):\n")
    w("| bucket | cells (% of scored) | % of multi points | skill vs c | REL | WLS slope |")
    w("|---|---|---|---|---|---|")
    for pop in ("n2", "n3", "n4+", "multi"):
        r = et("ec", f"{pop}|all=all")
        if not r:
            continue
        u = r["union"]
        w(f"| {pop} | {r['n_cells']:,} ({100 * r['n_cells'] / n_all:.1f}%) | "
          f"{100 * r['n_points_union'] / R['n_points_union']:.0f}% | "
          f"{fmt_ci(u['bss_vs_const'])} | {fmt_ci(u['murphy_reliability'], 4)} | "
          f"{fmt_ci(u['calib_slope_wls'], 2)} |")
    w(f"\nBinary rules (EC), dF1 vs the union; SESOI {sesoi}:\n")
    w("| bucket | macro dF1 p>=0.5 | micro dF1 p>=0.5 | macro dF1 p=1 | verdict (p>=0.5, macro) |")
    w("|---|---|---|---|---|")
    for pop in ("n2", "n3", "n4+", "multi", "all"):
        r = et("ec", f"{pop}|all=all")
        if not r:
            continue
        bb = r["binary"]
        w(f"| {pop} | {fmt_ci(bb['p_ge_0.5']['d_macro_F1_vs_union'], 4)} | "
          f"{fmt_ci(bb['p_ge_0.5']['d_micro_F1_vs_union'], 4)} | "
          f"{fmt_ci(bb['p_eq_1']['d_macro_F1_vs_union'], 4)} | "
          f"{'identical to the union' if pop == 'n2' else sesoi_verdict(bb['p_ge_0.5']['d_macro_F1_vs_union'], sesoi)} |")
    if rel_df is not None:
        sub = rel_df[(rel_df.run == "exact_tie|d100") & (rel_df.space == "ec") &
                     (rel_df.stratum_kind == "all") &
                     (((rel_df["pop"].isin(["n2", "n3", "n4+", "single"])) &
                       (rel_df.bin_scheme == "exact_or_decile")) |
                      ((rel_df["pop"] == "multi") & (rel_df.bin_scheme == "p0")))]
        sub = sub[sub.bin.str.startswith("n=") | (sub.bin_scheme == "p0")]
        w("\nObserved presence at each exact p (EC; n >= 5 and the decile rows are in "
          "reliability.csv):\n")
        w("| pop | bin | points | cells | observed [group CI] | genus CI |")
        w("|---|---|---|---|---|---|")
        for r in sub.itertuples(index=False):
            gci = (f"[{r.ci_lo_genus:.3f}, {r.ci_hi_genus:.3f}]"
                   if r.ci_lo_genus == r.ci_lo_genus else "")
            w(f"| {r.pop} | {r.bin} | {int(r.n_points):,} | {int(r.n_cells):,} | "
              f"{_relfmt(r)} | {gci} |")
    w("\nEvery primary run (EC, multi, union):\n")
    w("| run | multi cells | skill vs c [group] | genus | source | REL |")
    w("|---|---|---|---|---|---|")
    for rk in P:
        r = g(rk, "ec", "multi|all=all")
        if not r:
            continue
        d = r["union"]["bss_vs_const"]
        gs = d.get("ci_genus") or [None, None]
        ss = d.get("ci_src") or [None, None]
        f2 = (lambda v: f"[{v[0]:.3f}, {v[1]:.3f}]" if v[0] is not None else "n/a")  # noqa
        w(f"| {rk} | {r['n_cells']:,} | {fmt_ci(d)} | {f2(gs)} | {f2(ss)} | "
          f"{fmt_ci(r['union']['murphy_reliability'], 4)} |")
    if ta:
        w("\nExploratory top-k (EC, all exact_tie cells; dF1 = macro F1 of q >= 0.5 minus the "
          "shipped union):\n")
        w("| k | universe Brier | d vs k/n | dF1, all records | dF1, non-saturated records |")
        w("|---|---|---|---|---|")
        for k in TOPK_K:
            d = ta.get(f"k={k}")
            dn = (tn or {}).get(f"k={k}", {})
            if not d:
                continue
            w(f"| {k} | {fmt_ci(d['brier_universe'], 4)} | "
              f"{fmt_ci(d['d_brier_universe_vs_exact_tie'], 4)} | "
              f"{fmt_ci(d['d_macro_F1_q_ge_0.5_vs_exact_tie_union'], 4)} | "
              f"{fmt_ci(dn.get('d_macro_F1_q_ge_0.5_vs_exact_tie_union'), 4)} |")
    w("")

    # ------------------------------------------------------------------ method
    pu = assertions.get("context", {}).get("product_unit", {})
    w("## Method\n")
    w("- **Population**: score_folds' scored cells (`cells.csv.gz`, space ec, outcome "
      "scored) of exact_tie, legacy, graduated, exact_tie_guards (d100), exact_tie d050 / "
      "d025 and exact_tie forced. Reps = `predictions.json pred[mode][md5].reps`; truth = "
      "the cell's `src_genome_id` (`_cache_panel_ec` through norm_ec); p is recomputed from "
      "the reps' function sets. Abstentions have no p. Strata: all, single, multi, n2, n3, "
      "n4+, crossed with novelty_bin, rank1_rank, region and a twin-free one-cell-per-md5 "
      "population.")
    w("- **Units**: `ec` (2,235 scored ECs); `rxn_ec` (score_folds' rxn space, w(R) = max "
      "over member ECs); `rxn_role` (PGFam role route, w(R) = max over member PGFams); "
      "`pgfam`; `product` (PGFam family_product strings minus empty and generic "
      f"'hypothetical protein' variants, {pu.get('distinct_generic_hypothetical_strings_dropped', '?')} "
      "distinct strings dropped; a proxy for the RAST strings the merge counts). Union "
      "denominator = the union's own functions; universe = the fixed vocabulary with p = 0 "
      "as a bin (ec, rxn_ec, rxn_role)." + _maxrule_note(P))
    w("- **Scores**: Brier of p, of the union (= 1 - precision) and of the constant c = "
      "pooled precision (in-sample; an out-of-fold c is also reported); skill = 1 - "
      "Brier(p)/Brier(c); Murphy decomposition over distinct p (exact; residual reported); "
      "WLS calibration line; logistic slope on logit(clip(p, eps)); ECE; out-of-fold "
      "isotonic / Platt recalibration (4 folds fit, 5th scored; refitting noise is not in "
      "the CI); binary rules union, p >= 0.5, p > 0.5, p = 1, micro and macro (= the "
      "benchmark's per-cell mean) F1.")
    w("- **Uncertainty**: score_folds.cluster_bootstrap's draw verbatim "
      "(`default_rng(seed).integers(0, K, size=(B, K))` over the stratum's sorted cluster "
      "keys), percentile 95%. Primary unit = holdout group (`group_id`); sensitivities = "
      "genus (`ci_genus`) and source genome (`ci_src`) in scores.json and reliability.csv. "
      "Sources nest in genera and genera in groups (asserted). One set of draws per "
      "stratum serves every unit, so unit contrasts (`_unit_contrasts`) are paired. The "
      "as-shipped secondary has no holdout groups; there the primary unit is the source "
      "genus.")
    w("- **Leakage** (asserted, every cell of every run): the source genome is never a rep "
      "or a top-k neighbour, and no rep or neighbour shares its genus, species, taxon_id or "
      "holdout group.\n")
    w("## Caveats\n")
    w("- Genus-blocked folds: every primary number is cross-genus transfer and pessimistic "
      "for deployment; single-reference p = 1 is not the manuscript's conspecific estimate. "
      "The panel is species-dereplicated, so multi sets are multi-species.")
    ratios = []
    for s in args.spaces:
        x = et(s)
        for k in ("bss_vs_const", "murphy_reliability", "calib_slope_wls", "brier_p",
                  "c_pooled_precision", "ece_distinct_p"):
            d = (x or {}).get("union", {}).get(k) or {}
            cs = d.get("ci_src") or [None, None]
            if d.get("lo") is not None and cs[0] is not None and cs[1] > cs[0]:
                ratios.append((d["hi"] - d["lo"]) / (cs[1] - cs[0]))
    rtxt = (f"; across the headline metrics of the five units they are "
            f"{min(ratios):.1f}-{max(ratios):.1f}x as wide as source-genome CIs"
            if ratios else "")
    w(f"- One holdout group (grp00000, all of fold 0) carries "
      f"{100 * (cl.get('group', {}).get('largest_share_points') or 0):.0f}% of the multi "
      f"points, so group CIs mostly reflect how often that group is drawn{rtxt}.")
    w("- Absence = absence of BV-BRC PGFam annotation on both sides; partial ECs dropped.")
    w("- The as-shipped secondary uses pseudo-folds for OOF recalibration and predates the "
      "20 -> 100 retained-hit rerun; include-self is a ceiling.")
    w("- Top-k is exploratory; its truncation evidence is the non-saturated-records "
      "population only, and that check is itself a heuristic, not a bound: it assumes "
      "identity follows align_score order beyond rank 200, which fails near the floor in "
      "a large share of those cells.\n")
    w("## Outputs and assertions\n")
    w("`scores.json` (per run x unit x stratum; `lo`/`hi` = primary unit, `ci_genus`, "
      "`ci_src`; `clusters` gives the cluster counts and the largest cluster's share), "
      "`reliability.csv` (long format, one stratification at a time), `assertions.json` "
      f"(every join, invariant and leakage check; all required checks zero: "
      f"**{assertions.get('all_required_zero')}**).\n")
    with open(path, "w") as fh:
        fh.write("\n".join(L))


# ============================================================================ main
def git_commit():
    import subprocess
    try:
        return subprocess.check_output(["git", "-C", SR, "rev-parse", "HEAD"],
                                       text=True).strip()
    except Exception:
        return "unknown"


def tojson(o):
    if isinstance(o, dict):
        return {str(k): tojson(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [tojson(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        f = float(o)
        return float(f"{f:.7g}") if math.isfinite(f) else None
    if isinstance(o, np.bool_):
        return bool(o)
    if isinstance(o, np.ndarray):
        return tojson(o.tolist())
    return o


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--modes", default="exact_tie,legacy,graduated,exact_tie_guards")
    ap.add_argument("--densities", default="1.0,0.5,0.25")
    ap.add_argument("--extra-density-modes", default="exact_tie",
                    help="modes run at the non-1.0 densities (spec: exact_tie)")
    ap.add_argument("--forced", action="store_true")
    ap.add_argument("--forced-modes", default="exact_tie")
    ap.add_argument("--spaces", default=",".join(ALL_SPACES))
    ap.add_argument("--bootstrap", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1729)
    ap.add_argument("--eps", default="0.01,0.001")
    ap.add_argument("--outdir", default=os.path.join(D, "calibration"))
    ap.add_argument("--cache-dir", default=None,
                    help="scratch directory for intermediate pickles (never the outdir)")
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--limit-cells", type=int, default=0, help="smoke test: cells per fold")
    ap.add_argument("--no-secondary", action="store_true")
    ap.add_argument("--no-topk", action="store_true")
    ap.add_argument("--reuse-layer-a", action="store_true",
                    help="reuse <cache-dir>/layerA.pkl if present")
    ap.add_argument("--write-only", action="store_true",
                    help="regenerate the four outputs from <cache-dir>/final.pkl")
    args = ap.parse_args()
    args.modes = [m for m in args.modes.split(",") if m]
    args.densities = [d for d in args.densities.split(",") if d]
    args.extra_density_modes = [m for m in args.extra_density_modes.split(",") if m]
    args.forced_modes = [m for m in args.forced_modes.split(",") if m]
    args.spaces = [s for s in args.spaces.split(",") if s]
    eps_list = [float(e) for e in args.eps.split(",") if e]
    folds = [int(f) for f in args.folds.split(",") if f]
    for s in args.spaces:
        if s not in ALL_SPACES:
            raise SystemExit(f"unknown space {s}")
    if args.cache_dir and os.path.abspath(args.cache_dir) == os.path.abspath(args.outdir):
        raise SystemExit("--cache-dir must not be the outdir")
    T0 = time.time()
    if args.write_only:
        with open(os.path.join(args.cache_dir, "final.pkl"), "rb") as fh:
            rel_all, scores, assertions = pickle.load(fh)
        write_outputs(args, rel_all, scores, assertions, T0)
        return
    runs = build_runs(args)
    G["runs_by_key"] = {r["key"]: r for r in runs}
    log(f"runs: {[r['key'] for r in runs]}; spaces {args.spaces}; workers {args.workers}")
    ctxm = mp.get_context("fork")

    # ---------------------------------------------------------------- secondary (bg)
    sec_pool = None
    sec_async = None
    if not args.no_secondary:
        sec_pool = ctxm.Pool(1)
        sec_async = sec_pool.apply_async(secondary_prepare, (args,))

    load_context(args, runs)
    load_cells(runs)

    # ---------------------------------------------------------------- layer A
    cache_a = os.path.join(args.cache_dir, "layerA.pkl") if args.cache_dir else None
    if args.reuse_layer_a and cache_a and os.path.exists(cache_a):
        with open(cache_a, "rb") as fh:
            merged, dir_asserts, jobstats = pickle.load(fh)
        log(f"layer A reused from {cache_a}")
    else:
        jobs = []
        by_dir = collections.defaultdict(list)
        for r in runs:
            by_dir[r["dir"]].append(r["key"])
        for dtag, rks in by_dir.items():
            for f in folds:
                jobs.append((dtag, f, rks, args.spaces, args.limit_cells))
        # biggest first (d100 carries four modes)
        jobs.sort(key=lambda j: -len(j[2]))
        log(f"layer A: {len(jobs)} (directory, fold) jobs on {args.workers} workers")
        results = []
        with ctxm.Pool(args.workers, maxtasksperchild=1) as pool:
            for res in pool.imap_unordered(layerA_job, jobs):
                results.append(res)
                log(f"  layer A f{res['fold']:02d}_{res['dir']}: {res['seconds']}s, "
                    f"worker peak rss {res['peak_rss_gb']} GB")
        merged, dir_asserts = merge_layerA(results, runs, args.spaces)
        jobstats = [{"job": f"f{r['fold']:02d}_{r['dir']}", "seconds": r["seconds"],
                     "peak_rss_gb": r["peak_rss_gb"]} for r in results]
        del results
        if cache_a:
            os.makedirs(args.cache_dir, exist_ok=True)
            with open(cache_a, "wb") as fh:
                pickle.dump((merged, dir_asserts, jobstats), fh, protocol=5)
    G["A_merged"] = merged
    log(f"layer A done ({time.time()-T0:.0f}s)")

    # ---------------------------------------------------------------- secondary result
    if sec_async is not None:
        SEC = sec_async.get()
        sec_pool.close()
        sec_pool.join()
        G["SEC"] = SEC
        G["AS_NU"] = SEC["AS_NU"]
        log(f"secondary prepared ({SEC['asserts'].get('seconds')}s)")

    # ---------------------------------------------------------------- top-k streaming
    topk_parts = None
    if not args.no_topk and "exact_tie|d100" in merged:
        log("exploratory top-k: streaming hits files")
        with ctxm.Pool(min(args.workers, len(folds)), maxtasksperchild=1) as pool:
            topk_parts = list(pool.imap_unordered(topk_job,
                                                  [(f, args.limit_cells) for f in folds]))
        for p in sorted(topk_parts, key=lambda x: x["fold"]):
            log(f"  top-k fold {p['fold']}: stream {p['seconds_stream']}s total "
                f"{p['seconds']}s peak rss {p['peak_rss_gb']} GB; {p['asserts']}")

    # ---------------------------------------------------------------- layer B
    log("layer B: statistics")
    jobsB = [(r["key"], args.spaces, args.bootstrap, args.seed, eps_list)
             for r in runs if r["key"] in merged]
    secjobs = ([(reg, args.bootstrap, args.seed, eps_list) for reg in G["SEC"]["runs"]]
               if "SEC" in G else [])
    primary, rel_all, secondary, bstats = {}, [], {}, []
    with ctxm.Pool(args.workers, maxtasksperchild=1) as pool:
        asyncs = [pool.apply_async(layerB_run, (j,)) for j in jobsB]
        sasyncs = [pool.apply_async(secondary_layerB, (j,)) for j in secjobs]
        for a in asyncs:
            res = a.get()
            primary[res["run"]] = res["scores"]
            rel_all.extend(res["rel"])
            bstats.append({"job": res["run"], "seconds": res["seconds"],
                           "peak_rss_gb": res["peak_rss_gb"]})
            log(f"  layer B {res['run']}: {res['seconds']}s")
        for a in sasyncs:
            res = a.get()
            secondary[res["regime"]] = res["scores"]
            rel_all.extend(res["rel"])
            bstats.append({"job": f"as_shipped|{res['regime']}", "seconds": res["seconds"]})
            log(f"  secondary {res['regime']}: {res['seconds']}s")
    topk_scores = {}
    if topk_parts:
        topk_scores, trel = topk_analyze(sorted(topk_parts, key=lambda x: x["fold"]),
                                         args.bootstrap, args.seed, eps_list)
        rel_all.extend(trel)
        log("  top-k analysed")

    # ---------------------------------------------------------------- assertions
    ctxA = G["A_ctx"]
    req = collections.OrderedDict()
    req["a_union_of_reps_EC_equals_pred_ec"] = sum(merged[r]["asserts"].get("union_ec_ne_pred_ec", 0) for r in merged)
    req["b_recomputed_f1_equals_cells_f1_within_1e-3"] = sum(merged[r]["asserts"].get("f1_mismatch", 0) for r in merged)
    req["c_recomputed_n_truth_equals_cells_n_truth"] = sum(merged[r]["asserts"].get("n_truth_mismatch", 0) for r in merged)
    req["d_len_reps_equals_cells_and_pred_n_reps"] = sum(merged[r]["asserts"].get("n_reps_mismatch", 0) for r in merged)
    req["e_ec_subset_cache_equals_panel_ec_cap_universe"] = sum(v["ec_subset_cache"]["mismatches"] for v in dir_asserts.values())
    req["f_pgfam_cache_cross_fold_conflicts"] = ctxA["pgfam_cache_conflicts_total"]
    req["g_pgfam_ec_projection_reproduces_panel_ec_universe"] = ctxA["pgfam_ec_projection"]["universe_restricted_mismatches"]
    req["h_recomputed_rxn_f1_equals_cells_rxn_f1_within_1e-3"] = sum(merged[r]["asserts"].get("rxn_f1_mismatch", 0) for r in merged)
    req["i_recomputed_rxn_n_truth_equals_cells_rxn_n_truth"] = sum(merged[r]["asserts"].get("rxn_n_truth_mismatch", 0) for r in merged)
    req["j_md5_missing_or_reps_outside_panel"] = sum(merged[r]["asserts"].get("md5_missing_in_predictions", 0) + merged[r]["asserts"].get("reps_or_src_outside_panel_or_empty", 0) for r in merged)
    req["k_reps_outside_fold_reference"] = sum(merged[r]["asserts"].get("reps_outside_fold_reference", 0) for r in merged)
    req["l_twin_md5_recount_equals_meta"] = sum(int(v["twin_md5"]["recomputed"] != v["twin_md5"]["meta_n_twin_md5"]) for v in dir_asserts.values())
    # ---- leakage, explicit (every primary run, every scored cell)
    def runsum(k):
        return sum(merged[r]["asserts"].get(k, 0) for r in merged)
    req["m_source_never_in_reps"] = runsum("leak_source_in_reps")
    req["n_no_rep_shares_genus_with_source"] = runsum("leak_rep_shares_genus_with_source")
    req["o_no_rep_shares_species_with_source"] = runsum("leak_rep_shares_species_with_source")
    req["p_no_rep_shares_taxon_id_with_source"] = runsum("leak_rep_shares_taxon_with_source")
    req["q_no_rep_shares_holdout_group_with_source"] = runsum("leak_rep_shares_group_with_source")
    # ---- cluster labels (bootstrap units)
    req["r_cells_group_id_equals_fold_plan_genome_group"] = runsum("group_id_ne_fold_plan_genome_group")
    req["u_panel_genomes_without_fold_plan_group"] = ctxA["panel_genomes_without_fold_plan_group"]
    req["v_genera_spanning_more_than_one_group"] = ctxA["panel_taxonomy"]["genera_spanning_more_than_one_group"]
    if "SEC" in G:
        sa = G["SEC"]["asserts"]
        req["s1_shared_ec_size_equals_2132"] = int(sa["shared_ec_size"] != 2132)
        req["s2_as_shipped_genomes_missing_function_sets"] = sa["genomes_missing_from_function_sets"]
        for reg, R in G["SEC"]["runs"].items():
            a = R["asserts"]
            req[f"s3_{reg}_f1_vs_comparison_cells"] = a.get("f1_mismatch_vs_comparison_cells", 0)
            req[f"s4_{reg}_n_truth_vs_comparison_cells"] = a.get("n_truth_mismatch_vs_comparison_cells", 0)
            req[f"s5_{reg}_scored_set_vs_comparison_ok"] = a.get("scored_here_but_not_cmp_ok", 0) + a.get("cmp_ok_but_not_scored_here", 0)
    if topk_parts:
        req["t1_topk_identity_tie_set_equals_exact_tie_reps"] = sum(p["asserts"].get("tie_set_ne_exact_tie_reps", 0) for p in topk_parts)
        req["t2_topk_hits_record_missing"] = sum(p["asserts"].get("record_missing", 0) for p in topk_parts)

        def tsum(k):
            return sum(p["asserts"].get(k, 0) for p in topk_parts)
        req["t4_source_never_in_topk_neighbours"] = tsum("leak_source_in_topk_neighbours")
        req["t5a_no_topk_neighbour_shares_genus_with_source"] = tsum("leak_topk_neighbour_shares_genus_with_source")
        req["t5b_no_topk_neighbour_shares_species_with_source"] = tsum("leak_topk_neighbour_shares_species_with_source")
        req["t5c_no_topk_neighbour_shares_taxon_id_with_source"] = tsum("leak_topk_neighbour_shares_taxon_with_source")
        req["t5d_no_topk_neighbour_shares_holdout_group_with_source"] = tsum("leak_topk_neighbour_shares_group_with_source")
        req["t6_topk_neighbours_all_in_fold_reference"] = tsum("topk_neighbours_outside_fold_reference")
        # rank-1 (score) copy F1 vs cells.csv base_nn1
        R = merged["exact_tie|d100"]
        F = {k: np.concatenate([p["feats"][k] for p in topk_parts]) for k in ("cell", "rank1_score_copy.F1")}
        nmis = nchk = nmiss = 0
        for ci_, f1 in zip(F["cell"], F["rank1_score_copy.F1"]):
            amp = R["meta"]["amplicon"][ci_]
            ref = G["nn1_f1"].get(amp)
            if ref is None:
                nmiss += 1
                continue
            nchk += 1
            nmis += int(abs(f1 - ref) > 1e-3)
        req["t3_rank1_score_copy_f1_equals_cells_base_nn1"] = nmis
        topk_nn1 = {"checked": nchk, "base_nn1_row_missing": nmiss}
    assertions = {
        "all_required_zero": all(v == 0 for v in req.values()),
        "required": req,
        "per_run": {r: merged[r]["asserts"] for r in merged},
        "per_directory": dir_asserts,
        "context": ctxA,
        "secondary": G.get("SEC", {}).get("asserts"),
        "secondary_per_regime": {k: v["asserts"] for k, v in G.get("SEC", {}).get("runs", {}).items()},
        "topk": ({"per_fold": {p["fold"]: p["asserts"] for p in topk_parts},
                  "base_nn1_crosscheck": topk_nn1} if topk_parts else None),
        "notes": {
            "a": "union of the reps' universe-restricted EC sets equals predictions.json ec, per cell",
            "b,c": "recomputed union F1 / n_truth equal cells.csv.gz (space ec) for every scored cell",
            "d": "len(reps) == pred n_reps == cells.csv.gz n_reps",
            "e": "_cache_ec_subset[g] == _cache_panel_ec[g].ec cap universe, every fold directory",
            "f": "PGFam caches agree across every fold directory (d100, d100_forced, d050, d025)",
            "g": "projecting the PGFam caches through function_maps.pgfam_ec reproduces _cache_panel_ec cap universe, all panel genomes (raw-set mismatches reported in context)",
            "h,i": "the rxn_ec unit reproduces score_folds' rxn space (cells.csv.gz space rxn)",
            "l": "twin md5s recounted from amplicon_md5_fanout.json equal predictions meta n_twin_md5",
            "m-q": "LEAKAGE, counted in cells, every primary run: the source genome is never one of the reps, and no rep shares its genus, species, taxon_id (panel_eligible) or holdout group (fold_plan genome_group)",
            "r": "cells.csv.gz group_id == fold_plan genome_group[src_genome_id] for every scored cell (the primary bootstrap unit)",
            "u,v": "every panel genome has a holdout group; no genus spans two groups (source genome nests in genus nests in group, so the three bootstrap schemes reweight the same per-source sums)",
            "t4-t6": "LEAKAGE, top-k: across every neighbour any k <= 100 uses (stage-0 survivors, deduplicated by genome), the source genome never appears, no neighbour shares its genus, species, taxon_id or holdout group, and every neighbour is in the fold's reference (_cache_ec_subset of f0X_d100)",
            "norm_ec": "context.panel_ec_entries_changed_by_norm_ec counts truth EC entries norm_ec rewrote or dropped (truth now goes through norm_ec, as score_folds.load_functions)",
            "s*": "as-shipped secondary reproduces compare.py (comparison_cells.csv.gz, method prfba_include / prfba_exclude, space ec)",
            "t1": "stage-0 + identity ranking + genome dedup reproduces the exact_tie reps as the identity-1 tie set (validates the top-k ranking)",
            "t3": "the rank-1 (align_score) copy comparator equals score_folds' base_nn1 arm"}}

    # ---------------------------------------------------------------- write
    meta = {
        "script": os.path.abspath(__file__), "git_commit": git_commit(),
        "generated": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "args": {k: v for k, v in vars(args).items()},
        "runs": runs, "spaces": args.spaces, "eps": eps_list,
        "bootstrap": {"replicates": args.bootstrap, "seed": args.seed,
                      "draw": "score_folds.cluster_bootstrap verbatim: default_rng(seed)."
                              "integers(0, K, size=(B, K)) over the stratum's sorted "
                              "cluster keys, percentile 95% CI",
                      "unit": "PRIMARY holdout group (cells group_id == fold_plan "
                              "genome_group) -> fields lo/hi; sensitivity source genome "
                              "(ci_src) and genus (ci_genus). Secondary (as shipped): "
                              "primary genus (truth_lineage), sensitivity source (ci_src).",
                      "sesoi_delta_f1": G.get("SESOI")},
        "population": "score_folds scored cells (cells.csv.gz space=ec, outcome=scored); "
                      "abstentions have no p",
        "metric_keys": {
            "union": "denominator = the union's own functions (p>0)",
            "universe": "denominator = full fixed vocabulary, p=0 as a real bin",
            "brier_union": "Brier of the binary union (= 1 - precision on the union denominator)",
            "brier_const_c": "every union member gets c = stratum pooled precision (in-sample)",
            "oof_*": "out-of-fold recalibration of p within the stratum (4 folds train, 5th test)",
            "calib_*": "logistic of y on logit(clip(p, eps, 1-eps)); perfect = intercept 0, slope 1",
            "murphy_*": "Brier = REL - RES + UNC over distinct p values (exact)",
            "binary": "rules union (p>0), p>=0.5, p>0.5, p=1; micro pools tp/fp/fn, macro "
                      "averages per-cell P/R/F1 (empty prediction: P = 0, F1 = 0)"},
        "layer_a_jobs": jobstats, "layer_b_jobs": bstats,
        "seconds_total": round(time.time() - T0, 1),
        "parent_peak_rss_gb": round(rss_gb(), 2),
        "children_peak_rss_gb": round(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1e6, 2)}
    scores = {"meta": meta, "primary_controlled": primary,
              "secondary_as_shipped": secondary, "exploratory_topk": topk_scores,
              "exploratory_topk_note": "EXPLORATORY, NOT PRE-REGISTERED: top-k identity-"
                                       "weighted prevalence from the controlled hits files"}
    if args.cache_dir:
        os.makedirs(args.cache_dir, exist_ok=True)
        with open(os.path.join(args.cache_dir, "final.pkl"), "wb") as fh:
            pickle.dump((rel_all, scores, assertions), fh, protocol=5)
    write_outputs(args, rel_all, scores, assertions, T0)


REL_COLS = ["analysis", "run", "mode", "density", "forced", "space", "pop", "n_bucket",
            "novelty_bin", "rank1_rank", "region", "stratum_kind", "stratum_value",
            "denominator", "bin_scheme", "bin", "p_lo", "p_hi", "n_points", "n_present",
            "n_cells", "n_src", "n_distinct_pairs", "mean_p", "obs_freq", "ci_unit", "ci_lo", "ci_hi",
            "ci_lo_src", "ci_hi_src", "ci_lo_genus", "ci_hi_genus"]


def write_outputs(args, rel_all, scores, assertions, T0):
    os.makedirs(args.outdir, exist_ok=True)
    rel_df = pd.DataFrame(rel_all)
    # the spec's explicit stratum columns (one crossing at a time; 'all' = not stratified)
    rel_df["n_bucket"] = rel_df["pop"]
    for k in ("novelty_bin", "rank1_rank", "region"):
        rel_df[k] = np.where(rel_df["stratum_kind"] == k, rel_df["stratum_value"], "all")
    for c in REL_COLS:
        if c not in rel_df.columns:
            rel_df[c] = np.nan
    rel_df = rel_df[REL_COLS]
    for c in ("p_lo", "p_hi", "mean_p", "obs_freq", "ci_lo", "ci_hi", "ci_lo_src",
              "ci_hi_src", "ci_lo_genus", "ci_hi_genus"):
        rel_df[c] = pd.to_numeric(rel_df[c], errors="coerce").round(6)
    for c in ("n_cells", "n_src", "n_distinct_pairs"):
        rel_df[c] = pd.to_numeric(rel_df[c], errors="coerce").astype("Int64")
    rel_df.to_csv(os.path.join(args.outdir, "reliability.csv"), index=False)
    scores["meta"]["seconds_total_at_write"] = round(time.time() - T0, 1)
    with open(os.path.join(args.outdir, "scores.json"), "wb") as fh:
        fh.write(orjson.dumps(tojson(scores), option=orjson.OPT_INDENT_2))
    with open(os.path.join(args.outdir, "assertions.json"), "wb") as fh:
        fh.write(orjson.dumps(tojson(assertions), option=orjson.OPT_INDENT_2))
    write_readme(os.path.join(args.outdir, "README.md"), tojson(scores), tojson(assertions),
                 args, rel_df)
    log(f"wrote {args.outdir}: reliability.csv ({len(rel_df)} rows), scores.json, "
        f"assertions.json, README.md; all required assertions zero: "
        f"{assertions['all_required_zero']}")
    for k, v in assertions["required"].items():
        log(f"  assert {k}: {v}")
    if not assertions["all_required_zero"]:
        raise SystemExit("[calibrate] a required assertion is non-zero: the outputs above are "
                         "not trustworthy (see assertions.json)")

if __name__ == "__main__":
    main()
