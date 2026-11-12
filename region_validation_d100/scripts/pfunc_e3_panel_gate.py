#!/usr/bin/env python
"""Validation gate for the production E3 estimator (pfunc_estimator.py) on the 10k panel.

PASS CRITERION (pre-declared, INTEGRATION_SPEC.md "Validation gate"; decided before any
result is computed): under self-exclusion at 0.987 (exclude_0.987, the primary mode), in the
n >= 2 stratum (n = number of reps of the answered re-selection), production E3 with species
dedup (e3_species, the production setting) must have
  (a) a LOWER pooled union Brier than k/n, with the 95% cluster-bootstrap CI of
      Brier(e3_species) - Brier(kn) entirely below 0, AND
  (b) a HIGHER pooled union AUROC than k/n, with the 95% CI of
      AUROC(e3_species) - AUROC(kn) entirely above 0.
Both must hold. Otherwise FAIL (E3 ships opt-in only). exclude_0.995, n = 1, "all",
e3_genome and the per-region numbers are reported but do not enter the verdict.

What is scored (per amplicon of region_validation_d100, 67,991 amplicons):
  * self-exclusion re-selection exactly as score_genecap / genecap_tau0995:
    score_genecap.exclude_self_record(hit_record, src, src_taxid, tau) ->
    select_references.select_representatives(filtered, gene_provider=<slim PGFam store>,
    knobs=score_genecap._select_knobs()) -> score_genecap.score_one(...) for the outcome
    (the domain argument is left blank so taxopy is never loaded);
  * only outcome == 'ok' amplicons are scored (truth present and non-empty, every rep in
    the cache, >= 1 rep, non-empty union);
  * cells = the union U of the reps' PGFams (k/n's support); truth y = 1[f in PGF(source)];
  * k/n(f) = (# reps carrying f) / n, n = len(reps) (gene-less reps count in n);
  * E3(f) = pfunc_estimator.e3_probabilities({f: {f}}, neighbours, members, store) with
      neighbours = pfunc_estimator.neighbours_from_hits(FILTERED record, reps=selected dicts,
                   dedup='species' (e3_species) | 'genome' (e3_genome), reps_first=True),
      members    = pfunc_estimator.clade_members(selected[0], catalogue,
                   exclude_species={source's NCBI species taxid}, exclude_genomes={source}),
    and k/n kept for the amplicon when E3 returns no estimate (production behaviour).
Metrics per mode x stratum (n = 1, n >= 2, all) x estimator, pooled over cells:
  Brier, AUROC (ties 1/2; point from a 1,000,001-bin histogram, CI from 1,001-bin
  histograms as Stage E), logistic calibration slope (truth ~ logit(clip(q, 1e-3, 1-1e-3)),
  fitted on the 1,001-bin grouped data, point and CI), Stage E's OLS ("WLS") slope (exact),
  mean q, prevalence. Deltas vs k/n: cluster bootstrap over source genomes, B = 2,000,
  seed 1729, calibrate_pfunc.boot_weights' draw verbatim; all estimators of a stratum share
  the draws (paired). Per region (n >= 2): the same, clusters = amplicons (one per source).

Memory/IO: selection.json and the hits JSON are streamed with ijson (never json.load); the
shared PGFam cache is streamed once into ONE concatenated int32 array (PGFam strings interned
to ints, ~1.3 GB for all ~160k genomes) and never written; workers fork after the load and
read it copy-on-write. The ~130M cells per mode are never materialised: workers return
per-amplicon sufficient statistics and histograms.

Outputs: data/pfunc_e3_gate.json (summary), data/pfunc_e3_gate/per_amplicon.csv.gz,
data/pfunc_e3_gate/run_meta.json; with --limit N everything goes to data/pfunc_e3_gate/smoke/.

Interpreter: ~/Documents/py_venv/bin/python. Run under systemd-run (MemoryMax=40G).
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import multiprocessing as mp
import os
import resource
import sys
import threading
import time
from collections import Counter, OrderedDict

import ijson
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import PATHS  # noqa: E402
import score_genecap as S  # noqa: E402

sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import select_references as SR  # noqa: E402
import pfunc_estimator as P  # noqa: E402

TAUS = (0.987, 0.995)
MODES = ("exclude_0.987", "exclude_0.995")
PRIMARY_MODE = "exclude_0.987"
ESTS = ("kn", "e3_species", "e3_genome")
DEDUP = {"e3_species": "species", "e3_genome": "genome"}
STRATA = ("n1", "n2+", "all")
NBIN = 1001            # bootstrap histogram bins: round(q * 1000)
FINE = 1_000_000       # point-estimate histogram: floor(q * 1e6), q = 1 -> 1e6
B_BOOT = 2000
SEED = 1729
EPS_LOGIT = 1e-3
SET_CACHE = 1500       # per-worker LRU of genome -> frozenset(int PGFam ids)

PASS_CRITERION = (
    "PASS iff, under exclude_0.987 in the n >= 2 stratum, e3_species (production E3, species "
    "dedup) has a lower pooled union Brier than k/n with the 95% cluster-bootstrap CI of "
    "Brier(e3_species) - Brier(kn) entirely below 0, AND a higher pooled union AUROC than k/n "
    "with the 95% CI of AUROC(e3_species) - AUROC(kn) entirely above 0 (bootstrap over source "
    "genomes, B = 2000, seed 1729). Otherwise FAIL: E3 ships opt-in only. Declared before the "
    "run (INTEGRATION_SPEC.md 'Validation gate').")

DATA = os.path.dirname(PATHS.selection_out)
OUT_JSON = os.path.join(DATA, "pfunc_e3_gate.json")
OUT_DIR = os.path.join(DATA, "pfunc_e3_gate")
REF_CELLS = {"exclude_0.987": PATHS.genecap_cells_csv,
             "exclude_0.995": os.path.join(DATA, "gene_capture_cells_tau0995.csv")}


def log(m):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576
    print("[e3-gate %s | main maxrss %.1f GB] %s" % (time.strftime("%H:%M:%S"), rss, m), flush=True)


def cgroup_mem():
    """memory.current / memory.peak of this process's cgroup (the systemd unit), GB."""
    try:
        rel = open("/proc/self/cgroup").read().strip().split("::", 1)[1]
        base = "/sys/fs/cgroup" + rel
        out = {}
        for f in ("memory.current", "memory.peak"):
            p = os.path.join(base, f)
            if os.path.exists(p):
                out[f] = round(int(open(p).read().strip()) / 1e9, 2)
        return out
    except Exception:
        return {}


# =========================================================================== PGFam store
class IntStore:
    """genome_id -> sorted unique int32 PGFam ids, all genomes in one array."""

    def __init__(self):
        self.vocab = {}
        self.off = {}
        self.data = None
        self.n_empty = 0

    def fill(self, path):
        chunks, pos = [], 0
        vocab = self.vocab
        t0 = time.time()
        with open(path, "rb") as fh:
            for i, (gid, fams) in enumerate(ijson.kvitems(fh, "", use_float=True), 1):
                if fams:
                    a = np.fromiter((vocab.setdefault(f, len(vocab)) for f in fams),
                                    dtype=np.int32, count=len(fams))
                    a = np.unique(a)
                else:
                    a = np.zeros(0, np.int32)
                    self.n_empty += 1
                self.off[gid] = (pos, pos + len(a))
                pos += len(a)
                chunks.append(a)
                if i % 20000 == 0:
                    log("  pgfam cache: %d genomes, %d ids, %d distinct PGFams (%.0fs)"
                        % (i, pos, len(vocab), time.time() - t0))
        self.data = np.concatenate(chunks)
        del chunks

    def arr(self, gid):
        o = self.off.get(str(gid))
        return None if o is None else self.data[o[0]:o[1]]

    def __contains__(self, gid):
        return str(gid) in self.off


class LRUSets:
    """store.get(gid) -> frozenset of int ids (None when absent); per-worker LRU. Called as a
    gene provider (SR / score_one) it returns set() for an absent genome and records it in
    .misses; E3's .get lookups of absent neighbours are not provider misses."""

    def __init__(self, store, maxsize=SET_CACHE):
        self.store, self.maxsize = store, maxsize
        self.c = OrderedDict()
        self.misses = set()

    def get(self, gid):
        gid = str(gid)
        s = self.c.get(gid)
        if s is not None:
            self.c.move_to_end(gid)
            return s
        a = self.store.arr(gid)
        if a is None:              # E3 neighbour absent from the cache: counted by the caller
            return None
        s = frozenset(a.tolist())
        self.c[gid] = s
        if len(self.c) > self.maxsize:
            self.c.popitem(last=False)
        return s

    def __call__(self, gid):           # gene_provider protocol (SR / score_one): set or empty
        s = self.get(gid)
        if s is None:                  # a provider miss (genecap's fixed-point signal)
            self.misses.add(str(gid))
            return set()
        return s


# =========================================================================== worker
G = {}   # worker globals (set before fork)


def _sums(q, y):
    yf = y.astype(np.float64)
    return np.array([((q - yf) ** 2).sum(), q.sum(), (q * q).sum(), (q * yf).sum()])


def _hist(q, y):
    b = np.clip(np.floor(q * 1000.0 + 0.5), 0, 1000).astype(np.int64)
    return np.stack([np.bincount(b[y], minlength=NBIN), np.bincount(b[~y], minlength=NBIN)])


def _fine(q, y):
    f = np.minimum(np.floor(q * FINE), FINE).astype(np.int32)
    out = []
    for m in (y, ~y):
        u, c = np.unique(f[m], return_counts=True)
        out.append((u, c.astype(np.int64)))
    return out


def auroc_exact(q, y):
    if len(q) == 0 or y.all() or not y.any():
        return float("nan")
    uq, inv = np.unique(q, return_inverse=True)
    Hp = np.bincount(inv, weights=y.astype(np.float64), minlength=len(uq))
    Hn = np.bincount(inv, weights=(~y).astype(np.float64), minlength=len(uq))
    return float(auroc_from_hist(Hp, Hn))


def q_independent(U, nb, mem, store, alpha=P.ALPHA, min_clade=P.MIN_CLADE):
    """Independent numpy recomputation of the E3 formula (PE:742-745) for the smoke check."""
    num = np.zeros(len(U))
    den = 0.0
    for g, w in nb:
        a = store.arr(g)
        if a is None or len(a) == 0:
            continue
        num += w * np.isin(U, a, assume_unique=True)
        den += w
    used = [store.arr(m) for m in mem]
    used = [a for a in used if a is not None and len(a)]
    if len(used) >= min_clade:
        pi = np.zeros(len(U))
        for a in used:
            pi += np.isin(U, a, assume_unique=True)
        num += alpha * pi / len(used)
        den += alpha
    if den <= 0:
        return None
    q = num / den
    q[q >= 1 - 1e-9] = 1.0
    q[q <= 1e-12] = 0.0
    return q


def work(item):
    key, hr = item
    store, sets, cat = G["store"], G["sets"], G["cat"]
    cache_keys = G["cache_keys"]
    src_gid, region = S.split_key(key)
    src_tid = S._src_taxon_id(src_gid)
    out = []
    for tau, mode in zip(TAUS, MODES):
        t0 = time.perf_counter()
        excl_rec, n_removed = S.exclude_self_record(hr, src_gid, src_tid, tau)
        sel = SR.select_representatives(excl_rec, gene_provider=sets, knobs=S._select_knobs())
        cell = S.score_one(key, mode, sel, src_gid, src_tid, region, "", sets, cache_keys,
                           available=True)
        row = dict(key=key, src_genome_id=src_gid, region=region, mode=mode,
                   outcome=cell["outcome"], tier=sel.get("tier") or "", n=cell["n_reps"],
                   n_selected=cell["n_selected"], n_geneless=cell["n_reps_geneless"],
                   n_removed_self=n_removed, truth_size=cell["truth_size"],
                   union_size=cell["pred_size"], tp=cell["intersection"], scored=False)
        res = {"row": row}
        if cell["outcome"] != S.OUTCOME_OK:
            out.append(res)
            continue
        n = cell["n_reps"]
        if G["n1_mod"] and n == 1 and int(hashlib.md5(key.encode()).hexdigest(), 16) % G["n1_mod"]:
            row["sampled_out"] = True
            out.append(res)
            continue
        reps = S.selected_rep_gids(sel)
        arrs = [store.arr(g) for g in reps]
        U, cnt = np.unique(np.concatenate(arrs), return_counts=True)
        T = store.arr(src_gid)
        y = np.isin(U, T, assume_unique=True)
        row["check_union"] = int(len(U) == cell["pred_size"] and int(y.sum()) == cell["intersection"]
                                 and len(T) == cell["truth_size"] and len(reps) == n)
        kn = cnt / float(n)
        sp_src, _f, _o = cat.ranks_for(src_gid)
        rank, mem = P.clade_members(reps[0], cat,
                                    exclude_species={sp_src} if sp_src is not None else None,
                                    exclude_genomes={src_gid})
        a_sp, a_fam, a_ord = cat.ranks_for(reps[0])
        row.update(anchor=reps[0], anchor_in_catalogue=int(reps[0] in cat.genomes),
                   anchor_family=a_fam or "", anchor_order=a_ord or "",
                   src_species_known=int(sp_src is not None),
                   clade_rank=rank or "", n_members=len(mem))
        Ul = U.tolist()
        feats = {f: {f} for f in Ul}
        qs = {"kn": kn}
        missing = set()
        t1 = time.perf_counter()
        for est in ("e3_species", "e3_genome"):
            info = {}
            nb = P.neighbours_from_hits(excl_rec, reps=sel["selected"], dedup=DEDUP[est],
                                        reps_first=True, info=info)
            q, diag = P.e3_probabilities(feats, nb, mem, sets)
            sfx = "_sp" if est == "e3_species" else "_g"
            if q:
                qa = np.fromiter((q[f] for f in Ul), dtype=np.float64, count=len(Ul))
            else:
                qa = kn.copy()
            for g, _w in nb:
                if g not in store:
                    missing.add(g)
            row.update({
                "n_stage0" + sfx: info.get("n_stage0"), "n_nb" + sfx: diag["n_neighbours"],
                "n_nb_used" + sfx: diag["n_neighbours_used"],
                "n_nb_missing" + sfx: diag["n_neighbours_missing"],
                "n_nb_empty" + sfx: diag["n_neighbours_empty"],
                "w_total" + sfx: diag["weight_total"], "w_missing" + sfx: diag["weight_missing"],
                "fallback_reps" + sfx: int(bool(info.get("fallback_reps"))),
                "clade_used" + sfx: int(diag["clade_used"]),
                "n_members_used" + sfx: diag["n_members_used"],
                "no_estimate" + sfx: int(diag["no_estimate"]),
            })
            if G["verify_q"]:
                qi = q_independent(U, nb, mem, store)
                d = (float(np.abs(qi - qa).max()) if (qi is not None and q) else
                     (0.0 if (qi is None and not q) else float("inf")))
                row["verify_maxdiff" + sfx] = d
            qs[est] = qa
        t2 = time.perf_counter()
        sums = np.stack([_sums(qs[e], y) for e in ESTS])
        hist = np.stack([_hist(qs[e], y) for e in ESTS]).astype(np.int32)
        fine = [_fine(qs[e], y) for e in ESTS]
        for e, tag in zip(ESTS, ("kn", "e3s", "e3g")):
            row["brier_" + tag] = sums[ESTS.index(e), 0] / len(U)
            row["meanq_" + tag] = sums[ESTS.index(e), 1] / len(U)
            row["auroc_" + tag] = auroc_exact(qs[e], y)
        row.update(scored=True, n_missing_gids=len(missing),
                   t_select_ms=round((t1 - t0) * 1000, 2), t_e3_ms=round((t2 - t1) * 1000, 2))
        res.update(sums=sums, hist=hist, fine=fine, ncell=len(U), npos=int(y.sum()),
                   missing=sorted(missing))
        if G["keep_cells"]:
            res["cells"] = (np.stack([qs[e] for e in ESTS]), y)
        out.append(res)
    if sets.misses:
        out[0]["provider_misses"] = sorted(sets.misses)
        sets.misses = set()
    return out


def _init_worker():
    G["sets"] = LRUSets(G["store"])


# =========================================================================== streaming
def stream_items(path):
    with open(path, "rb") as fh:
        for k, v in ijson.kvitems(fh, "", use_float=True):
            yield k, v


def first_hit_keys(n):
    out = []
    for k, v in stream_items(PATHS.hits_json):
        if k != "_meta" and isinstance(v, dict):
            out.append(k)
            if len(out) >= n:
                break
    return out


def selection_keys(limit_keys):
    keys, n_sel_include = set(), {}
    want = set(limit_keys) if limit_keys is not None else None
    for k, v in stream_items(PATHS.selection_out):
        if k == "_meta":
            continue
        if want is not None and k not in want:
            continue
        keys.add(k)
        n_sel_include[k] = (v or {}).get("n_selected") if isinstance(v, dict) else None
        if want is not None and len(keys) == len(want):
            break
    return keys, n_sel_include


def load_ref_cells(keys):
    ref = {}
    for mode, path in REF_CELLS.items():
        with open(path) as fh:
            for r in csv.DictReader(fh):
                if r["self_mode"] != mode or (keys is not None and r["ampliconKey"] not in keys):
                    continue
                ref[(r["ampliconKey"], mode)] = (r["outcome"], int(r["n_reps"]), int(r["truth_size"]),
                                                 int(r["pred_size"]), int(r["intersection"]))
    return ref


# =========================================================================== statistics
def boot_weights(K, B, seed):
    """calibrate_pfunc.boot_weights verbatim (score_folds' cluster-bootstrap draw)."""
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, K, size=(B, K))
    flat = (idx + (np.arange(B, dtype=np.int64) * K)[:, None]).ravel()
    return np.bincount(flat, minlength=B * K).reshape(B, K).astype(np.float64)


def auroc_from_hist(Hp, Hn):
    """pfunc_estimators.auroc_from_hist verbatim: last axis = bins, ascending q; ties 1/2."""
    above = np.cumsum(Hp[..., ::-1], -1)[..., ::-1] - Hp
    num = (Hn * (above + 0.5 * Hp)).sum(-1)
    with np.errstate(all="ignore"):
        return num / (Hp.sum(-1) * Hn.sum(-1))


_XL = np.log(np.clip(np.arange(NBIN) / 1000.0, EPS_LOGIT, 1 - EPS_LOGIT)
             / (1 - np.clip(np.arange(NBIN) / 1000.0, EPS_LOGIT, 1 - EPS_LOGIT)))


def logistic_slope(Hp, Hn, x=_XL, iters=60):
    """Grouped logistic regression y ~ a + s*x per row (vectorised Newton). Returns (a, s);
    NaN where x has no weighted variance (e.g. k/n at n = 1, every cell at q = 1)."""
    Hp = np.atleast_2d(np.asarray(Hp, np.float64))
    Hn = np.atleast_2d(np.asarray(Hn, np.float64))
    Nt = Hp + Hn
    W = Nt.sum(-1)
    mx = (Nt * x).sum(-1) / W
    vx = (Nt * (x - mx[:, None]) ** 2).sum(-1) / W
    ok = vx > 1e-12
    pbar = np.clip(Hp.sum(-1) / W, 1e-12, 1 - 1e-12)
    a = np.log(pbar / (1 - pbar))
    s = np.zeros(len(W))
    for _ in range(iters):
        eta = a[:, None] + s[:, None] * x
        with np.errstate(over="ignore"):
            p = 1 / (1 + np.exp(-eta))
        r = Hp - Nt * p
        ga, gs = r.sum(-1), (r * x).sum(-1)
        v = Nt * p * (1 - p)
        haa, has, hss = v.sum(-1), (v * x).sum(-1), (v * x * x).sum(-1)
        det = haa * hss - has * has
        with np.errstate(all="ignore"):
            da = np.where(ok, (hss * ga - has * gs) / det, 0.0)
            ds = np.where(ok, (haa * gs - has * ga) / det, 0.0)
        a, s = a + da, s + ds
        if np.nanmax(np.abs(np.r_[da, ds])) < 1e-10:
            break
    s = np.where(ok, s, np.nan)
    a = np.where(ok, a, np.nan)
    return a, s


def ci(v):
    v = np.asarray(v, np.float64)
    v = v[np.isfinite(v)]
    if len(v) == 0:
        return [None, None]
    return [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]


def fnum(x, nd=6):
    if x is None:
        return None
    x = float(x)
    return None if not math.isfinite(x) else round(x, nd)


def aggregate(rows_idx, cluster_ids, *arrays):
    """Sum per-amplicon arrays into per-cluster rows (clusters sorted by id)."""
    cid = cluster_ids[rows_idx]
    uc, inv = np.unique(cid, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    starts = np.flatnonzero(np.r_[True, inv[order][1:] != inv[order][:-1]])
    out = []
    for A in arrays:
        Ar = A[rows_idx][order].astype(np.float64)
        out.append(np.add.reduceat(Ar, starts, axis=0))
    return uc, out


def stratum_block(idx, clusters, SUMS, HIST, NCELL, NPOS, fine_pos, fine_neg, B, seed):
    """Metrics + paired cluster-bootstrap deltas for one set of amplicon slots."""
    if len(idx) == 0:
        return {"n_amplicons": 0}
    uc, (sums_c, hist_c, ncell_c, npos_c) = aggregate(idx, clusters, SUMS, HIST, NCELL, NPOS)
    K = len(uc)
    N, POS = ncell_c.sum(), npos_c.sum()
    W = boot_weights(K, B, seed) if (B and K >= 2) else None
    blk = {"n_amplicons": int(len(idx)), "n_sources": int(K), "n_cells": int(N),
           "n_pos": int(POS), "prevalence": fnum(POS / N), "estimators": {}, "delta_vs_kn": {}}
    if W is not None:
        Nb, POSb = W @ ncell_c, W @ npos_c
    boot = {}
    for ei, e in enumerate(ESTS):
        s = sums_c[:, ei, :].sum(0)
        Hp_c, Hn_c = hist_c[:, ei, 0, :], hist_c[:, ei, 1, :]
        Hp, Hn = Hp_c.sum(0), Hn_c.sum(0)
        brier = s[0] / N
        sxx = s[2] - s[1] ** 2 / N
        sxy = s[3] - s[1] * POS / N
        ols = sxy / sxx if sxx > 1e-12 * N else float("nan")
        a_l, s_l = logistic_slope(Hp, Hn)
        fp, fn_ = fine_pos[ei], fine_neg[ei]
        blk["estimators"][e] = {
            "brier": fnum(brier),
            "auroc": fnum(auroc_from_hist(np.asarray(fp, np.float64), np.asarray(fn_, np.float64))),
            "auroc_1001bin": fnum(auroc_from_hist(Hp, Hn)),
            "logistic_slope": fnum(s_l[0], 4), "logistic_intercept": fnum(a_l[0], 4),
            "ols_slope_stageE": fnum(ols, 4), "mean_q": fnum(s[1] / N),
        }
        if W is not None:
            sb = W @ sums_c[:, ei, :]
            Pb, Nn = W @ Hp_c, W @ Hn_c
            br = sb[:, 0] / Nb
            sxx = sb[:, 2] - sb[:, 1] ** 2 / Nb
            sxy = sb[:, 3] - sb[:, 1] * POSb / Nb
            with np.errstate(all="ignore"):
                olsb = np.where(sxx > 1e-12 * Nb, sxy / sxx, np.nan)
            boot[e] = {"brier": br, "auroc": auroc_from_hist(Pb, Nn),
                       "logistic_slope": logistic_slope(Pb, Nn)[1], "ols_slope_stageE": olsb,
                       "mean_q": sb[:, 1] / Nb}
            for m in ("brier", "auroc", "logistic_slope", "ols_slope_stageE", "mean_q"):
                blk["estimators"][e][m + "_ci95"] = [fnum(x) for x in ci(boot[e][m])]
    if W is not None:
        for e in ESTS[1:]:
            d = {}
            for m in ("brier", "auroc", "logistic_slope", "ols_slope_stageE", "mean_q"):
                pt_e, pt_k = blk["estimators"][e][m], blk["estimators"]["kn"][m]
                pt = None if (pt_e is None or pt_k is None) else pt_e - pt_k
                if m == "auroc":   # the bootstrap delta is on the 1,001-bin AUROC
                    pt_bin = blk["estimators"][e]["auroc_1001bin"] - blk["estimators"]["kn"]["auroc_1001bin"]
                    d["auroc_1001bin"] = fnum(pt_bin)
                lo, hi = ci(boot[e][m] - boot["kn"][m])
                d[m] = fnum(pt)
                d[m + "_ci95"] = [fnum(lo), fnum(hi)]
            d["brier_ci_excludes_0_below"] = bool(d["brier_ci95"][1] is not None and d["brier_ci95"][1] < 0)
            d["auroc_ci_excludes_0_above"] = bool(d["auroc_ci95"][0] is not None and d["auroc_ci95"][0] > 0)
            blk["delta_vs_kn"][e] = d
    blk["bootstrap"] = {"B": B, "seed": seed, "unit": "source genome" if clusters is not None else "",
                        "n_clusters": int(K)}
    return blk


# =========================================================================== main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--limit", type=int, default=0, help="smoke test on the first N hit records")
    ap.add_argument("--workers", type=int, default=12)
    ap.add_argument("--n1-sample-mod", type=int, default=0,
                    help="score only n=1 amplicons with md5(key) %% M == 0 (0 = all)")
    ap.add_argument("--bootstrap", type=int, default=B_BOOT)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--verify-q", action="store_true",
                    help="recompute every E3 q with an independent numpy formula (smoke)")
    args = ap.parse_args()

    out_json, out_dir = OUT_JSON, OUT_DIR
    if args.limit:
        out_dir = os.path.join(OUT_DIR, "smoke")
        out_json = os.path.join(out_dir, "pfunc_e3_gate.json")
    os.makedirs(out_dir, exist_ok=True)
    t_start = time.time()
    print("PASS CRITERION: " + PASS_CRITERION, flush=True)
    meta = OrderedDict(pass_criterion=PASS_CRITERION)
    meta["started"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    meta["inputs"] = {p: os.path.getsize(p) for p in (PATHS.selection_out, PATHS.hits_json, PATHS.pgfam_cache)}
    meta["args"] = vars(args)

    # ---------------------------------------------------------------- catalogue
    log("loading the E3 clade catalogue")
    cat = P.load_clade_catalogue()
    cat_current = P.catalogue_is_current(cat)
    meta["catalogue"] = {k: v for k, v in cat.meta.items() if k != "empty_pgfam_genomes"}
    meta["catalogue"]["is_current"] = cat_current
    if not cat_current:
        raise SystemExit("clade catalogue is stale against the shared PGFam cache; rebuild first")
    P.clade_members(next(iter(cat.genomes)), cat)          # build the index before fork

    # ---------------------------------------------------------------- keys
    limit_keys = None
    if args.limit:
        limit_keys = first_hit_keys(args.limit)
        log("smoke: %d hit keys" % len(limit_keys))
    log("streaming selection.json for the amplicon key set")
    sel_keys, n_sel_inc = selection_keys(limit_keys)
    log("%d selection keys" % len(sel_keys))

    # ---------------------------------------------------------------- PGFam store
    log("streaming the shared PGFam cache into an int32 store (read-only)")
    store = IntStore()
    store.fill(PATHS.pgfam_cache)
    cache_keys = set(store.off)
    log("store: %d genomes (%d empty), %d PGFam ids, %d distinct PGFams, %.2f GB"
        % (len(store.off), store.n_empty, len(store.data), len(store.vocab), store.data.nbytes / 1e9))
    meta["store"] = {"n_genomes": len(store.off), "n_empty": store.n_empty,
                     "n_ids": int(len(store.data)), "n_distinct_pgfams": len(store.vocab)}
    store.vocab = None   # not needed by the workers

    ref = load_ref_cells(set(limit_keys) if limit_keys else None)
    log("reference gene-capture cells: %d (key, mode) rows" % len(ref))

    # ---------------------------------------------------------------- accumulators
    nmax = len(sel_keys)
    acc = {}
    for mode in MODES:
        acc[mode] = dict(
            SUMS=np.zeros((nmax, 3, 4)), HIST=np.zeros((nmax, 3, 2, NBIN), np.int32),
            NCELL=np.zeros(nmax, np.int64), NPOS=np.zeros(nmax, np.int64),
            N=np.zeros(nmax, np.int32), SRC=np.zeros(nmax, np.int64), REG=np.zeros(nmax, np.int16),
            FINE=np.zeros((2, 3, 2, FINE + 1), np.int64), n=0)
    src_index, reg_index = {}, {}
    rows, missing_gids, prov_misses = [], Counter(), set()
    repro = Counter()
    repro_examples = []
    keep_cells = bool(args.limit)
    cells_keep = {m: {s: [] for s in ("n1", "n2+")} for m in MODES}

    G.update(store=store, cat=cat, cache_keys=cache_keys, n1_mod=args.n1_sample_mod,
             verify_q=args.verify_q, keep_cells=keep_cells)
    sem = threading.BoundedSemaphore(max(64, args.workers * 16))

    def feed():
        n = 0
        want = set(limit_keys) if limit_keys else None
        for k, v in stream_items(PATHS.hits_json):
            if k == "_meta" or not isinstance(v, dict) or k not in sel_keys:
                continue
            if want is not None and k not in want:
                continue
            sem.acquire()
            yield (k, v)
            n += 1
            if want is not None and n >= len(want):
                return

    log("scoring with %d workers" % args.workers)
    ctx = mp.get_context("fork")
    n_done, t_loop = 0, time.time()
    seen_keys = set()
    with ctx.Pool(args.workers, initializer=_init_worker) as pool:
        for out in pool.imap_unordered(work, feed(), chunksize=2):
            sem.release()
            n_done += 1
            if "provider_misses" in out[0]:
                prov_misses.update(out[0]["provider_misses"])
            for res in out:
                row = res["row"]
                key, mode = row["key"], row["mode"]
                seen_keys.add(key)
                r = ref.get((key, mode))
                mine = (row["outcome"], row["n"], row["truth_size"], row["union_size"], row["tp"])
                if r is None:
                    repro["no_reference_row"] += 1
                elif r == mine:
                    repro["match"] += 1
                else:
                    repro["mismatch"] += 1
                    if len(repro_examples) < 10:
                        repro_examples.append({"key": key, "mode": mode, "ref": r, "mine": mine})
                if row.get("scored"):
                    if not row.get("check_union"):
                        repro["union_check_failed"] += 1
                    a = acc[mode]
                    i = a["n"]
                    a["n"] += 1
                    row["slot"] = i
                    a["SUMS"][i] = res["sums"]
                    a["HIST"][i] = res["hist"]
                    a["NCELL"][i], a["NPOS"][i], a["N"][i] = res["ncell"], res["npos"], row["n"]
                    a["SRC"][i] = src_index.setdefault(row["src_genome_id"], len(src_index))
                    a["REG"][i] = reg_index.setdefault(row["region"], len(reg_index))
                    st = 0 if row["n"] == 1 else 1
                    for ei in range(3):
                        for pn in range(2):
                            u, c = res["fine"][ei][pn]
                            a["FINE"][st, ei, pn][u] += c
                    for g in res["missing"]:
                        missing_gids[g] += 1
                    if keep_cells:
                        cells_keep[mode]["n1" if row["n"] == 1 else "n2+"].append(res["cells"])
                rows.append(row)
            if n_done % 2000 == 0:
                el = time.time() - t_loop
                log("%d amplicons (%.1f/s), slots %s, cgroup %s"
                    % (n_done, n_done / el, {m: acc[m]["n"] for m in MODES}, cgroup_mem()))
    log("scoring done: %d amplicons in %.0fs" % (n_done, time.time() - t_loop))
    meta["n_selection_keys"] = len(sel_keys)
    meta["n_amplicons_streamed"] = n_done
    meta["selection_keys_without_hit_record"] = len(sel_keys - seen_keys)
    meta["provider_misses"] = len(prov_misses)
    meta["reproduction_vs_gene_capture_cells"] = {"counts": dict(repro), "examples": repro_examples,
                                                  "fields": "outcome, n_reps, truth_size, pred_size, intersection"}
    meta["neighbours_absent_from_cache"] = {"distinct_genomes": len(missing_gids),
                                            "amplicon_mode_occurrences": int(sum(missing_gids.values()))}
    if repro["mismatch"] or repro["union_check_failed"] or prov_misses:
        log("WARNING: reproduction mismatches %s, provider misses %d" % (dict(repro), len(prov_misses)))

    # ---------------------------------------------------------------- per-amplicon table
    cols = ["key", "src_genome_id", "region", "mode", "outcome", "tier", "n", "n_selected",
            "n_geneless", "n_removed_self", "truth_size", "union_size", "tp", "scored",
            "sampled_out", "check_union", "anchor", "anchor_in_catalogue", "anchor_family",
            "anchor_order", "src_species_known", "clade_rank", "n_members"]
    for sfx in ("_sp", "_g"):
        cols += [c + sfx for c in ("n_stage0", "n_nb", "n_nb_used", "n_nb_missing", "n_nb_empty",
                                   "w_total", "w_missing", "fallback_reps", "clade_used",
                                   "n_members_used", "no_estimate")]
    cols += ["brier_kn", "brier_e3s", "brier_e3g", "meanq_kn", "meanq_e3s", "meanq_e3g",
             "auroc_kn", "auroc_e3s", "auroc_e3g", "n_missing_gids", "t_select_ms", "t_e3_ms"]
    if args.verify_q:
        cols += ["verify_maxdiff_sp", "verify_maxdiff_g"]
    rows.sort(key=lambda r: (r["key"], r["mode"]))
    csv_path = os.path.join(out_dir, "per_amplicon.csv.gz")
    with gzip.open(csv_path, "wt", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for r in rows:
            vals = []
            for c in cols:
                v = r.get(c, "")
                if isinstance(v, float):
                    v = "" if not math.isfinite(v) else ("%.6g" % v)
                elif isinstance(v, bool):
                    v = int(v)
                vals.append(v)
            w.writerow(vals)
    log("wrote %s (%d rows)" % (csv_path, len(rows)))

    # ---------------------------------------------------------------- metrics
    # Deterministic cluster order: results arrive in imap_unordered order, so re-index the
    # sources (bootstrap clusters) and regions by sorted name before any draw is applied;
    # with a fixed seed the CIs are then reproducible run to run.
    src_rank = np.empty(len(src_index), np.int64)
    for r_, g_ in enumerate(sorted(src_index)):
        src_rank[src_index[g_]] = r_
    reg_rank = np.empty(len(reg_index), np.int16)
    for r_, g_ in enumerate(sorted(reg_index)):
        reg_rank[reg_index[g_]] = r_
    for mode in MODES:
        a = acc[mode]
        a["SRC"][:a["n"]] = src_rank[a["SRC"][:a["n"]]]
        a["REG"][:a["n"]] = reg_rank[a["REG"][:a["n"]]]
    log("metrics + cluster bootstrap (B=%d, seed %d)" % (args.bootstrap, args.seed))
    regions = sorted(reg_index)
    results, per_region, diag = OrderedDict(), OrderedDict(), OrderedDict()
    for mode in MODES:
        a = acc[mode]
        n_sc = a["n"]
        Nn = a["N"][:n_sc]
        masks = {"n1": np.flatnonzero(Nn == 1), "n2+": np.flatnonzero(Nn >= 2),
                 "all": np.arange(n_sc)}
        results[mode] = OrderedDict()
        for st in STRATA:
            idx = masks[st]
            if st == "all":
                fp = a["FINE"][:, :, 0].sum(0)
                fn_ = a["FINE"][:, :, 1].sum(0)
            else:
                k = 0 if st == "n1" else 1
                fp, fn_ = a["FINE"][k, :, 0], a["FINE"][k, :, 1]
            results[mode][st] = stratum_block(idx, a["SRC"], a["SUMS"], a["HIST"], a["NCELL"], a["NPOS"],
                                              fp, fn_, args.bootstrap, args.seed)
            log("  %s %s done" % (mode, st))
        # per region, n >= 2 (clusters = amplicons: one per source within a region)
        per_region[mode] = OrderedDict()
        for ri, reg in enumerate(regions):
            idx = masks["n2+"][a["REG"][masks["n2+"]] == ri]
            if len(idx) == 0:
                continue
            Hs = a["HIST"][idx].astype(np.int64)
            fp = Hs[:, :, 0, :].sum(0)
            fn_ = Hs[:, :, 1, :].sum(0)
            assert len(np.unique(a["SRC"][idx])) == len(idx), "one amplicon per source per region"
            blk = stratum_block(idx, a["SRC"], a["SUMS"], a["HIST"], a["NCELL"], a["NPOS"],
                                fp, fn_, args.bootstrap, args.seed)
            blk["note"] = "AUROC point estimates here use the 1,001-bin histogram"
            per_region[mode][reg] = blk
        # diagnostics over scored rows
        sc = [r for r in rows if r["mode"] == mode and r.get("scored")]
        d = OrderedDict()
        d["outcomes"] = dict(Counter(r["outcome"] for r in rows if r["mode"] == mode))
        d["n_scored"] = len(sc)
        d["n_sampled_out_n1"] = sum(1 for r in rows if r["mode"] == mode and r.get("sampled_out"))
        for st, f in (("n1", lambda r: r["n"] == 1), ("n2+", lambda r: r["n"] >= 2)):
            ss = [r for r in sc if f(r)]
            if not ss:
                continue
            dd = OrderedDict()
            dd["clade_rank_before_store"] = dict(Counter(r["clade_rank"] or "none" for r in ss))
            dd["anchor_in_catalogue"] = sum(r["anchor_in_catalogue"] for r in ss)
            dd["anchor_without_family_and_order"] = sum(1 for r in ss if not r["anchor_family"] and not r["anchor_order"])
            dd["source_species_unknown"] = sum(1 for r in ss if not r["src_species_known"])
            dd["amplicons_with_geneless_rep"] = sum(1 for r in ss if r["n_geneless"])
            for sfx, e in (("_sp", "e3_species"), ("_g", "e3_genome")):
                wt = sum(r["w_total" + sfx] for r in ss)
                wm = sum(r["w_missing" + sfx] for r in ss)
                dd[e] = {
                    "clade_used": sum(r["clade_used" + sfx] for r in ss),
                    "fallback_reps": sum(r["fallback_reps" + sfx] for r in ss),
                    "no_estimate_kept_kn": sum(r["no_estimate" + sfx] for r in ss),
                    "median_neighbours": float(np.median([r["n_nb" + sfx] for r in ss])),
                    "median_neighbours_used": float(np.median([r["n_nb_used" + sfx] for r in ss])),
                    "pooled_weight_missing_share": fnum(wm / wt if wt else None),
                    "amplicons_with_missing_neighbour": sum(1 for r in ss if r["n_nb_missing" + sfx]),
                }
                if args.verify_q:
                    vals = [r["verify_maxdiff" + sfx] for r in ss]
                    dd[e]["verify_q_maxdiff"] = max(vals) if vals else None
            dd["median_t_e3_ms_both_dedups"] = float(np.median([r["t_e3_ms"] for r in ss]))
            d[st] = dd
        diag[mode] = d

    # ---------------------------------------------------------------- exact AUROC check (smoke)
    if keep_cells:
        chk = {}
        for mode in MODES:
            for st in ("n1", "n2+"):
                L = cells_keep[mode][st]
                if not L:
                    continue
                Q = np.concatenate([c[0] for c in L], axis=1)
                Y = np.concatenate([c[1] for c in L])
                for ei, e in enumerate(ESTS):
                    ex = auroc_exact(Q[ei], Y)
                    rep = results[mode][st]["estimators"][e]
                    chk[f"{mode}|{st}|{e}"] = {"exact": fnum(ex), "fine_hist": rep["auroc"],
                                               "bin1001": rep["auroc_1001bin"],
                                               "brier_exact": fnum(float(((Q[ei] - Y) ** 2).mean())),
                                               "brier_reported": rep["brier"]}
        meta["smoke_exact_checks"] = chk

    # ---------------------------------------------------------------- verdict
    prim = results[PRIMARY_MODE]["n2+"]
    dsp = prim.get("delta_vs_kn", {}).get("e3_species", {})
    a_ok = bool(dsp.get("brier_ci_excludes_0_below"))
    b_ok = bool(dsp.get("auroc_ci_excludes_0_above"))
    valid = not (repro["mismatch"] or repro["union_check_failed"] or prov_misses)
    verdict = "PASS" if (a_ok and b_ok and valid) else "FAIL"
    if not valid:
        verdict += " (INVALID RUN: reproduction mismatch or provider misses)"
    if args.limit:
        verdict = "SMOKE ONLY (not a gate decision): " + verdict
    meta["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    meta["elapsed_s"] = round(time.time() - t_start)
    meta["main_maxrss_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576, 2)
    meta["cgroup_memory_gb"] = cgroup_mem()
    meta["definitions"] = {
        "cells": "(answered amplicon, PGFam in the union of the re-selected reps' PGFams)",
        "truth": "PGFam in the source genome's cached PGFam set",
        "kn": "share of the n reps carrying the PGFam (n = len(reps), gene-less reps included)",
        "e3_species": "pfunc_estimator.e3_probabilities, neighbours_from_hits(filtered record, reps "
                      "first, dedup species), clade_members(anchor, exclude source species taxid "
                      "+ source genome); empty-PGFam genomes excluded (empty_as_missing)",
        "e3_genome": "as e3_species with dedup genome (the Stage E/O validated setting)",
        "hit_exclusion": "score_genecap.exclude_self_record: genome == source OR prefix taxid == "
                         "source prefix taxid OR identity >= tau",
        "prior_exclusion": "every catalogue genome of the source's NCBI species taxid, and the source",
        "brier": "pooled mean (q - y)^2 over cells",
        "auroc": "pooled over cells, ties 1/2; point from a 1,000,001-bin histogram "
                 "(floor(q*1e6)); bootstrap CIs and deltas from 1,001-bin histograms "
                 "(round(q*1000)); delta point = fine-histogram difference",
        "logistic_slope": "y ~ a + s*logit(clip(q_bin, 1e-3, 1-1e-3)) on the 1,001-bin grouped "
                          "data (q rounded to 0.001), Newton; NaN when q has no variance",
        "ols_slope_stageE": "Stage E's union WLS slope = point-level OLS of y on q (exact)",
        "bootstrap": "calibrate_pfunc.boot_weights draw verbatim, B = %d, seed %d, unit = source "
                     "genome (per region: amplicon = source); percentile 95%% CI; paired "
                     "across estimators" % (args.bootstrap, args.seed),
        "strata": "n1: n = 1; n2+: n >= 2; all",
    }
    summary = OrderedDict()
    summary["_meta"] = meta
    summary["verdict"] = verdict
    summary["verdict_detail"] = {
        "mode": PRIMARY_MODE, "stratum": "n2+", "estimator": "e3_species",
        "delta_brier": dsp.get("brier"), "delta_brier_ci95": dsp.get("brier_ci95"),
        "delta_auroc": dsp.get("auroc"), "delta_auroc_ci95": dsp.get("auroc_ci95"),
        "criterion_a_brier": a_ok, "criterion_b_auroc": b_ok, "run_valid": valid}
    summary["results"] = results
    summary["per_region_n2+"] = per_region
    summary["diagnostics"] = diag
    with open(out_json, "w") as fh:
        json.dump(summary, fh, indent=1)
    with open(os.path.join(out_dir, "run_meta.json"), "w") as fh:
        json.dump(meta, fh, indent=1)
    log("wrote %s" % out_json)
    log("VERDICT: %s  %s" % (verdict, json.dumps(summary["verdict_detail"])))
    for mode in MODES:
        for st in STRATA:
            blk = results[mode][st]
            if not blk.get("n_amplicons"):
                continue
            for e in ESTS:
                v = blk["estimators"][e]
                log("  %s %-4s %-10s Brier %s AUROC %s slope %s meanq %s" % (
                    mode, st, e, v["brier"], v["auroc"], v["logistic_slope"], v["mean_q"]))
            for e in ESTS[1:]:
                dd = blk.get("delta_vs_kn", {}).get(e, {})
                log("  %s %-4s %-10s dBrier %s %s dAUROC %s %s" % (
                    mode, st, e, dd.get("brier"), dd.get("brier_ci95"), dd.get("auroc"), dd.get("auroc_ci95")))
    if not valid:
        sys.exit(3)


if __name__ == "__main__":
    main()
