#!/usr/bin/env python
"""Adversarial review of pfunc_e3_panel_gate.py: independent cell-level recompute of the
headline cell (exclude_0.987, n >= 2) plus leakage diagnostics. READ-ONLY on every repo input.

Independent of the gate's code for: self-exclusion filter (asserted equal to
score_genecap.exclude_self_record), stage 0, sort, species dedup with reps first, weights,
clade-prior membership (from the catalogue JSON), the E3 formula, k/n, Brier, and AUROC
(exact midrank Mann-Whitney over all cells, no histogram). Reused (production code, not under
review): select_references.select_representatives and _norm, score_genecap._select_knobs.

Leakage diagnostics per amplicon: neighbours (incl. reps) that are conspecific with the source
by NCBI species taxid (catalogue genomes, else catalogue taxa by prefix taxid) or by normalised
BV-BRC species string; neighbours whose genome had a hit removed by the self-exclusion filter
(multi-copy 16S reappearance); hard assertions (no neighbour == source genome, no neighbour with
the source prefix taxid, no neighbour hit identity >= tau, no prior member of the source species,
source genome not a prior member).

Variants: e3s (main), e3s_nocs (non-rep NCBI-conspecific neighbours dropped), e3s_nomc
(non-rep neighbours whose genome had a hit removed dropped), e3s_clean (both dropped).
"""
import gzip, csv, hashlib, json, math, multiprocessing as mp, os, sys, time, resource, threading
import numpy as np
import ijson

sys.path.insert(0, "/home/freiburger/Documents/prFBA/region_validation_d100/scripts")
sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import select_references as SR  # noqa
import score_genecap as S  # noqa

# Run 2026-11-12 by the gate's adversarial reviewer (outputs: indep_recompute.json, committed as
# data/pfunc_e3_gate/review_recompute.json). Set PFUNC_GATE_REVIEW_OUT to a scratch dir to rerun.
OUT = os.environ.get("PFUNC_GATE_REVIEW_OUT", "/home/freiburger/Documents/prFBA/region_validation_d100/data/pfunc_e3_gate/review")
CSV = "/home/freiburger/Documents/prFBA/region_validation_d100/data/pfunc_e3_gate/per_amplicon.csv.gz"
HITS = "/home/freiburger/Documents/prFBA/region_validation_d100/data/hits/asv_top20_alignment_hits.json"
CACHE = "/home/freiburger/Documents/prFBA/bvbrc_cache/genome_gene_families.json"
CATF = "/home/freiburger/Documents/prFBA/bvbrc_cache/e3_clade_catalogue.json"
TAU_X = 0.987
K, TAU_W, ALPHA, MIN_CLADE, MAXSP = 100, 0.02, 2.0, 3, 100
LIMIT = int(os.environ.get("LIMIT", "0"))
WORKERS = int(os.environ.get("WORKERS", "8"))
VARS = ("kn", "e3s", "e3s_nocs", "e3s_nomc", "e3s_clean")


def log(m):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576
    print("[review %s | %.1f GB] %s" % (time.strftime("%H:%M:%S"), rss, m), flush=True)


G = {}


def norm_sp(s):
    return SR._norm(s)


def taxid_prefix(g):
    try:
        return int(str(g).split(".")[0])
    except ValueError:
        return None


def ncbi_sp(g):
    v = G["cat_g"].get(str(g))
    if v is None:
        t = taxid_prefix(g)
        v = G["cat_t"].get(t) if t is not None else None
    return None if v is None else v[0]


def ncbi_ranks(g):
    v = G["cat_g"].get(str(g))
    if v is None:
        t = taxid_prefix(g)
        v = G["cat_t"].get(t) if t is not None else None
    return (None, None, None) if v is None else tuple(v)


def arr(g):
    o = G["off"].get(str(g))
    return None if o is None else G["data"][o[0]:o[1]]


class Prov:
    def __init__(self):
        self.c = {}

    def __call__(self, g):
        g = str(g)
        s = self.c.get(g)
        if s is None:
            a = arr(g)
            s = frozenset() if a is None else frozenset(a.tolist())
            if len(self.c) > 3000:
                self.c.clear()
            self.c[g] = s
        return s


def my_exclude(top, src, src_tid, tau):
    kept, removed = [], set()
    for h in top:
        g = str(h["genome_id"])
        t = h.get("taxon_id")
        if g == src or (src_tid is not None and t is not None and int(t) == src_tid) or h["identity"] >= tau:
            removed.add(g)
            continue
        kept.append(h)
    kept.sort(key=lambda h: (-h["identity"], -h["align_score"]))
    return kept, removed


def my_stage0(top, L):
    L = L or 1
    out = []
    for h in top:
        if round(min(h["aligned_len"] / L, 1.0), 4) < 0.90:
            continue
        e = h.get("edlib_identity")
        if e is not None and round(abs(h["identity"] - e), 4) > 0.02:
            continue
        if h["identity"] < 0.865:
            continue
        out.append(h)
    out.sort(key=lambda h: (-h["identity"], -h["align_score"]))
    return out


def work(item):
    key, hr = item
    src, region = key.split("__", 1)
    src_tid = taxid_prefix(src)
    top = hr.get("top20") or []
    kept, removed = my_exclude(top, src, src_tid, TAU_X)
    excl, _nrem = S.exclude_self_record(hr, src, S._src_taxon_id(src), TAU_X)
    same_filter = [str(h["genome_id"]) for h in kept] == [str(h["genome_id"]) for h in excl["top20"]]
    prov = G["prov"]
    sel = SR.select_representatives(excl, gene_provider=prov, knobs=S._select_knobs())
    reps, seen = [], set()
    for s in sel.get("selected") or []:
        g = str(s["genome_id"])
        if g not in seen:
            seen.add(g)
            reps.append((g, s))
    n = len(reps)
    rec = {"key": key, "src": src, "n": n, "same_filter": int(same_filter)}
    if n < 2:
        rec["skip"] = "n<2 (%d)" % n
        return rec, None
    ra = [arr(g) for g, _ in reps]
    if any(a is None for a in ra):
        rec["skip"] = "rep missing"
        return rec, None
    U, cnt = np.unique(np.concatenate(ra), return_counts=True)
    Tsrc = arr(src)
    y = np.isin(U, Tsrc, assume_unique=True)
    kn = cnt / float(n)

    # ---- neighbours: stage 0 on the filtered list (first 100), reps first, species dedup
    surv = my_stage0(excl["top20"][:100], hr.get("asv_len"))
    best_surv, best_any, sp_hit = {}, {}, {}
    for h in surv:
        g = str(h["genome_id"])
        best_surv.setdefault(g, float(h["identity"]))
        sp_hit.setdefault(g, (h.get("lineage") or {}).get("Species"))
    for h in excl["top20"]:
        g = str(h["genome_id"])
        x = float(h["identity"])
        if g not in best_any or x > best_any[g]:
            best_any[g] = x
        sp_hit.setdefault(g, (h.get("lineage") or {}).get("Species"))
    nb, keys, gs = [], set(), set()   # nb entries: (gid, identity, is_rep, bvbrc species)
    for g, s in reps:
        x = best_surv.get(g, best_any.get(g, s.get("identity")))
        sp = s.get("species") if s.get("species") is not None else sp_hit.get(g)
        gs.add(g)
        keys.add(norm_sp(sp) or "genome:" + g)
        nb.append((g, x, True, sp))
    for h in surv:
        if len(nb) >= K:
            break
        g = str(h["genome_id"])
        if g in gs:
            continue
        sp = (h.get("lineage") or {}).get("Species")
        k = norm_sp(sp) or "genome:" + g
        if k in keys:
            continue
        keys.add(k)
        gs.add(g)
        nb.append((g, float(h["identity"]), False, sp))
    nb = nb[:K]
    ids = [x for _, x, _, _ in nb if x is not None]
    id1 = max(ids)
    W = np.array([math.exp(-(id1 - (id1 if x is None else x)) / TAU_W) for _, x, _, _ in nb])

    # ---- source species (NCBI and BV-BRC string)
    sp_src, fam_src, ord_src = ncbi_ranks(src)
    src_str = None
    for h in top:
        if str(h["genome_id"]) == src or (src_tid is not None and h.get("taxon_id") is not None
                                           and int(h["taxon_id"]) == src_tid):
            src_str = (h.get("lineage") or {}).get("Species")
            if src_str:
                break
    src_key = norm_sp(src_str) if src_str else ""
    # ---- hard assertions on the neighbour list
    nb_ids = {g: x for g, x, _, _ in nb}
    hard = {"nb_is_source": sum(1 for g in nb_ids if g == src),
            "nb_src_prefix_taxid": sum(1 for g in nb_ids if src_tid is not None and taxid_prefix(g) == src_tid),
            "nb_identity_ge_tau": sum(1 for g, x, _, _ in nb if x is not None and x >= TAU_X)}
    # ---- usable neighbours and flags
    X, Wu, f_cs, f_mc, f_rep = [], [], [], [], []
    n_cs_ncbi = n_cs_str = n_mc = 0
    w_cs = w_mc = w_cs_str = 0.0
    rep_cs = 0
    for (g, x, isrep, sp), w in zip(nb, W):
        a = arr(g)
        cs = sp_src is not None and ncbi_sp(g) == sp_src
        cs_str = bool(src_key) and norm_sp(sp) == src_key
        mc = g in removed
        if isrep and cs:
            rep_cs += 1
        if a is None or len(a) == 0:
            continue
        X.append(np.isin(U, a, assume_unique=True))
        Wu.append(w)
        f_cs.append(cs and not isrep)
        f_mc.append(mc and not isrep)
        f_rep.append(isrep)
        if cs:
            n_cs_ncbi += 1
            w_cs += w
        if cs_str:
            n_cs_str += 1
            w_cs_str += w
        if mc:
            n_mc += 1
            w_mc += w
    Wu = np.array(Wu)
    X = np.array(X) if X else np.zeros((0, len(U)), bool)
    f_cs, f_mc = np.array(f_cs, bool), np.array(f_mc, bool)

    # ---- clade prior (independent of pfunc_estimator.clade_members)
    a_sp, a_fam, a_ord = ncbi_ranks(reps[0][0])
    rank, members = None, []
    for r, lab in (("family", a_fam), ("order", a_ord)):
        if lab is None:
            continue
        mem = []
        for sp, gl in G["idx"][r].get(lab, ()):
            if sp_src is not None and sp == sp_src:
                continue
            g = next((x for x in gl if x != src), None)
            if g is None:
                continue
            mem.append(g)
            if len(mem) >= MAXSP:
                break
        if len(mem) >= MIN_CLADE:
            rank, members = r, mem
            break
    hard["member_src_species"] = sum(1 for m in members if sp_src is not None and ncbi_sp(m) == sp_src)
    hard["member_is_source"] = int(src in members)
    mu = [arr(m) for m in members]
    mu = [a for a in mu if a is not None and len(a)]
    has = len(mu) >= MIN_CLADE
    pi = np.zeros(len(U))
    if has:
        for a in mu:
            pi += np.isin(U, a, assume_unique=True)
        pi /= len(mu)

    def e3(mask):
        w = Wu * mask
        num = (w[:, None] * X).sum(0) if len(w) else np.zeros(len(U))
        den = w.sum()
        if has:
            num = num + ALPHA * pi
            den = den + ALPHA
        if den <= 0:
            return kn.copy()
        q = num / den
        q[q >= 1 - 1e-9] = 1.0
        q[q <= 1e-12] = 0.0
        return q

    ones = np.ones(len(Wu), bool)
    Q = {"kn": kn, "e3s": e3(ones), "e3s_nocs": e3(~f_cs), "e3s_nomc": e3(~f_mc), "e3s_clean": e3(~(f_cs | f_mc))}
    yf = y.astype(np.float64)
    wt = float(Wu.sum())
    rec.update(N=int(len(U)), pos=int(y.sum()), clade_rank=rank or "", n_members=len(members),
               n_members_used=len(mu) if has else 0, n_nb=len(nb), n_nb_used=int(len(Wu)),
               sp_src_known=int(sp_src is not None), src_str_known=int(bool(src_key)),
               n_cs_ncbi=n_cs_ncbi, w_cs_ncbi=w_cs, n_cs_str=n_cs_str, w_cs_str=w_cs_str,
               n_mc=n_mc, w_mc=w_mc, w_used=wt, rep_cs=rep_cs,
               nonrep_cs=int(f_cs.sum()), nonrep_mc=int(f_mc.sum()), **hard)
    for v in VARS:
        rec["sb_" + v] = float(((Q[v] - yf) ** 2).sum())
        rec["sq_" + v] = float(Q[v].sum())
    cells = (np.stack([Q[v] for v in VARS]), y, cnt.astype(np.int16), n)
    return rec, cells


def _init():
    G["prov"] = Prov()


def auc_midrank(q, y):
    """Exact AUROC (Mann-Whitney with midranks for ties)."""
    u, inv = np.unique(q, return_inverse=True)
    c = np.bincount(inv, minlength=len(u)).astype(np.float64)
    p = np.bincount(inv, weights=y.astype(np.float64), minlength=len(u))
    below = np.cumsum(c) - c
    mid = below + (c + 1.0) / 2.0
    R = float((p * mid).sum())
    P = float(y.sum())
    Nn = float(len(y)) - P
    return (R - P * (P + 1) / 2.0) / (P * Nn)


def logit_slope(q, y, eps, grid=10000):
    """Grouped logistic regression y ~ a + s*logit(clip(q, eps, 1-eps)), q rounded to 1/grid."""
    b = np.rint(q * grid).astype(np.int64)
    nt = np.bincount(b, minlength=grid + 1).astype(np.float64)
    npos = np.bincount(b, weights=y.astype(np.float64), minlength=grid + 1)
    keep = nt > 0
    x = np.clip(np.arange(grid + 1) / grid, eps, 1 - eps)
    x = np.log(x / (1 - x))[keep]
    nt, npos = nt[keep], npos[keep]
    if np.var(x) == 0 or len(x) < 2:
        return float("nan")
    a, s = math.log(npos.sum() / (nt.sum() - npos.sum())), 0.0
    for _ in range(100):
        p = 1 / (1 + np.exp(-(a + s * x)))
        r = npos - nt * p
        v = nt * p * (1 - p)
        ga, gs = r.sum(), (r * x).sum()
        haa, has_, hss = v.sum(), (v * x).sum(), (v * x * x).sum()
        det = haa * hss - has_ ** 2
        da, ds = (hss * ga - has_ * gs) / det, (haa * gs - has_ * ga) / det
        a, s = a + da, s + ds
        if abs(da) + abs(ds) < 1e-12:
            break
    return s


def main():
    t0 = time.time()
    # ---- target keys + gate per-amplicon values
    gate = {}
    with gzip.open(CSV, "rt") as fh:
        for r in csv.DictReader(fh):
            if r["mode"] == "exclude_0.987" and r["scored"] == "1" and int(r["n"]) >= 2:
                gate[r["key"]] = r
    keys = set(gate)
    if LIMIT:   # smoke: the first LIMIT target keys in hits-file order (cheap prescan)
        first = []
        with open(HITS, "rb") as fh:
            for k in ijson.kvitems(fh, "", use_float=True):
                if k[0] in keys:
                    first.append(k[0])
                    if len(first) >= LIMIT:
                        break
        keys = set(first)
    log("target keys: %d" % len(keys))
    cat = json.load(open(CATF))
    G["cat_g"] = {g: v for g, v in cat["genomes"].items()}
    G["cat_t"] = {int(t): v for t, v in cat.get("taxa", {}).items()}
    idx = {"family": {}, "order": {}}
    for g, (sp, fam, od) in G["cat_g"].items():
        if sp is None:
            continue
        for r, lab in (("family", fam), ("order", od)):
            if lab is not None:
                idx[r].setdefault(lab, {}).setdefault(sp, []).append(g)
    for r in idx:
        for lab, d in list(idx[r].items()):
            idx[r][lab] = [(sp, sorted(d[sp])) for sp in sorted(d, key=lambda s: hashlib.sha1(str(s).encode()).hexdigest())]
    G["idx"] = idx
    log("catalogue: %d genomes, %d taxa" % (len(G["cat_g"]), len(G["cat_t"])))
    # ---- PGFam store (stream, read-only)
    vocab, off, chunks, pos = {}, {}, [], 0
    with open(CACHE, "rb") as fh:
        for gid, fams in ijson.kvitems(fh, "", use_float=True):
            if fams:
                a = np.unique(np.fromiter((vocab.setdefault(f, len(vocab)) for f in fams), np.int32, len(fams)))
            else:
                a = np.zeros(0, np.int32)
            off[gid] = (pos, pos + len(a))
            pos += len(a)
            chunks.append(a)
    G["data"] = np.concatenate(chunks)
    del chunks, vocab
    G["off"] = off
    log("store: %d genomes, %d ids (%.0fs)" % (len(off), len(G["data"]), time.time() - t0))
    total = sum(int(gate[k]["union_size"]) for k in keys)
    QA = np.zeros((len(VARS), total), np.float64)
    YA = np.zeros(total, bool)
    CA = np.zeros(total, np.int16)
    NA = np.zeros(total, np.int16)
    AA = np.zeros(total, np.int32)
    recs, ptr = [], 0
    sem = threading.BoundedSemaphore(WORKERS * 16)

    def feed():
        got = 0
        with open(HITS, "rb") as fh:
            for k, v in ijson.kvitems(fh, "", use_float=True):
                if k in keys:
                    sem.acquire()
                    yield (k, v)
                    got += 1
                    if got == len(keys):
                        return

    with mp.get_context("fork").Pool(WORKERS, initializer=_init) as pool:
        for i, (rec, cells) in enumerate(pool.imap_unordered(work, feed(), chunksize=2), 1):
            sem.release()
            if cells is not None:
                q, y, c, n = cells
                m = len(y)
                QA[:, ptr:ptr + m] = q
                YA[ptr:ptr + m] = y
                CA[ptr:ptr + m] = c
                NA[ptr:ptr + m] = n
                AA[ptr:ptr + m] = len(recs)
                rec["slot"] = len(recs)
                ptr += m
            recs.append(rec)
            if i % 2000 == 0:
                log("%d amplicons, %d cells" % (i, ptr))
    log("loop done: %d records, %d cells (expected %d)" % (len(recs), ptr, total))
    QA, YA, CA, NA, AA = QA[:, :ptr], YA[:ptr], CA[:ptr], NA[:ptr], AA[:ptr]
    res = {"n_records": len(recs), "n_cells": int(ptr), "n_cells_expected_from_csv": total}
    res["skipped"] = [r for r in recs if "skip" in r][:20]
    res["n_skipped"] = sum(1 for r in recs if "skip" in r)
    res["filter_equal_to_exclude_self_record"] = sum(r["same_filter"] for r in recs)
    ok = [r for r in recs if "skip" not in r]
    # ---- per-amplicon agreement with the gate's CSV
    agree = {"n_mismatch": 0, "union_size_mismatch": 0, "pos_mismatch": 0, "clade_rank_mismatch": 0,
             "n_members_mismatch": 0, "n_nb_mismatch": 0, "examples": []}
    maxrel = {"kn": 0.0, "e3s": 0.0}
    for r in ok:
        gr = gate[r["key"]]
        if int(gr["n"]) != r["n"]:
            agree["n_mismatch"] += 1
        if int(gr["union_size"]) != r["N"]:
            agree["union_size_mismatch"] += 1
        if int(gr["tp"]) != r["pos"]:
            agree["pos_mismatch"] += 1
        if (gr["clade_rank"] or "") != r["clade_rank"]:
            agree["clade_rank_mismatch"] += 1
        if int(gr["n_members"]) != r["n_members"]:
            agree["n_members_mismatch"] += 1
        if int(gr["n_nb_sp"]) != r["n_nb"]:
            agree["n_nb_mismatch"] += 1
        for v, c in (("kn", "brier_kn"), ("e3s", "brier_e3s")):
            mine = r["sb_" + v] / r["N"]
            theirs = float(gr[c])
            rel = abs(mine - theirs) / max(abs(theirs), 1e-12)
            if rel > maxrel[v]:
                maxrel[v] = rel
            if rel > 1e-5 and len(agree["examples"]) < 10:
                agree["examples"].append({"key": r["key"], "v": v, "mine": mine, "gate": theirs})
    agree["max_rel_diff_per_amplicon_brier"] = maxrel
    res["agreement_with_gate_csv"] = agree
    # ---- hard leakage assertions
    res["hard_leakage_counts"] = {k: int(sum(r[k] for r in ok)) for k in
                                  ("nb_is_source", "nb_src_prefix_taxid", "nb_identity_ge_tau",
                                   "member_src_species", "member_is_source")}
    # ---- soft leakage prevalence
    wt = sum(r["w_used"] for r in ok)
    res["conspecific"] = {
        "source_species_known": sum(r["sp_src_known"] for r in ok),
        "source_bvbrc_species_string_known": sum(r["src_str_known"] for r in ok),
        "amplicons_with_ncbi_conspecific_neighbour_incl_reps": sum(1 for r in ok if r["n_cs_ncbi"]),
        "amplicons_with_ncbi_conspecific_rep": sum(1 for r in ok if r["rep_cs"]),
        "amplicons_with_ncbi_conspecific_nonrep_neighbour": sum(1 for r in ok if r["nonrep_cs"]),
        "pooled_weight_share_ncbi_conspecific": sum(r["w_cs_ncbi"] for r in ok) / wt,
        "amplicons_with_bvbrc_string_conspecific_neighbour": sum(1 for r in ok if r["n_cs_str"]),
        "pooled_weight_share_bvbrc_string_conspecific": sum(r["w_cs_str"] for r in ok) / wt,
        "amplicons_with_removed_genome_reappearing": sum(1 for r in ok if r["n_mc"]),
        "amplicons_with_removed_genome_reappearing_nonrep": sum(1 for r in ok if r["nonrep_mc"]),
        "pooled_weight_share_removed_genome_reappearing": sum(r["w_mc"] for r in ok) / wt,
    }
    log("leakage tallies done")

    # ---- metrics
    def block(mask_cells, label):
        y = YA[mask_cells]
        out = {"cells": int(mask_cells.sum()) if mask_cells.dtype == bool else int(len(y)),
               "prevalence": float(y.mean())}
        for vi, v in enumerate(VARS):
            q = QA[vi][mask_cells]
            out[v] = {"brier": float(((q - y) ** 2).mean()), "auroc": auc_midrank(q, y), "mean_q": float(q.mean())}
        for v in VARS[1:]:
            out["d_" + v] = {"brier": out[v]["brier"] - out["kn"]["brier"], "auroc": out[v]["auroc"] - out["kn"]["auroc"]}
        log("metrics %s" % label)
        return out
    allm = np.ones(ptr, bool)
    res["n2plus"] = block(allm, "all n>=2")
    for lab, m in (("n=2", NA == 2), ("n=3", NA == 3), ("n>=4", NA >= 4)):
        res["by_n_" + lab] = block(m, lab)
    # subsets by amplicon
    def amp_mask(pred):
        good = np.zeros(len(recs), bool)
        for r in ok:
            if pred(r):
                good[r["slot"]] = True
        return good[AA]
    res["subset_no_ncbi_conspecific_incl_reps"] = block(amp_mask(lambda r: r["n_cs_ncbi"] == 0), "no conspecific")
    res["subset_no_conspecific_any_no_reappear"] = block(
        amp_mask(lambda r: r["n_cs_ncbi"] == 0 and r["n_cs_str"] == 0 and r["n_mc"] == 0), "clean")
    res["subset_with_ncbi_conspecific"] = block(amp_mask(lambda r: r["n_cs_ncbi"] > 0), "with conspecific")
    res["subset_amplicons"] = {
        "no_ncbi_conspecific_incl_reps": sum(1 for r in ok if r["n_cs_ncbi"] == 0),
        "clean": sum(1 for r in ok if r["n_cs_ncbi"] == 0 and r["n_cs_str"] == 0 and r["n_mc"] == 0),
        "with_ncbi_conspecific": sum(1 for r in ok if r["n_cs_ncbi"] > 0)}
    # saturated recalibration of k/n: per-(k, n) empirical prevalence (in-sample, generous)
    code = CA.astype(np.int64) * 10000 + NA.astype(np.int64)
    u, inv = np.unique(code, return_inverse=True)
    pg = np.bincount(inv, weights=YA.astype(np.float64)) / np.bincount(inv)
    qsat = pg[inv]
    res["kn_saturated_recalibration"] = {
        "n_levels": int(len(u)), "brier": float(((qsat - YA) ** 2).mean()), "auroc": auc_midrank(qsat, YA),
        "e3s_minus_knsat_brier": res["n2plus"]["e3s"]["brier"] - float(((qsat - YA) ** 2).mean()),
        "e3s_minus_knsat_auroc": res["n2plus"]["e3s"]["auroc"] - auc_midrank(qsat, YA)}
    log("saturated k/n done")
    res["logistic_slope"] = {v: {"eps1e-3": logit_slope(QA[vi], YA, 1e-3), "eps1e-2": logit_slope(QA[vi], YA, 1e-2)}
                             for vi, v in enumerate(VARS[:2])}
    # cluster bootstrap of dBrier on the leakage subsets (source = cluster; one amplicon per
    # source per region, sources unique per species in this panel)
    def boot_dbrier(pred, v, B=1000, seed=7):
        rows = [r for r in ok if pred(r)]
        srcs = sorted({r["src"] for r in rows})
        si = {s: i for i, s in enumerate(srcs)}
        Kc = len(srcs)
        Nc, Sk, Se = np.zeros(Kc), np.zeros(Kc), np.zeros(Kc)
        for r in rows:
            i = si[r["src"]]
            Nc[i] += r["N"]
            Sk[i] += r["sb_kn"]
            Se[i] += r["sb_" + v]
        rng = np.random.default_rng(seed)
        Wb = np.stack([np.bincount(rng.integers(0, Kc, Kc), minlength=Kc) for _ in range(B)]).astype(float)
        d = (Wb @ Se - Wb @ Sk) / (Wb @ Nc)
        return {"amplicons": len(rows), "sources": Kc, "dbrier": float((Se.sum() - Sk.sum()) / Nc.sum()),
                "ci95": [float(x) for x in np.percentile(d, [2.5, 97.5])]}
    res["boot_dbrier"] = {
        "all_e3s": boot_dbrier(lambda r: True, "e3s"),
        "all_e3s_clean": boot_dbrier(lambda r: True, "e3s_clean"),
        "no_ncbi_conspecific_e3s": boot_dbrier(lambda r: r["n_cs_ncbi"] == 0, "e3s"),
        "clean_subset_e3s": boot_dbrier(lambda r: r["n_cs_ncbi"] == 0 and r["n_cs_str"] == 0 and r["n_mc"] == 0, "e3s")}
    res["elapsed_s"] = round(time.time() - t0)
    with open(os.path.join(OUT, "indep_recompute%s.json" % ("_smoke" if LIMIT else "")), "w") as fh:
        json.dump(res, fh, indent=1, default=float)
    with gzip.open(os.path.join(OUT, "indep_per_amplicon%s.json.gz" % ("_smoke" if LIMIT else "")), "wt") as fh:
        json.dump(recs, fh, default=float)
    log("DONE")
    print(json.dumps(res, indent=1, default=float))


if __name__ == "__main__":
    main()
