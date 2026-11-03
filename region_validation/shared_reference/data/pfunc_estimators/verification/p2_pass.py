"""Phase 2: independent per-cell estimator pass, nested selection, composites, E4, leakage.

Everything is recomputed from the raw inputs (vctx + the phase-1 stream cache).
Outputs: cache/passA.npz (per cell x variant statistics), cache/passB.pkl (union-support
points of the composites + E4 maps), p2_summary.json (selection, reproduction, leakage).
"""
import collections
import json
import math
import multiprocessing as mp
import os
import pickle
import sys
import time

import numpy as np

import vctx

EPS = vctx.UNIVERSE_EPS
TAU_GRID1 = (0.0025, 0.005, 0.01, 0.02, 0.05)
K_GRID1 = (20, 50, 100)
TAU_X1 = (0.1, 0.2, math.inf)               # out-of-grid probes (exploratory)
ALPHA2 = (0.5, 1.0, 2.0, 4.0, 8.0)
TAU_GRID3 = (0.005, 0.01, 0.02)
ALPHA3 = (0.5, 1.0, 2.0, 4.0)
TAU_X3 = (0.05, 0.1, math.inf)              # out-of-grid probes (exploratory)
ALPHA_X3 = (1.0, 2.0, 4.0)
STARTS = ("genus", "family")
CHAIN = ("genus", "family", "order")
MINC = 3
NBINS_ISO = 10000                           # my binning (author used 4,000)

G = {}


def variants():
    V = []

    def add(**kw):
        kw.setdefault("tau", None)
        kw.setdefault("K", None)
        kw.setdefault("start", None)
        kw.setdefault("alpha", None)
        kw.setdefault("grid", False)
        V.append(kw)
    add(name="E0", fam="E0")
    for tau in TAU_GRID1:
        for K in K_GRID1:
            add(name=f"E1|K={K}|tau={tau}", fam="E1", tau=tau, K=K, grid=True)
    for tau in TAU_X1:
        for K in (50, 100):
            add(name=f"E1x|K={K}|tau={tau}", fam="E1x", tau=tau, K=K)
    for s in STARTS:
        for a in ALPHA2:
            add(name=f"E2|{s}|alpha={a}", fam="E2", start=s, alpha=a, grid=True)
    for tau in TAU_GRID3:
        for s in STARTS:
            for a in ALPHA3:
                add(name=f"E3|tau={tau}|{s}|alpha={a}", fam="E3", tau=tau, K=100, start=s,
                    alpha=a, grid=True)
    for tau in TAU_X3:
        for a in ALPHA_X3:
            add(name=f"E3x|tau={tau}|family|alpha={a}", fam="E3x", tau=tau, K=100,
                start="family", alpha=a)
    for K in (5, 10, 20, 50, 100):
        add(name=f"R200|K={K}|tau=0.01", fam="R200", tau=0.01, K=K)
    add(name="PREV_foldref", fam="BASE")                     # fold-reference prevalence
    add(name="CLADE_family_only", fam="BASE", start="family")  # pi alone (0 if no clade)
    return V


# --------------------------------------------------------------------------- setup
def setup():
    t0 = time.time()
    C = vctx.Ctx(load_truth=True)
    cells = vctx.load_cells("1.0")
    NG, NU, E = C.NG, C.NU, C.E
    A = collections.Counter()
    streams = {}
    for f in range(5):
        with open(os.path.join(vctx.CACHE, f"stream_f{f}.pkl"), "rb") as fh:
            streams[f] = pickle.load(fh)
    # per-fold reference EC matrices from _cache_ec_subset.json (the design's source)
    Ef, inref = {}, {}
    for f in range(5):
        sub = vctx.jl(os.path.join(vctx.DATA, "prfba", f"f{f:02d}_d100",
                                   "_cache_ec_subset.json"))
        M = np.zeros((NG, NU), np.uint8)
        m = np.zeros(NG, bool)
        for g, ecl in sub.items():
            i = C.gi[g]
            m[i] = True
            idx = [C.ei[n] for n in (vctx.norm_ec(x) for x in ecl) if n in C.ei]
            M[i, idx] = 1
            if not np.array_equal(M[i], E[i]):
                A["ec_subset_row_ne_panel_truth"] += 1
        Ef[f], inref[f] = M, m
    # per-cell joins
    nC = len(cells)
    src = np.array([C.gi[s] for s in cells.src_genome_id], np.int64)
    fold = cells.fold.to_numpy(np.int64)
    reps, l100, l200, s100same = [], [], [], []
    rec_n_surv = np.zeros(nC, np.int64)
    for i, (f, m) in enumerate(zip(fold, cells.md5)):
        st = streams[int(f)]
        pr = st["reps"][m]
        rc = st["recs"][m]
        r = np.array([C.gi[x] for x in pr["reps"]], np.int64)
        reps.append(r)
        g1 = np.array([C.gi[x[0]] for x in rc["L100"]], np.int64)
        i1 = np.array([x[1] for x in rc["L100"]], np.float64)
        g2 = np.array([C.gi[x[0]] for x in rc["L200"]], np.int64)
        i2 = np.array([x[1] for x in rc["L200"]], np.float64)
        l100.append((g1, i1))
        l200.append((g2, i2))
        rec_n_surv[i] = len(g1)
        # species dedup (hit-lineage species key) of the 100-hit list
        seen, keep = set(), []
        for x in rc["L100"]:
            if x[3] in seen:
                continue
            seen.add(x[3])
            keep.append(x[0])
        s100same.append(keep == [x[0] for x in rc["L100"]])
        if len(r) != int(cells.n_reps.iat[i]) or len(r) != pr["n_reps"]:
            A["n_reps_mismatch"] += 1
        # tie sets vs reps
        if len(g1) == 0:
            A["no_survivor_in_first_100_hits"] += 1
        else:
            if set(g1[i1 == i1[0]].tolist()) != set(r.tolist()):
                A["tie_set_100_ne_reps"] += 1
        if len(g2) and set(g2[i2 == i2[0]].tolist()) != set(r.tolist()):
            A["tie_set_200_ne_reps"] += 1
        if r[0] != (g1[0] if len(g1) else -1):
            A["reps0_ne_first_survivor_100"] += 1   # informative only (anchor order)
        # union EC == predictions' EC set; union F1 == cells f1
        u = E[r].max(0).astype(bool)
        pe = set(n for n in (vctx.norm_ec(x) for x in pr["ec"]) if n in C.ei)
        if set(np.flatnonzero(u).tolist()) != {C.ei[x] for x in pe}:
            A["union_ne_pred_ec"] += 1
        y = E[src[i]].astype(bool)
        f1 = 2 * (u & y).sum() / (u.sum() + y.sum())
        if abs(f1 - float(cells.f1.iat[i])) > 1e-3:
            A["union_f1_ne_cells_f1"] += 1
        if not inref[int(f)][r].all():
            A["reps_outside_fold_reference"] += 1
    s100same = np.array(s100same)
    # records (fold, md5) level count of species-dedup differences
    recs_seen = {}
    for i, (f, m) in enumerate(zip(fold, cells.md5)):
        recs_seen[(int(f), m)] = s100same[i]
    A["records"] = len(recs_seen)
    A["records_species_dedup_differs"] = int(sum(1 for v in recs_seen.values() if not v))

    # ---- clades (anchor = reps[0]; panel lineage; fold reference members) -------------
    lab = {rk: C.lab[rk] for rk in CHAIN}
    cnt = {}
    for f in range(5):
        for rk in CHAIN:
            cnt[(f, rk)] = collections.Counter(lab[rk][j] for j in np.flatnonzero(inref[f]))
    clade_key = {}
    cl_of = {s: np.full(nC, -1, np.int64) for s in STARTS}
    cl_rank = {s: np.full(nC, -1, np.int64) for s in STARTS}
    for i in range(nC):
        f = int(fold[i])
        a = reps[i][0]
        for s in STARTS:
            for ri in range(CHAIN.index(s), len(CHAIN)):
                rk = CHAIN[ri]
                L = lab[rk][a]
                if not L or L.startswith("unclassified_"):
                    continue
                if cnt[(f, rk)][L] < MINC:
                    continue
                key = (f, rk, L)
                if key not in clade_key:
                    clade_key[key] = len(clade_key)
                cl_of[s][i] = clade_key[key]
                cl_rank[s][i] = ri
                break
    members = [None] * len(clade_key)
    PI = np.zeros((len(clade_key), NU))
    for (f, rk, L), j in clade_key.items():
        mem = np.array([k for k in np.flatnonzero(inref[f]) if lab[rk][k] == L], np.int64)
        members[j] = mem
        PI[j] = Ef[f][mem].mean(0)
    PREV = {f: Ef[f][inref[f]].mean(0) for f in range(5)}
    A["setup_seconds"] = round(time.time() - t0, 1)
    G.update(C=C, cells=cells, E=E, Ef=Ef, inref=inref, src=src, fold=fold, reps=reps,
             l100=l100, l200=l200, rec_n_surv=rec_n_surv, cl_of=cl_of, cl_rank=cl_rank,
             members=members, PI=PI, PREV=PREV, clade_key=clade_key, A=A, NU=NU)


# --------------------------------------------------------------------------- per cell
def cell_Q(i, V):
    E, NU = G["E"], G["NU"]
    f = int(G["fold"][i])
    Ef = G["Ef"][f]
    r = G["reps"][i]
    n = len(r)
    k = Ef[r].sum(0).astype(np.float64)
    y = E[G["src"][i]].astype(bool)
    u = k > 0
    g1, i1 = G["l100"][i]
    X1 = Ef[g1].astype(np.float64)
    g2, i2 = G["l200"][i]
    Q = np.zeros((len(V), NU))
    cache_num = {}
    pi = {}
    for s in STARTS:
        j = G["cl_of"][s][i]
        pi[s] = (G["PI"][j], True) if j >= 0 else (None, False)
    for vi, v in enumerate(V):
        fam = v["fam"]
        if fam == "E0":
            Q[vi] = k / n
        elif fam in ("E1", "E1x", "E3", "E3x"):
            key = (v["tau"], v["K"])
            if key not in cache_num:
                d = i1[0] - i1[:v["K"]]
                w = np.ones_like(d) if math.isinf(v["tau"]) else np.exp(-d / v["tau"])
                cache_num[key] = (w @ X1[:v["K"]], w.sum())
            num, den = cache_num[key]
            if fam in ("E1", "E1x"):
                Q[vi] = num / den
            else:
                p, ok = pi[v["start"]]
                Q[vi] = (num + v["alpha"] * p) / (den + v["alpha"]) if ok else num / den
        elif fam == "E2":
            p, ok = pi[v["start"]]
            Q[vi] = (k + v["alpha"] * p) / (n + v["alpha"]) if ok else k / n
        elif fam == "R200":
            kk = min(v["K"], len(g2), 100)
            w = np.exp(-(i2[0] - i2[:kk]) / v["tau"])
            Q[vi] = (w @ G["Ef"][f][g2[:kk]].astype(np.float64)) / w.sum()
        elif v["name"] == "PREV_foldref":
            Q[vi] = G["PREV"][f]
        elif v["name"] == "CLADE_family_only":
            p, ok = pi["family"]
            Q[vi] = p if ok else 0.0
    Q[Q >= 1 - 1e-9] = 1.0
    Q[Q <= 1e-12] = 0.0
    return Q, y, u


STATS = ("SQ", "LL", "SQs", "Sqs", "Sqqs", "Sqys", "TP", "PRED", "F1", "NSUP", "TPu",
         "PREDu", "F1u")


def stats(Q, y, u):
    yf = y.astype(np.float64)
    T = y.sum()
    D = Q - yf
    SQ = (D * D).sum(1)
    Qc = np.clip(Q, EPS, 1 - EPS)
    LL = -(np.log(Qc[:, y]).sum(1) + np.log(1 - Qc[:, ~y]).sum(1))
    Qs = Q[:, u]
    ys = y[u]
    SQs = ((Qs - ys) ** 2).sum(1)
    Sqs = Qs.sum(1)
    Sqqs = (Qs * Qs).sum(1)
    Sqys = Qs[:, ys].sum(1)
    P = Q >= 0.5
    TP = P[:, y].sum(1)
    PRED = P.sum(1)
    F1 = np.where(PRED + T > 0, 2 * TP / np.maximum(PRED + T, 1), 0.0)
    Pu = P & u
    TPu = Pu[:, y].sum(1)
    PREDu = Pu.sum(1)
    F1u = np.where(PREDu + T > 0, 2 * TPu / np.maximum(PREDu + T, 1), 0.0)
    NSUP = (Q > 0).sum(1)
    return np.stack([SQ, LL, SQs, Sqs, Sqqs, Sqys, TP, PRED, F1, NSUP, TPu, PREDu, F1u], 1)


def passA_job(rng):
    lo, hi = rng
    V = G["V"]
    out = np.zeros((hi - lo, len(V), len(STATS)))
    for j, i in enumerate(range(lo, hi)):
        Q, y, u = cell_Q(i, V)
        out[j] = stats(Q, y, u)
    return lo, hi, out


def passB_job(rng):
    """Composite rows (per fold row choice), with optional isotonic maps: stats, union
    points, and per-fold universe histograms (for OOF isotonic fitting)."""
    lo, hi = rng
    V, comps = G["V"], G["COMPS"]
    nb = NBINS_ISO + 2
    st = np.zeros((hi - lo, len(comps), len(STATS)))
    H = np.zeros((len(comps), 3, nb))       # points, positives, sum q (this chunk = 1 fold)
    uq = [[] for _ in comps]
    uy, uc = [], []
    for j, i in enumerate(range(lo, hi)):
        f = int(G["fold"][i])
        rows = sorted({c["rows"][f] for c in comps})
        Vsub = [V[r] for r in rows]
        Qs, y, u = cell_Q(i, Vsub)
        pos = {r: t for t, r in enumerate(rows)}
        Qc = np.zeros((len(comps), G["NU"]))
        for c, comp in enumerate(comps):
            q = Qs[pos[comp["rows"][f]]].copy()
            if comp.get("maps") is not None:
                x_, y_ = comp["maps"][f]
                q = np.interp(q, x_, y_)
                q[q >= 1 - 1e-9] = 1.0
                q[q <= 1e-12] = 0.0
            Qc[c] = q
            b = np.where(q <= 0, 0, np.where(q >= 1, nb - 1,
                                             np.minimum((q * NBINS_ISO).astype(np.int64),
                                                        NBINS_ISO - 1) + 1))
            H[c, 0] += np.bincount(b, minlength=nb)
            H[c, 1] += np.bincount(b[y], minlength=nb)
            H[c, 2] += np.bincount(b, weights=q, minlength=nb)
            uq[c].append(q[u].astype(np.float32))
        st[j] = stats(Qc, y, u)
        uy.append(y[u])
        uc.append(np.full(int(u.sum()), i, np.int32))
    return lo, hi, st, H, [np.concatenate(x) for x in uq], np.concatenate(uy), np.concatenate(uc)


def chunks(size=400):
    fold = G["fold"]
    out = []
    for f in range(5):
        idx = np.flatnonzero(fold == f)
        lo, hi = int(idx[0]), int(idx[-1]) + 1
        assert hi - lo == len(idx)
        for a in range(lo, hi, size):
            out.append((a, min(a + size, hi)))
    return out


def run(fn, jobs, workers=10):
    res = []
    with mp.get_context("fork").Pool(workers) as pool:
        for r in pool.imap_unordered(fn, jobs):
            res.append(r)
    return sorted(res, key=lambda r: r[0])


def fit_iso(Hs):
    from sklearn.isotonic import IsotonicRegression
    M, P, S = Hs
    k = M > 0
    x = S[k] / M[k]
    yv = P[k] / M[k]
    ir = IsotonicRegression(y_min=0.0, y_max=1.0, increasing=True, out_of_bounds="clip")
    ir.fit(x, yv, sample_weight=M[k])
    return np.asarray(ir.X_thresholds_, float), np.asarray(ir.y_thresholds_, float)


# --------------------------------------------------------------------------- main
def main():
    T0 = time.time()
    V = variants()
    G["V"] = V
    setup()
    A = G["A"]
    print(f"[verify] setup done {time.time()-T0:.0f}s; asserts {dict(A)}", flush=True)
    jobs = chunks()
    res = run(passA_job, jobs)
    nC = len(G["cells"])
    ST = np.zeros((nC, len(V), len(STATS)))
    for lo, hi, o in res:
        ST[lo:hi] = o
    print(f"[verify] pass A done {time.time()-T0:.0f}s", flush=True)
    np.savez_compressed(os.path.join(vctx.CACHE, "passA.npz"), ST=ST)

    # ---------------- nested selection (declared grids only) ------------------------
    fold = G["fold"]
    NU = G["NU"]
    LLf = np.stack([ST[fold == f][:, :, 1].sum(0) for f in range(5)], 1)   # V x 5
    NCf = np.array([(fold == f).sum() for f in range(5)], float)

    def tiekey(v):
        return ((v["K"] or 0), -(v["tau"] or 0.0) if v["tau"] is not None else 0.0,
                0 if v["start"] in (None, "genus") else 1, v["alpha"] or 0.0)
    sel = {}
    for fam in ("E1", "E2", "E3"):
        rows = [j for j, v in enumerate(V) if v["fam"] == fam and v["grid"]]
        per = {}
        for f in range(5):
            tr = [x for x in range(5) if x != f]
            ll = LLf[rows][:, tr].sum(1) / (NCf[tr].sum() * NU)
            best = ll.min()
            cand = [rows[j] for j in range(len(rows)) if ll[j] <= best * (1 + 1e-12)]
            ch = min(cand, key=lambda r: tiekey(V[r]))
            srt = np.argsort(ll)
            per[f] = {"row": int(ch), "variant": V[ch]["name"],
                      "train_ll": float(ll[rows.index(ch)]),
                      "outer_ll": float(LLf[ch, f] / (NCf[f] * NU)),
                      "runner_up": V[rows[srt[1]]]["name"],
                      "runner_up_train_ll": float(ll[srt[1]]),
                      "margin_to_runner_up": float(ll[srt[1]] - ll[srt[0]])}
        sel[fam] = per
    # train-fold optimum of E3 over the extended tau probes too (exploratory)
    ext = {}
    for f in range(5):
        tr = [x for x in range(5) if x != f]
        rows = [j for j, v in enumerate(V) if v["fam"] in ("E3", "E3x")]
        ll = LLf[rows][:, tr].sum(1) / (NCf[tr].sum() * NU)
        ext[f] = V[rows[int(np.argmin(ll))]]["name"]
    print("[verify] selection:", {fam: {f: p["variant"] for f, p in d.items()}
                                  for fam, d in sel.items()}, flush=True)

    # ---------------- composites + pass B -------------------------------------------
    ix = {v["name"]: j for j, v in enumerate(V)}
    comps = [{"name": "E0", "rows": {f: 0 for f in range(5)}}]
    for fam in ("E1", "E2", "E3"):
        comps.append({"name": fam, "rows": {f: sel[fam][f]["row"] for f in range(5)}})
    fixed = ["E1|K=100|tau=0.02", "E3x|tau=0.05|family|alpha=2.0", "E3x|tau=inf|family|alpha=2.0",
             "E1x|K=100|tau=inf", "PREV_foldref", "CLADE_family_only"]
    for nm in fixed:
        comps.append({"name": nm, "rows": {f: ix[nm] for f in range(5)}})
    G["COMPS"] = comps
    resB = run(passB_job, jobs)
    print(f"[verify] pass B done {time.time()-T0:.0f}s", flush=True)
    STB = np.zeros((nC, len(comps), len(STATS)))
    Hf = {f: np.zeros((len(comps), 3, NBINS_ISO + 2)) for f in range(5)}
    UQ, UY, UC = [[] for _ in comps], [], []
    for lo, hi, st, H, uq, uy, uc in resB:
        STB[lo:hi] = st
        Hf[int(fold[lo])] += H
        for c in range(len(comps)):
            UQ[c].append(uq[c])
        UY.append(uy)
        UC.append(uc)
    UQ = np.stack([np.concatenate(x) for x in UQ], 0)
    UY = np.concatenate(UY)
    UC = np.concatenate(UC)

    # ---------------- E4 (OOF isotonic) of E0 and E3, and the all-fold E3 map -----------
    maps = {}
    comps4 = []
    for base in ("E0", "E3"):
        c = [x["name"] for x in comps].index(base)
        mf = {}
        for f in range(5):
            Hs = sum(Hf[f2][c] for f2 in range(5) if f2 != f)
            mf[f] = fit_iso(Hs)
        maps[f"E4({base})"] = mf
        comps4.append({"name": f"E4({base})", "rows": comps[c]["rows"], "maps": mf})
    c3 = [x["name"] for x in comps].index("E3")
    iso_all = fit_iso(sum(Hf[f][c3] for f in range(5)))
    G["COMPS"] = comps4
    resC = run(passB_job, jobs)
    print(f"[verify] pass C done {time.time()-T0:.0f}s", flush=True)
    STC = np.zeros((nC, len(comps4), len(STATS)))
    UQ4 = [[] for _ in comps4]
    for lo, hi, st, H, uq, uy, uc in resC:
        STC[lo:hi] = st
        for c in range(len(comps4)):
            UQ4[c].append(uq[c])
    UQ4 = np.stack([np.concatenate(x) for x in UQ4], 0)
    allcomps = [c["name"] for c in comps] + [c["name"] for c in comps4]
    STALL = np.concatenate([STB, STC], 1)
    UQALL = np.concatenate([UQ, UQ4], 0)
    with open(os.path.join(vctx.CACHE, "passB.pkl"), "wb") as fh:
        pickle.dump({"comps": allcomps, "ST": STALL, "UQ": UQALL, "UY": UY, "UC": UC,
                     "maps": {k: {f: (a.tolist(), b.tolist()) for f, (a, b) in v.items()}
                              for k, v in maps.items()},
                     "iso_all_E3": (iso_all[0].tolist(), iso_all[1].tolist()),
                     "rows": {c["name"]: c["rows"] for c in comps + comps4}}, fh, protocol=5)

    # ---------------- leakage (every cell) -----------------------------------------
    C = G["C"]
    code = C.code
    spn = np.array([hash(vctx.sp_norm(s)) for s in C.lab["species"]])
    L = collections.Counter()
    for i in range(nC):
        s = G["src"][i]
        f = int(fold[i])
        a = G["reps"][i][0]
        L["anchor_genus_eq_source_genus"] += int(code["genus"][a] == code["genus"][s])
        lists = {"L100": G["l100"][i][0], "L200": G["l200"][i][0], "reps": G["reps"][i]}
        for st_ in STARTS:
            j = G["cl_of"][st_][i]
            if j >= 0:
                lists[f"clade_{st_}"] = G["members"][j]
        for nm, g in lists.items():
            L[f"{nm}:source_member"] += int((g == s).any())
            for c_ in ("genus", "species", "taxon", "group"):
                L[f"{nm}:shares_{c_}"] += int((code[c_][g] == code[c_][s]).any())
            L[f"{nm}:shares_species_norm"] += int((spn[g] == spn[s]).any())
            L[f"{nm}:outside_fold_reference"] += int((~G["inref"][f][g]).sum())
    # clade fold == cell fold
    inv = {v: k for k, v in G["clade_key"].items()}
    for st_ in STARTS:
        j = G["cl_of"][st_]
        L[f"clade_{st_}:fold_ne_cell_fold"] = int(sum(1 for i in range(nC) if j[i] >= 0 and
                                                      inv[j[i]][0] != int(fold[i])))
        L[f"clade_{st_}:rank_counts"] = dict(collections.Counter(
            ("none" if r < 0 else CHAIN[r]) for r in G["cl_rank"][st_]))
    # why none (family start)
    why = collections.Counter()
    for i in np.flatnonzero(G["cl_of"]["family"] < 0):
        a = G["reps"][i][0]
        o = C.lab["order"][a]
        why["order_placeholder" if (not o or o.startswith("unclassified_")) else
            "order_lt3"] += 1
    L["family_start_none_reason"] = dict(why)

    # ---------------- reproduction targets -------------------------------------------
    cs = vctx.jl(os.path.join(vctx.DATA, "calibration", "scores.json"))
    t = cs["exploratory_topk"]["all|all_records"]
    et = t["exact_tie_k_over_n"]
    targ = {"E0_universe_brier": et["brier_universe"]["est"],
            "E0_universe_logloss": et["logloss_universe_eps0.01"]["est"],
            "E0_union_macro_F1": et["union.macro_F1"]["est"],
            "E0_union_micro_F1": et["union.micro_F1"]["est"],
            "E0_union_support_brier": et["brier_union_support"]["est"],
            "E0_p_ge_0.5_macro_F1": et["p_ge_0.5.macro_F1"]["est"]}
    for K in (5, 10, 20, 50, 100):
        targ[f"R200_K{K}_universe_brier"] = t[f"k={K}"]["brier_universe"]["est"]
        targ[f"R200_K{K}_q_ge_0.5_macro_F1"] = t[f"k={K}"]["q_ge_0.5.macro_F1"]["est"]
    S = {nm: j for j, nm in enumerate(STATS)}
    e0 = ST[:, 0]
    # union quantities (estimator independent)
    Ns = np.array([int((G["Ef"][int(fold[i])][G["reps"][i]].max(0) > 0).sum())
                   for i in range(nC)], float)
    Tt = G["E"][G["src"]].sum(1).astype(float)
    POS = np.array([int(((G["Ef"][int(fold[i])][G["reps"][i]].max(0) > 0) &
                         G["E"][G["src"][i]].astype(bool)).sum()) for i in range(nC)], float)
    mine = {"E0_universe_brier": e0[:, S["SQ"]].sum() / (nC * NU),
            "E0_universe_logloss": e0[:, S["LL"]].sum() / (nC * NU),
            "E0_union_macro_F1": float(np.mean(2 * POS / (Ns + Tt))),
            "E0_union_micro_F1": float(2 * POS.sum() / (Ns.sum() + Tt.sum())),
            "E0_union_support_brier": e0[:, S["SQs"]].sum() / Ns.sum(),
            "E0_p_ge_0.5_macro_F1": e0[:, S["F1"]].mean()}
    for K in (5, 10, 20, 50, 100):
        b = ST[:, ix[f"R200|K={K}|tau=0.01"]]
        mine[f"R200_K{K}_universe_brier"] = b[:, S["SQ"]].sum() / (nC * NU)
        mine[f"R200_K{K}_q_ge_0.5_macro_F1"] = b[:, S["F1"]].mean()
    b = ST[:, ix["E1|K=20|tau=0.01"]]
    mine["E1_100hit_K20_universe_brier"] = b[:, S["SQ"]].sum() / (nC * NU)
    # E0 universe Brier decomposition: union SQ + truth ECs outside the union
    mine["E0_universe_SQ_minus_unionSQ_minus_FN_outside_union"] = float(
        e0[:, S["SQ"]].sum() - e0[:, S["SQs"]].sum() - (Tt - POS).sum())
    repro = {k: {"mine": float(v), "target": targ.get(k),
                 "abs_diff": (abs(float(v) - targ[k]) if k in targ else None)}
             for k, v in mine.items()}
    np.savez_compressed(os.path.join(vctx.CACHE, "union_quantities.npz"), Ns=Ns, T=Tt, POS=POS)

    # ---------------- compare with selected.json -------------------------------------
    so = vctx.jl(os.path.join(vctx.DATA, "pfunc_estimators", "selected.json"))
    cmp_sel = {}
    tr_their = so["train_universe_logloss_by_variant_and_outer_fold"]

    def their_name(v):
        if v["fam"] == "E1":
            return f"E1|K={v['K']}|tau={v['tau']}|genome"
        if v["fam"] == "E2":
            return f"E2|{v['start']}|alpha={v['alpha']:g}"
        if v["fam"] == "E3":
            return f"E3|K=100|tau={v['tau']}|genome|{v['start']}|alpha={v['alpha']:g}"
        return None
    maxdiff = 0.0
    for j, v in enumerate(V):
        tn = their_name(v)
        if tn is None or not v["grid"]:
            continue
        for f in range(5):
            tr = [x for x in range(5) if x != f]
            mine_ll = LLf[j, tr].sum() / (NCf[tr].sum() * NU)
            maxdiff = max(maxdiff, abs(mine_ll - tr_their[tn][str(f)]))
    cmp_sel["max_abs_diff_train_logloss_all_grid_variants_all_folds"] = maxdiff
    for fam in ("E1", "E2", "E3"):
        cmp_sel[fam] = {f: {"mine": sel[fam][f]["variant"],
                            "theirs": so["families"][fam]["per_fold"][str(f)]["variant"],
                            "mine_train_ll": sel[fam][f]["train_ll"],
                            "theirs_train_ll": so["families"][fam]["per_fold"][str(f)]["train_universe_logloss"],
                            "mine_outer_ll": sel[fam][f]["outer_ll"],
                            "theirs_outer_ll": so["families"][fam]["per_fold"][str(f)]["outer_universe_logloss"],
                            "runner_up": sel[fam][f]["runner_up"],
                            "margin": sel[fam][f]["margin_to_runner_up"]}
                        for f in range(5)}
    # their all-fold isotonic map vs mine
    sow = vctx.jl(os.path.join(vctx.DATA, "pfunc_estimators", "stage_o_weights.json"))
    iso_t = sow["meta"]["isotonic"]
    grid = np.linspace(0, 1, 2001)
    cmp_iso = {"their_knots": len(iso_t["x"]), "my_knots": len(iso_all[0]),
               "their_r0_r05_r1": [float(np.interp(z, iso_t["x"], iso_t["y"])) for z in (0, .5, 1)],
               "my_r0_r05_r1": [float(np.interp(z, *iso_all)) for z in (0, .5, 1)],
               "max_abs_diff_on_grid": float(np.abs(np.interp(grid, iso_t["x"], iso_t["y"]) -
                                                    np.interp(grid, *iso_all)).max()),
               "their_max_y": float(max(iso_t["y"])),
               "my_oof_r1_by_fold": {f: float(np.interp(1.0, *maps["E4(E3)"][f])) for f in range(5)}}
    summ = {"asserts_setup": dict(A), "leakage": dict(L), "reproduction": repro,
            "selection": {fam: {str(f): p for f, p in d.items()} for fam, d in sel.items()},
            "selection_extended_tau_probe_E3": ext,
            "selection_vs_selected_json": cmp_sel, "isotonic_vs_stage_o": cmp_iso,
            "variants": [v["name"] for v in V], "stats": list(STATS),
            "seconds": round(time.time() - T0, 1)}
    with open(os.path.join(vctx.OUT, "p2_summary.json"), "w") as fh:
        json.dump(summ, fh, indent=1, default=float)
    print(json.dumps({"reproduction": repro, "selection_vs": cmp_sel, "iso": cmp_iso,
                      "leak": dict(L), "asserts": dict(A)}, indent=1, default=float))


if __name__ == "__main__":
    main()
