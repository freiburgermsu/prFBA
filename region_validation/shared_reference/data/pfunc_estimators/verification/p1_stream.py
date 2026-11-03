"""Phase 1: stream each d100 fold's predictions.json (exact-tie reps) and hits file once.

Per record (md5) keeps the stage-0 / identity-ranked / genome-deduplicated neighbour lists
built from the FIRST 100 hits ("L100") and from all stored hits ("L200"), plus checks.
Stage 0 follows select_references.select_representatives: coverage = round(min(aligned_len /
asv_len, 1), 4) >= 0.90; |SW - edlib| = round(abs(identity - edlib_identity), 4) <= 0.02 (if
edlib_identity present); raw identity >= 0.865. Survivors are sorted by (-identity,
-align_score) (stable on rank order) and deduplicated by genome_id (first kept).
"""
import collections
import multiprocessing as mp
import os
import pickle
import sys
import time

import ijson

import vctx


def survivors(hits, L):
    keep = []
    for h in hits:
        cov = round(min(h["aligned_len"] / L, 1.0), 4)
        if cov < vctx.COV_MIN:
            continue
        ei = h.get("edlib_identity")
        if ei is not None and round(abs(h["identity"] - ei), 4) > vctx.SWED_MAX:
            continue
        if h["identity"] < vctx.FLOOR:
            continue
        keep.append(h)
    keep = sorted(keep, key=lambda h: (-h["identity"], -h["align_score"]))   # stable
    out, seen = [], set()
    for h in keep:
        g = str(h["genome_id"])
        if g in seen:
            continue
        seen.add(g)
        lin = h.get("lineage") or {}
        sp = vctx.sp_norm(lin.get("Species"))
        if not sp:
            sp = f"taxon:{h['taxon_id']}" if h.get("taxon_id") else g
        out.append((g, float(h["identity"]), float(h["align_score"]), sp))
    return out


def job(fold):
    t0 = time.time()
    cells = vctx.load_cells("1.0")
    want = set(cells[cells.fold == fold].md5)
    ddir = os.path.join(vctx.DATA, "prfba", f"f{fold:02d}_d100")
    ref = set(vctx.jl(os.path.join(ddir, "_cache_ec_subset.json")).keys())
    A = collections.Counter()
    reps = {}
    with open(os.path.join(ddir, "predictions.json"), "rb") as fh:
        for m, rec in ijson.kvitems(fh, "pred.exact_tie", use_float=True):
            if m in want:
                reps[m] = {"reps": [str(x) for x in rec["reps"]], "n_reps": int(rec["n_reps"]),
                           "ec": list(rec["ec"]), "abstained": rec.get("abstained")}
    recs = {}
    with open(os.path.join(ddir, "hits", "asv_top20_alignment_hits.json"), "rb") as fh:
        for m, rec in ijson.kvitems(fh, "", use_float=True):
            if m not in want:
                continue
            if not isinstance(rec, dict):
                A["record_not_dict"] += 1
                continue
            top = rec.get("top20") or []
            L = rec.get("asv_len")
            if not L:
                A["asv_len_missing_or_zero"] += 1
                L = 1
            if [h["rank"] for h in top] != list(range(1, len(top) + 1)):
                A["rank_order_violation"] += 1
            # align_score monotone non-increasing along rank?
            sc = [h["align_score"] for h in top]
            if any(sc[i] < sc[i + 1] for i in range(len(sc) - 1)):
                A["align_score_not_monotone_in_rank"] += 1
            A["hit_genomes_outside_fold_reference_all200"] += sum(
                1 for h in top if str(h["genome_id"]) not in ref)
            l100 = survivors(top[:100], L)
            l200 = survivors(top, L)
            # coverage-qualified hits among the first 100 (unrounded, as a diagnostic)
            q = [h for h in top[:100] if h["aligned_len"] / L >= vctx.COV_MIN]
            recs[m] = {"n_hits": len(top), "asv_len": L, "L100": l100, "L200": l200,
                       "last_cov_qualified_identity_100": (q[-1]["identity"] if q else None),
                       "n_cov_qualified_100": len(q)}
    A["records_wanted"] = len(want)
    A["records_found"] = len(recs)
    A["preds_found"] = len(reps)
    out = {"fold": fold, "reps": reps, "recs": recs, "asserts": dict(A),
           "seconds": round(time.time() - t0, 1)}
    with open(os.path.join(vctx.CACHE, f"stream_f{fold}.pkl"), "wb") as fh:
        pickle.dump(out, fh, protocol=5)
    return fold, dict(A), out["seconds"]


if __name__ == "__main__":
    folds = [int(x) for x in sys.argv[1:]] or [0, 1, 2, 3, 4]
    with mp.get_context("fork").Pool(min(5, len(folds))) as pool:
        for r in pool.imap_unordered(job, folds):
            print(r, flush=True)
