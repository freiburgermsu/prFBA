#!/usr/bin/env python
"""Score the synthetic-MAG arm against the SAME truth and metric definitions the
16S amplicon arm uses, so the two markers are comparable organism-by-organism.

Comparability ("one vocabulary") rules carried over verbatim from
region_validation/scripts/score_taxacc.py and score_genecap.py:
  * Truth lineage = benchmark_truth.json's truth_lineage (NCBI via taxopy), which
    is None at unnamed ranks. A None truth rank is UNSCOREABLE: excluded from the
    denominator, never counted as wrong.
  * The PREDICTED lineage is re-resolved through the SAME _common.lineage_for()
    that produced the truth, from taxon_id = int(genome_id.split('.')[0]). Names
    are compared after _common._norm().
  * Per-rank cell is three-state: None = no-call (abstention / ambiguous
    prediction / truth undefined), True = correct, False = active misassignment.
  * `correct` = correct/called (CONDITIONAL on the pipeline answering).
    `correct_population` = correct x coverage over all scoreable queries, with the
    coverage denominator fixed to the FULL include-regime population so that
    exclusion abstentions cannot shrink it (the anti-inflation rule).
  * Gene capture: truth = the SOURCE genome's PGFam set; prediction = the union of
    the selected references' PGFam sets. Conditional macro = mean over queries
    that resolved; population macro additionally scores an abstention as
    recall/F1/Jaccard = 0 (precision undefined). truth-empty queries are excluded
    as undefined, not as misses.
  * Wilson 95% CI for rank proportions; normal-approx mean CI for means of
    proportions (gene capture) -- matching _macro() in score_genecap.py.

Regimes:
  include    -- hits as returned by skani.
  exclude    -- drop any hit that is the source genome, shares the source taxon,
                or sits at ANI >= TAU_ANI (95, the ANI species boundary and the
                structural analogue of the amplicon arm's 0.987 identity).
"""
import argparse, collections, csv, json, math, os, sys

ROOT = os.path.expanduser("~/Documents/prFBA")
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "region_validation", "scripts"))

import _common
from _common import lineage_for, _norm, wilson_ci
from mag_skani_hits import ANI_KNOBS, best_hit, classify
import select_references as SR

RANKS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
SCORED = ["Phylum", "Class", "Order", "Family", "Genus", "Species"]
# Exclusion stringency. 95 is the ANI species boundary and the default; --tau sweeps it.
# NOTE 99 is NOT a species boundary -- it retains same-species references at 95-99 ANI,
# so it is a strictly weaker exclusion and is not comparable to the amplicon arm's
# tau = 0.987 without saying so.
TAU_ANI = 95.0
AF_GATE = ANI_KNOBS["cov_min"] * 100.0     # 15.0

OUT_OK, OUT_ABSTAIN, OUT_TRUTH_EMPTY = "ok", "abstained", "truth_empty"


def src_taxon(gid):
    try:
        return int(str(gid).split(".")[0])
    except Exception:
        return None


def load_truth(p=None):
    """benchmark_truth.json -> src_genome_id -> truth_lineage.

    Parameterized: the k-means panel has its OWN truth file, and silently reusing
    the original panel's would score the new panel against the wrong organisms --
    it covers only the 804 genomes the two panels share.
    """
    p = p or os.path.join(ROOT, "region_validation", "data", "benchmark_truth.json")
    out = {}
    for v in json.load(open(p)).values():
        out.setdefault(v["src_genome_id"], v["truth_lineage"])
    return out


def load_pgfams(p=None):
    """Panel-only PGFam index: genome_id -> frozenset(int family ids).

    Built once by the pgfam_panel.json extraction step. The full BV-BRC cache is
    111,789 genomes and expands to ~25 GB of Python string objects; the reference
    database here IS the 10,000-genome panel, so only panel genomes can ever be
    selected or be the truth, and an int-encoded panel-only index is ~106 MB.
    """
    p = p or os.path.join(ROOT, "mag_validation", "pgfam_panel.json")
    d = json.load(open(p))["genomes"]
    return {k: frozenset(v) for k, v in d.items()}


def build_record(mag, hits, mag_bp):
    """skani hit rows -> the ASV-schema record select_representatives consumes."""
    hs = sorted(hits, key=lambda h: -h["ani"])
    top = []
    for i, h in enumerate(hs, 1):
        tid = src_taxon(h["ref"])
        lin = lineage_for(tid)
        top.append(dict(
            rank=i, identity=round(h["ani"] / 100.0, 4), align_score=h["ani"],
            aligned_len=int(round(h["af_query"] / 100.0 * mag_bp)),
            edlib_identity=None, edlib_distance=None, n_matches=None, md5=None,
            feature_id=None, taxon_id=tid, genome_id=h["ref"],
            organism=lin.get("Species") or lin.get("Genus") or h["ref"],
            lineage={k: lin.get(k) for k in RANKS},
            ani=h["ani"], af_ref=h["af_ref"], af_query=h["af_query"],
            af_pass=h["af_query"] >= AF_GATE))
    # classify() formats the novelty string with best["_dbword"], so carry it
    b = best_hit([dict(ani=h["ani"], af_query=h["af_query"], _dbword="bvbrc")
                  for h in hs], AF_GATE) if hs else None
    return dict(asv_len=mag_bp,
                best_identity=round(b["ani"] / 100.0, 4) if (b and b["af_pass"]) else None,
                midas_taxonomy=None, rel_ab=None, top20=top), b


def exclude_hits(hits, src_gid, src_tid):
    out = []
    for h in hits:
        if h["ref"] == src_gid:
            continue
        if src_tid is not None and src_taxon(h["ref"]) == src_tid:
            continue
        if h["ani"] >= TAU_ANI:
            continue
        out.append(h)
    return out


def rank_cell(truth_v, pred_v):
    if not pred_v:
        return None
    if not truth_v:
        return None
    return _norm(pred_v) == _norm(truth_v)


def score_query(mag, src_gid, hits, mag_bp, truth, pgf):
    """Return {regime: {cells, outcome, genecap, n_selected, ani_class}}."""
    src_tid = src_taxon(src_gid)
    tl = truth.get(src_gid) or {}
    truth_genes = pgf.get(src_gid) or set()
    res = {}
    for regime, hh in (("include", hits), ("exclude", exclude_hits(hits, src_gid, src_tid))):
        rec, b = build_record(mag, hh, mag_bp)
        ani_class, _ = classify(b or {}, len(hh), AF_GATE)
        try:
            sel = SR.select_representatives(rec, knobs=ANI_KNOBS,
                                            gene_provider=lambda g: set(pgf.get(g) or ()))
        except Exception:
            sel = {"selected": [], "tier": "error", "n_selected": 0}
        picks = [s["genome_id"] for s in (sel.get("selected") or [])]
        anchor = picks[0] if picks else None
        pred_lin = lineage_for(src_taxon(anchor)) if anchor else {k: None for k in RANKS}
        cells = {r: rank_cell(tl.get(r), pred_lin.get(r)) for r in SCORED}
        if not truth_genes:
            outcome, gc = OUT_TRUTH_EMPTY, None
        elif not picks:
            outcome, gc = OUT_ABSTAIN, dict(recall=0.0, precision=None, f1=0.0, jaccard=0.0)
        else:
            pred_genes = set()
            for g in picks:
                pred_genes |= (pgf.get(g) or frozenset())
            tp = len(truth_genes & pred_genes)
            rec_v = tp / len(truth_genes)
            pre_v = tp / len(pred_genes) if pred_genes else None
            f1 = (2 * rec_v * pre_v / (rec_v + pre_v)) if (pre_v and (rec_v + pre_v)) else 0.0
            jac = tp / len(truth_genes | pred_genes) if (truth_genes | pred_genes) else 0.0
            outcome, gc = OUT_OK, dict(recall=rec_v, precision=pre_v, f1=f1, jaccard=jac)
        res[regime] = dict(cells=cells, outcome=outcome, genecap=gc,
                           n_selected=len(picks), anchor=anchor,
                           tier=sel.get("tier"), ani_class=ani_class,
                           best_ani=(b or {}).get("ani"), n_hits=len(hh))
    return res


def _macro(vals):
    v = [x for x in vals if x is not None and not (isinstance(x, float) and math.isnan(x))]
    if not v:
        return (None, None, None, 0)
    m = sum(v) / len(v)
    if len(v) < 2:
        return (m, None, None, len(v))
    sd = math.sqrt(sum((x - m) ** 2 for x in v) / (len(v) - 1))
    h = 1.959964 * sd / math.sqrt(len(v))
    return (m, max(0.0, m - h), min(1.0, m + h), len(v))


def main():
    global TAU_ANI
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", action="append", required=True, help="skani search TSV (repeatable)")
    ap.add_argument("--manifest", action="append", required=True, help="synth manifest json (repeatable)")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--tau", type=float, default=TAU_ANI,
                    help="ANI exclusion threshold (default 95, the species boundary)")
    ap.add_argument("--truth", default=None,
                    help="benchmark_truth.json for THIS panel (default: the original panel's)")
    ap.add_argument("--pgfams", default=None,
                    help="compact PGFam index for THIS panel (default: the original panel's)")
    a = ap.parse_args()
    TAU_ANI = a.tau
    print("exclusion tau_ANI = %.1f" % TAU_ANI, flush=True)
    os.makedirs(a.outdir, exist_ok=True)

    print("loading truth + PGFams ...", flush=True)
    truth, pgf = load_truth(a.truth), load_pgfams(a.pgfams)
    print("  truth=%s" % (a.truth or "<original panel>"), flush=True)
    print("  pgfams=%s" % (a.pgfams or "<original panel>"), flush=True)
    meta = {}
    for mp in a.manifest:
        for r in json.load(open(mp))["rows"]:
            meta[r["mag"]] = r
    print("  truth=%d  pgfam=%d  synthetic MAGs=%d" % (len(truth), len(pgf), len(meta)), flush=True)

    hits = collections.defaultdict(list)
    nrow = 0
    for hp in a.hits:
        with open(hp) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                q = os.path.basename(row["Query_file"])[:-4]
                hits[q].append(dict(ref=os.path.basename(row["Ref_file"])[:-4],
                                    ani=float(row["ANI"]),
                                    af_query=float(row["Align_fraction_query"]),
                                    af_ref=float(row["Align_fraction_ref"])))
                nrow += 1
    print("  hit rows=%d over %d queries" % (nrow, len(hits)), flush=True)

    cells_w = csv.writer(open(os.path.join(a.outdir, "mag_cells.csv"), "w", newline=""))
    cells_w.writerow(["mag", "genome_id", "completeness", "contamination", "regime",
                      "outcome", "ani_class", "tier", "best_ani", "n_hits", "n_selected",
                      "anchor", "recall", "precision", "f1", "jaccard"] + SCORED)
    agg = collections.defaultdict(lambda: collections.defaultdict(list))
    n = 0
    for mag, m in sorted(meta.items()):
        gid = m["genome_id"]
        res = score_query(mag, gid, hits.get(mag, []), m["src_bp"] + m["contam_bp"], truth, pgf)
        for regime, r in res.items():
            key = (m["completeness_target"], m["contamination_target"], regime)
            gc = r["genecap"] or {}
            cells_w.writerow([mag, gid, m["completeness_target"], m["contamination_target"],
                              regime, r["outcome"], r["ani_class"], r["tier"], r["best_ani"],
                              r["n_hits"], r["n_selected"], r["anchor"] or "",
                              gc.get("recall"), gc.get("precision"), gc.get("f1"), gc.get("jaccard")]
                             + [("" if r["cells"][x] is None else int(r["cells"][x])) for x in SCORED])
            agg[key]["rows"].append((r, m))
        n += 1
        if n % 5000 == 0:
            print("  scored %d/%d" % (n, len(meta)), flush=True)

    summary = {}
    for key, d in sorted(agg.items()):
        comp, contam, regime = key
        rows = d["rows"]
        inc_n = len(rows)                      # same population in both regimes
        s = dict(completeness=comp, contamination=contam, regime=regime, n=inc_n)
        for rank in SCORED:
            scoreable = [r for r, m in rows if (truth.get(m["genome_id"]) or {}).get(rank)]
            called = [r for r in scoreable if r["cells"][rank] is not None]
            corr = [r for r in called if r["cells"][rank] is True]
            k, nn = len(corr), len(called)
            lo, hi = wilson_ci(k, nn)
            pk, pn = len(corr), len(scoreable)
            plo, phi = wilson_ci(pk, pn)
            s[rank] = dict(n_scoreable=pn, n_called=nn,
                           correct=(k / nn if nn else None), ci=[lo, hi],
                           coverage=(nn / pn if pn else None),
                           correct_population=(pk / pn if pn else None), ci_pop=[plo, phi])
        ok = [r for r, m in rows if r["outcome"] == OUT_OK]
        pop = [r for r, m in rows if r["outcome"] in (OUT_OK, OUT_ABSTAIN)]
        for lbl, pool in (("cond", ok), ("pop", pop)):
            for met in ("recall", "precision", "f1", "jaccard"):
                mm, lo, hi, nv = _macro([(r["genecap"] or {}).get(met) for r in pool])
                s["genecap_%s_%s" % (lbl, met)] = dict(mean=mm, ci=[lo, hi], n=nv)
        s["coverage_answered"] = len(ok) / len(pop) if pop else None
        s["outcomes"] = dict(collections.Counter(r["outcome"] for r, m in rows))
        s["ani_classes"] = dict(collections.Counter(r["ani_class"] for r, m in rows))
        s["median_n_selected"] = sorted(r["n_selected"] for r, m in rows)[len(rows) // 2]
        summary["c%g_x%g_%s" % (comp, contam, regime)] = s

    json.dump(dict(tau_ani=TAU_ANI, af_gate=AF_GATE, knobs=ANI_KNOBS,
                   n_queries=n, summary=summary),
              open(os.path.join(a.outdir, "mag_summary.json"), "w"), indent=1)
    print("wrote %s" % os.path.join(a.outdir, "mag_summary.json"), flush=True)


if __name__ == "__main__":
    main()
