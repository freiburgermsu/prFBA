#!/usr/bin/env python
"""Score PICRUSt2-via-extracted-16S against the prFBA MAG arm, in the benchmark's EC space.

A MAG reaches PICRUSt2 only through a 16S gene it carries; one locus has to be chosen
per MAG, and the rule here is LONGEST (ties: higher HMM bit score), applied identically
in both parts.

Part A -- the real thing, on a sample: MAGs regenerated with ``synth_mags.py``, 16S
extracted from the MAG's own contigs with nhmmer (``mag_extract_16s.py``), stock
PICRUSt2 run on the extracts. Reports: whether the extract is byte-identical to a
locus of the source genome (so Part B's locus-level predictions are what a user would
get), agreement with PICRUSt2's prediction for the benchmark's FullLength16S amplicon
of the same genome, whether a MAG's 16S copies yield one prediction, and, for
10%-contaminated MAGs, how often the chosen locus is the donor's and what that costs.

Part B -- every uncontaminated synthetic MAG the MAG arm searched (9,997 genomes x 14
completeness tiers, replayed by ``mag_16s_retention.py``):
  * ``picrust2_mag``       PICRUSt2's prediction for the longest 16S locus the MAG keeps
                           whole (``data/mag_arm/panel`` run, one prediction per distinct
                           locus); if none is whole but a >= 250 bp piece survives a
                           fragment boundary, the prediction for that piece's full locus
                           (counted in ``n_partial_proxy``); otherwise abstain;
  * ``picrust2_mag_full``  answers only when a >= 1200 bp locus is kept whole;
  * ``picrust2_ungated``   tier-independent reference: the genome's longest locus, as if
                           the 16S were always recovered;
  * ``prfba_include`` / ``prfba_exclude``  the MAG arm's re-derived selections
                           (``mag_prfba_ec.py``), union of the picks' EC sets;
  * ``base_constant``      calibration floor, as in shared_reference SPEC 8.3: the K most
                           prevalent ECs across the panel, K = median truth size (the
                           scored genome is 1 of ~9,760 contributing to prevalence).
Truth, vocabulary and metric are ``compare.py``'s: the source genome's BV-BRC EC set,
both sides restricted to the 2,132 ECs both systems express, empty truth excluded,
conditional (answered only) and population (abstention = F1 0) means. Strata: whether
the source assembly is in PICRUSt2's reference (its only analogue of self-inclusion; it
cannot self-exclude). Every method is scored on the SAME cells; a cell is dropped for
all methods when its truth or a prFBA pick has no EC projection (counted in ``skipped``).

Writes ``data/mag_score_summary.json`` and ``data/mag_arm/cells.csv.gz``.
"""
import collections
import csv
import gzip
import json
import os
import statistics as st
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mag_common as M  # noqa: E402

SAMPLE = os.path.join(M.DATA, "sample")
PANEL = os.path.join(M.DATA, "panel")
FULL_MIN = 1200


def jacc(a, b):
    return len(a & b) / len(a | b) if (a or b) else 1.0


def score_of(L):
    return max(L.get("bac") or 0.0, L.get("arc") or 0.0)


def longest(cands):
    return max(cands, key=lambda L: (L["len"], score_of(L)))


def load_runs(root):
    """Merge every finished batch of a run_picrust2.sh output tree."""
    pred, meta = {}, {}
    for b in sorted(os.listdir(root)):
        d = os.path.join(root, b)
        if b.startswith("batch_") and os.path.exists(os.path.join(d, "DONE")):
            pred.update(M.p2_ec(d))
            meta.update(M.p2_meta(d))
    return pred, meta


def load_fl16s():
    """{genome_id: (EC set, NSTI)} for the benchmark's FullLength16S amplicons."""
    want = {md5: g for g, md5 in M.fullength_md5().items()}
    pred, meta = load_runs(os.path.join(M.P2C, "data", "full_picrust2"))
    return {want[m]: (ec, meta[m][0]) for m, ec in pred.items() if m in want}


def in_reference():
    """{genome_id: source assembly is one of PICRUSt2's reference genomes} (compare.py)."""
    ref_ec, _ = M.load_reference_ec()
    bv = json.load(open(os.path.join(M.P2C, "data", "bvbrc_source_genomes.json")))
    return {g: (M.acc_core((r or {}).get("assembly_accession")) in ref_ec) for g, r in bv.items()}


def wilcoxon(d):
    try:
        from scipy.stats import wilcoxon as w
    except ImportError:
        return None
    d = [x for x in d if x != 0]
    return float(w(d).pvalue) if len(d) > 10 else None


def mean(vals, n_total=None):
    vals = list(vals)
    if n_total is not None:
        vals += [0.0] * (n_total - len(vals))
    return round(st.mean(vals), 4) if vals else None


# ----------------------------------------------------------------------------- Part A
def part_a(ec, shared, fl16s, panel_md5):
    loci = json.load(open(os.path.join(SAMPLE, "loci.json")))["mags"]
    pred, meta = load_runs(os.path.join(SAMPLE, "picrust2"))
    res = {"rule": "longest locus; ties -> higher HMM bit score"}
    multi = collections.Counter()
    for kind in ("clean", "contaminated"):
        acc = collections.defaultdict(list)
        n = with16 = dropped = 0
        for m, d in loci.items():
            if (d["contamination"] > 0) != (kind == "contaminated"):
                continue
            n += 1
            gid = d["genome_id"]
            truth = ec.get(gid, frozenset()) & shared
            if not d["loci"]:
                continue
            with16 += 1
            L = [x for x in d["loci"] if x["md5"] in pred]
            dropped += len(d["loci"]) - len(L)
            if not L:
                continue
            ch = longest(L)
            p = pred[ch["md5"]] & shared
            if ch["origin"] == "source":
                acc["chosen_is_genome_locus"].append(ch["md5"] in set(panel_md5.get(gid) or ()))
            if kind == "clean":
                src = [x for x in L if x["origin"] == "source"]
                if len(src) > 1:
                    multi["mags_with_multiple_16s"] += 1
                    multi["all_copies_same_prediction"] += (
                        len({pred[x["md5"]] & shared for x in src}) == 1)
                if gid in fl16s:
                    q = fl16s[gid][0] & shared
                    acc["vs_fl16s_" + ("full" if ch["len"] >= FULL_MIN else "partial")].append(
                        (p == q, jacc(p, q), abs(meta[ch["md5"]][0] - fl16s[gid][1])))
                    if truth:
                        acc["f1_pair"].append((M.metrics(truth, p)[3], M.metrics(truth, q)[3]))
            else:
                acc["chosen_" + ch["origin"]].append(1)
                if truth:
                    acc["f1_chosen_" + ch["origin"]].append(M.metrics(truth, p)[3])
                if ch["origin"] == "donor":
                    acc["donor_over_present_source"].append(
                        any(x["origin"] == "source" for x in L))
        out = dict(n_mags=n, with_16s=with16, loci_dropped_by_min_align=dropped)
        if acc["chosen_is_genome_locus"]:
            v = acc["chosen_is_genome_locus"]
            out["chosen_locus_identical_to_a_genome_locus"] = round(sum(v) / len(v), 4)
        for key in ("vs_fl16s_full", "vs_fl16s_partial"):
            v = acc[key]
            if v:
                out[key] = dict(n=len(v), identical_ec_set=round(sum(x[0] for x in v) / len(v), 4),
                                mean_jaccard=round(st.mean(x[1] for x in v), 4),
                                median_abs_nsti_diff=round(st.median(x[2] for x in v), 5))
        if acc["f1_pair"]:
            out["f1_extracted_vs_fl16s"] = dict(extracted=mean(x for x, _ in acc["f1_pair"]),
                                                fl16s=mean(y for _, y in acc["f1_pair"]),
                                                n=len(acc["f1_pair"]))
        for o in ("source", "donor"):
            if acc["chosen_" + o]:
                out["chosen_" + o] = len(acc["chosen_" + o])
            if acc["f1_chosen_" + o]:
                out["f1_when_chosen_" + o] = mean(acc["f1_chosen_" + o])
        if acc["donor_over_present_source"]:
            out["donor_chosen_although_source_16s_present"] = sum(acc["donor_over_present_source"])
        res[kind] = out
    res["multi_copy"] = dict(multi)
    return res


# ----------------------------------------------------------------------------- Part B
def part_b(ec, shared, fl16s, inref, panel_md5):
    loci = {}
    for line in open(os.path.join(M.DATA, "loci.jsonl")):
        r = json.loads(line)
        loci[r["id"]] = r["loci"]
    pred, meta = load_runs(os.path.join(PANEL, "picrust2"))
    retention = {}
    with gzip.open(os.path.join(M.DATA, "retention.csv.gz"), "rt") as fh:
        for r in csv.DictReader(fh):
            if float(r["contamination"]) == 0:
                retention[r["mag"]] = r
    picks = {}
    with gzip.open(os.path.join(M.DATA, "prfba_mag_selections.csv.gz"), "rt") as fh:
        for r in csv.DictReader(fh):
            picks[(r["mag"], r["regime"])] = [g for g in r["picks"].split(";") if g]

    def locus_pred(gid, idx):
        md5 = (panel_md5.get(gid) or [None] * (idx + 1))[idx]
        return pred.get(md5) if md5 else None

    def ungated(gid):
        cands = [(L, i) for i, L in enumerate(loci.get(gid, [])) if L["len"] >= 250]
        if not cands:
            return None
        _, i = max(cands, key=lambda t: (t[0]["len"], score_of(t[0])))
        return locus_pred(gid, i)

    sets = [ec[g] & shared for g in ec if ec[g] & shared]
    prevalence = collections.Counter(e for x in sets for e in x)
    k_const = int(st.median(len(x) for x in sets))
    constant = frozenset(e for e, _ in prevalence.most_common(k_const))

    cw = csv.writer(gzip.open(os.path.join(M.DATA, "cells.csv.gz"), "wt", newline=""))
    cw.writerow(["mag", "genome_id", "completeness", "in_picrust2_ref", "method", "outcome",
                 "n_truth", "n_pred", "recall", "precision", "f1"])
    per = collections.defaultdict(lambda: {"n": 0, "f1": [], "rec": [], "pre": []})
    paired = collections.defaultdict(list)
    skipped, notes = collections.Counter(), collections.Counter()
    for mag, r in sorted(retention.items()):
        gid, tier = r["genome_id"], float(r["completeness"])
        if gid not in ec:
            skipped["truth_unprojected"] += 1
            continue
        truth = ec[gid] & shared
        if not truth:
            skipped["truth_empty"] += 1
            continue
        preds, bad = {}, None
        for regime in ("include", "exclude"):
            pk = picks.get((mag, regime))
            if pk is None:
                bad = "no_prfba_record"
            elif any(g not in ec for g in pk):
                bad = "prfba_pick_unprojected"
            else:
                preds["prfba_" + regime] = (frozenset().union(*(ec[g] for g in pk)) & shared) if pk else None
        if bad:
            skipped[bad] += 1
            continue
        kept = [int(i) for i in r["kept_loci"].split(";") if i != ""]
        p = None
        if kept:
            p = locus_pred(gid, max(kept, key=lambda i: (loci[gid][i]["len"], score_of(loci[gid][i]))))
        elif r["status"] == "partial" and r["piece_locus"] != "":
            p = locus_pred(gid, int(r["piece_locus"]))
            notes["n_partial_proxy"] += 1
        if r["status"] != "none" and p is None:
            notes["kept_locus_without_prediction"] += 1          # dropped by --min_align
        preds["picrust2_mag"] = (p & shared) if p is not None else None
        preds["picrust2_mag_full"] = preds["picrust2_mag"] if r["status"] == "full" else None
        u = ungated(gid)
        preds["picrust2_ungated"] = (u & shared) if u is not None else None
        preds["base_constant"] = constant
        ir = bool(inref.get(gid, False))
        f1s = {}
        for meth, pr in preds.items():
            if pr:
                _, rec, pre, f1, _ = M.metrics(truth, pr)
                f1s[meth] = f1
                cw.writerow([mag, gid, tier, int(ir), meth, "ok", len(truth), len(pr), rec, pre, f1])
            else:
                cw.writerow([mag, gid, tier, int(ir), meth, "abstained", len(truth), 0, 0.0, "", 0.0])
            for stratum in ("all", "in_ref" if ir else "not_in_ref"):
                a = per[(meth, tier, stratum)]
                a["n"] += 1
                if pr:
                    a["f1"].append(f1s[meth])
                    a["rec"].append(M.metrics(truth, pr)[1])
                    a["pre"].append(M.metrics(truth, pr)[2])
        for a_, b_ in (("prfba_exclude", "picrust2_mag"), ("prfba_include", "picrust2_mag")):
            for stratum in ("all", "in_ref" if ir else "not_in_ref"):
                paired[(a_, b_, tier, stratum)].append((f1s.get(a_), f1s.get(b_)))
    table = [dict(method=m, completeness=t, stratum=s, n=a["n"],
                  coverage=round(len(a["f1"]) / a["n"], 4) if a["n"] else None,
                  cond_recall=mean(a["rec"]), cond_precision=mean(a["pre"]),
                  cond_f1=mean(a["f1"]), pop_f1=mean(a["f1"], a["n"]))
             for (m, t, s), a in sorted(per.items(), key=lambda kv: (kv[0][2], kv[0][0], -kv[0][1]))]
    pairs = []
    for (a_, b_, tier, stratum), v in sorted(paired.items(), key=lambda kv: (kv[0][3], kv[0][0], -kv[0][2])):
        both = [(x, y) for x, y in v if x is not None and y is not None]
        d = [x - y for x, y in both]
        pairs.append(dict(a=a_, b=b_, completeness=tier, stratum=stratum, n=len(v),
                          n_both=len(both),
                          a_only=sum(1 for x, y in v if x is not None and y is None),
                          b_only=sum(1 for x, y in v if x is None and y is not None),
                          neither=sum(1 for x, y in v if x is None and y is None),
                          cond_a=mean(x for x, _ in both), cond_b=mean(y for _, y in both),
                          cond_diff=round(st.mean(d), 4) if d else None,
                          a_wins=sum(1 for x in d if x > 0), b_wins=sum(1 for x in d if x < 0),
                          wilcoxon_p=wilcoxon(d),
                          pop_a=mean((x or 0.0) for x, _ in v), pop_b=mean((y or 0.0) for _, y in v)))
    # the benchmark's FullLength16S amplicon, per GENOME, as an external anchor
    fl = [(M.metrics(ec[g] & shared, fl16s[g][0] & shared)[3]) for g in sorted(fl16s)
          if g in ec and (ec[g] & shared) and (fl16s[g][0] & shared)]
    return dict(skipped=dict(skipped), notes=dict(notes), table=table, paired=pairs,
                base_constant_k=k_const,
                benchmark_fl16s_amplicon=dict(n_genomes=len(fl), mean_f1=mean(fl)),
                panel_run_sequences=len(pred))


def main():
    ec, shared = M.load_panel_ec()
    fl16s = load_fl16s()
    inref = in_reference()
    panel_md5 = json.load(open(os.path.join(PANEL, "loci_md5.json")))
    out = {"shared_ec": len(shared), "fl16s_predictions": len(fl16s),
           "panel_in_picrust2_ref": sum(1 for g in M.panel_ids() if inref.get(g))}
    if os.path.exists(os.path.join(SAMPLE, "picrust2", "batch_00", "DONE")):
        out["part_a_sample"] = part_a(ec, shared, fl16s, panel_md5)
    out["part_b_panel"] = part_b(ec, shared, fl16s, inref, panel_md5)
    json.dump(out, open(os.path.join(M.P2C, "data", "mag_score_summary.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "part_b_panel"}, indent=1))
    b = out["part_b_panel"]
    print(json.dumps({k: b[k] for k in ("skipped", "notes", "benchmark_fl16s_amplicon",
                                        "panel_run_sequences")}))
    for r in b["table"]:
        if r["stratum"] == "all" and r["completeness"] in (1.0, 0.5, 0.3, 0.1, 0.01):
            print("%-20s c%-4g n %5d cov %.3f condF1 %s popF1 %s" % (
                r["method"], r["completeness"] * 100, r["n"], r["coverage"], r["cond_f1"], r["pop_f1"]))


if __name__ == "__main__":
    sys.exit(main())
