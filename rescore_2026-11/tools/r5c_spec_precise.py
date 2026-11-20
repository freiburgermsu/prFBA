#!/usr/bin/env python
"""R5C: amend manuscript_numbers/claims_spec.json so number_diff's printed-precision values follow the
manuscript's rounding from the most precise source (idempotent; the pre-amendment spec is
$W/control/r5c/claims_spec.before_r5c.json).

  * `precise` sources for values whose stored artefact is rounded (4 dp) where a finer one exists:
    W6 recomputations (MAG pairings, tab:magval cells, controlled-benchmark bootstrap CIs and AURCs),
    cell means from integer counts (comparison and shared-reference cells), counts behind 4-dp shares,
    and three-panel differences of the unrounded standardized rates;
  * print formats the tex uses (PH-04 p/q at 2 dp, the Spearman p at one decimal of mantissa), the
    PH-02 q range's minimum, the bound rewritten 'within 0.016' -> 'within 0.017';
  * `tex_exception` for values the current tex deliberately prints otherwise;
  * new values printed since R5B with `check_new_print` (PICRUSt2 all-genome phenotype sensitivity,
    standardized surviving-uncovered rates, the unstratified Fisher p).
usage: r5c_spec_precise.py [SPEC]
"""
import json, re, sys

SPEC = sys.argv[1] if len(sys.argv) > 1 else "/home/freiburger/Documents/prFBA/manuscript_numbers/claims_spec.json"
CMP = "region_validation_d100/picrust2_comparison/data/comparison_cells.csv.gz"
CMP_COLS = ["amplicon", "region", "method", "space", "outcome", "n_truth", "n_pred", "tp", "nsti_bin",
            "src_in_picrust2_ref"]
SRG = "region_validation/shared_reference/data/cells.csv.gz"
SRF = "region_validation/shared_reference/data_family/cells.csv.gz"
SR_COLS = ["region", "density", "arm", "space", "outcome", "n_truth", "n_pred", "tp"]
REGIONS = ["FullLength16S", "V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V4-V5_944R", "V7-V9"]
SCOREABLE_SR = {"in": ["scored", "abstained", "pred_empty"]}
PH = "region_validation/shared_reference/data/phenotype_atp_c13/phenotype_summary.json"
EMP = "three_panel/emptied_strata_d100.json"


def cmp_mean(col_derive, flt):
    a = {"usecols": CMP_COLS, "filter": dict({"space": "ec"}, **flt), "column": "x", "fillna": 0.0}
    a["derive"] = {"x": col_derive}
    return {"artifact": CMP, "recipe": "csv_mean", "args": a}


def sr_pop(rel, arm, space="ec", region=None):
    flt = {"density": 1.0, "space": space, "arm": arm, "outcome": SCOREABLE_SR}
    if region:
        flt["region"] = region
    return {"artifact": rel, "recipe": "csv_mean",
            "args": {"usecols": SR_COLS, "derive": {"x": ["f1_counts", "tp", "n_truth", "n_pred"]},
                     "filter": flt, "column": "x", "fillna": 0.0}}


def precise_for(v):
    """A finer source for value v, or None."""
    art, sel = v.get("artifact"), v.get("selector") or []
    # tab:magvs16s columns (paired_comparison_d100.json, 4 dp) -> W6 magvs16s_* (compare_arms pairing)
    if art == "mag_validation/out/paired_comparison_d100.json" and len(sel) == 3:
        m = re.fullmatch(r"exclude\|c1_x0\|([^|]+)\|f1", str(sel[0]))
        if m:
            col = {("conditional", "mag"): "mag", ("conditional", "amp"): "amp",
                   ("conditional", "median_delta"): "median_delta",
                   ("population", "mag"): "pop_mag", ("population", "amp"): "pop_amp"}.get((sel[1], sel[2]))
            if col:
                return {"artifact": "derived:magvs16s_%s_%s" % (m.group(1), col)}
    # tab:magval cells (crossover_d100.json, 4 dp) -> W6 magval_c<comp>_<region>_<col>
    if art == "mag_validation/out/crossover_d100.json" and len(sel) == 2 and isinstance(sel[1], str):
        m = re.fullmatch(r"(FullLength16S|V4)_(mag|amp|delta)", sel[1])
        comp = (1.0 if sel[0] == 0 else 0.01 if sel[0] == -1 else
                sel[0]["where"]["completeness"] if isinstance(sel[0], dict) and "where" in sel[0] else None)
        if m and comp is not None:
            return {"artifact": "derived:magval_c%g_%s_%s" % (comp, m.group(1), m.group(2))}
    # controlled-benchmark EC bootstrap bounds (4 dp) -> W6 sr_ci_* (score_folds' own seeded bootstrap)
    for rank, rel in (("genus", "region_validation/shared_reference/data/summary_genus.json"),
                      ("family", "region_validation/shared_reference/data_family/summary_family.json")):
        if art == rel and len(sel) == 5 and sel[0] == "contrasts" and sel[2] == "ec" and sel[3] == "ci95":
            return {"artifact": "derived:sr_ci_%s_%s_%s" % (rank, sel[1], "lo" if sel[4] == 0 else "hi")}
    if art == "region_validation/shared_reference/data/summary_genus.json":
        if sel == ["arms", "picrust2", "ec", "aurc"]:
            return {"artifact": "derived:sr_aurc_genus_picrust2"}
        if sel == ["arms", "prfba_nestedcv", "ec", "aurc_cross_knob"]:
            return {"artifact": "derived:sr_aurc_genus_prfba_nestedcv_cross_knob"}
    # three-panel differences: est is round(phystd_x - phystd_y, 4); the unrounded rates are stored
    if art in ("three_panel/three_panel_comparison_d100.json", "three_panel/three_panel_mag_panelarm.json") \
            and len(sel) == 4 and sel[-1] == "est":
        m = re.fullmatch(r"delta_(cond_f1_)?(original|kmeans|random)_minus_(original|kmeans|random)", sel[2])
        if m:
            fld = "cond_f1_phystd" if m.group(1) else "phystd"
            return {"expr": "x - y", "vars": {"x": {"artifact": art, "selector": [sel[0], sel[1], m.group(2), fld]},
                                              "y": {"artifact": art, "selector": [sel[0], sel[1], m.group(3), fld]}}}
    # tab:emptied shares (4 dp) -> n / N
    if art == EMP and len(sel) == 5 and sel[0] == "panels" and sel[2] == "strata" and sel[4] == "share":
        return {"expr": "n / N", "vars": {"n": {"artifact": EMP, "selector": ["panels", sel[1], "strata", sel[3], "n"]},
                                         "N": {"artifact": EMP, "selector": ["panels", sel[1], "N"]}}}
    # as-shipped comparison, from the cells' integer counts
    if art == "region_validation_d100/picrust2_comparison/data/comparison_summary.json":
        if sel[:2] == ["by_region", "ec"] and len(sel) == 8 and sel[5:] == ["population", "f1", "value"] and sel[3] == "all":
            flt = {"method": sel[4], "outcome": {"in": ["ok", "abstained"]}}
            if sel[2] != "all":
                flt["region"] = sel[2]
            return cmp_mean(["f1_counts", "tp", "n_truth", "n_pred"], flt)
        if sel[:2] == ["by_region", "ec"] and len(sel) == 8 and sel[5:] == ["conditional", "precision", "value"] and sel[3] == "all":
            flt = {"method": sel[4], "outcome": "ok"}
            if sel[2] != "all":
                flt["region"] = sel[2]
            return dict(cmp_mean(["ratio", "tp", "n_pred"], flt), **{})
        if len(sel) == 8 and sel[:3] == ["by_stratum", "ec", "nsti_bin"] and sel[5:] == ["conditional", "precision", "value"]:
            return cmp_mean(["ratio", "tp", "n_pred"], {"method": sel[4], "outcome": "ok", "nsti_bin": sel[3]})
        if len(sel) == 6 and sel[:3] == ["by_stratum", "ec", "src_in_picrust2_ref"] and sel[-1] == "coverage":
            src = {"False": 0, "True": 1}[sel[3]]
            return {"artifact": CMP, "recipe": "csv_share",
                    "args": {"usecols": CMP_COLS, "filter": {"space": "ec", "method": sel[4], "src_in_picrust2_ref": src,
                                                            "outcome": {"in": ["ok", "abstained"]}},
                             "where": {"outcome": "ok"}}}
    # depth probe: the count behind a 4-dp share (unique for a denominator under 5,000)
    if art == "region_validation/picrust2_comparison/data/depth_probe/depth_probe_results.json" \
            and len(sel) == 3 and sel[0] == "anchor_accuracy_of_recovered" and sel[2] == "correct":
        base = ["anchor_accuracy_of_recovered", sel[1]]
        return {"expr": "round(c * n) / n", "vars": {"c": {"artifact": art, "selector": base + ["correct"]},
                                                     "n": {"artifact": art, "selector": base + ["n_scoreable"]}}}
    # deployment mean selection size (4 dp) -> total genomes / answered (8,597 / 3,755)
    if sel == ["mean_n_genomes_answered"] and str(art).startswith("deployment_d100/stats_d100_"):
        return {"expr": "round(m * a) / a", "vars": {"m": {"artifact": art, "selector": ["mean_n_genomes_answered"]},
                                                     "a": {"artifact": art, "selector": ["answered"]}}}
    return None


def main():
    spec = json.load(open(SPEC))
    items = spec["items"]
    n_prec = 0
    for iid, it in items.items():
        for v in it.get("values", []):
            if "precise" in v:
                continue
            p = precise_for(v)
            if p:
                v["precise"] = p
                n_prec += 1
    # crossover_d100.json stores <region>_delta as the difference of the two 4-dp means (0.814 - 0.7116 =
    # 0.1024 for 0.10245): allow a full unit of the 4th decimal against the unrounded W6 difference
    for iid, it in items.items():
        for v in it.get("values", []):
            if (v.get("artifact") == "mag_validation/out/crossover_d100.json" and len(v.get("selector") or []) == 2
                    and str(v["selector"][1]).endswith("_delta") and "precise" in v):
                v["precise_tol"] = 1e-4
    # old prints that were already wrong, or printed at another precision
    known = {
        ("picrust:SR-08", "  CI lo", "picrust2__vs__base_nn1"):
            ("already a rounding error before the re-score: the stored 4-dp bound 0.0295 was printed 0.030, but the "
             "unrounded bound (score_folds' bootstrap re-run, W6 sr_ci_*) is 0.02949, so 0.029; corrected in the tex "
             "2026-11-17"),
        ("picrust:PH-03", "Spearman p", None):
            "the earlier text printed one significant digit (2x10^-4 for 1.65e-4); the re-scored text prints two (1.5x10^-7)",
        ("panels:disc.refusal_credit", "Spearman p", None):
            "the earlier text printed one significant digit (2x10^-4 for 1.65e-4); the re-scored text prints two (1.5x10^-7)",
    }
    for (iid, label, contrast), why in known.items():
        for v in items[iid]["values"]:
            if v.get("label") == label and (contrast is None or contrast in json.dumps(v.get("selector"))):
                v["known_mismatch"] = why
    # shared-reference cell means (no summary-level finer source): from the cells' integer counts
    def setp(iid, label, p, occurrence=0):
        hits = [v for v in items[iid]["values"] if v.get("label") == label]
        hits[occurrence]["precise"] = p
    setp("picrust:FR-02", "PICRUSt2 pop F1", sr_pop(SRF, "picrust2"))
    setp("picrust:FR-03", "exact-tie vs nn1", {"expr": "a - b", "vars": {"a": sr_pop(SRF, "prfba_exact_tie"),
                                                                       "b": sr_pop(SRF, "base_nn1")}})
    reg_vars = {}
    for i, r in enumerate(REGIONS):
        reg_vars["e%d" % i] = sr_pop(SRF, "prfba_exact_tie", region=r)
        reg_vars["p%d" % i] = sr_pop(SRF, "picrust2", region=r)
    fr03 = [v for v in items["picrust:FR-03"]["values"] if v.get("label") == "exact-tie by region (largest)"][0]
    how = fr03.get("reduce")
    agg = {"min": "min", "max": "max", "maxabs": "min", "minabs": "max"}.get(how, "min")
    setp("picrust:FR-03", "exact-tie by region (largest)",
         {"expr": "%s([%s])" % (agg, ", ".join("e%d - p%d" % (i, i) for i in range(len(REGIONS)))), "vars": reg_vars})
    setp("picrust:SR-07", "exact-tie delta", {"expr": "a - b", "vars": {"a": sr_pop(SRG, "prfba_exact_tie", "rxn"),
                                                                      "b": sr_pop(SRG, "picrust2", "rxn")}})
    # print formats the tex uses
    for iid in ("picrust:PH-03", "panels:disc.refusal_credit"):
        for v in items[iid]["values"]:
            if v.get("label") == "Spearman p":
                v["fmt"] = "e1"
    for iid in ("picrust:PH-04", "panels:disc.atp_vs_nn1"):
        for v in items[iid]["values"]:
            if v.get("label") in ("p", "q"):
                v["fmt"] = "f2"
                v["note"] = "printed at 2 dp (p = 0.41, q = 0.51)"
    for iid in ("picrust:PH-02", "panels:disc.atp_vs_picrust2"):
        vals = items[iid]["values"]
        if not any(v.get("label") == "smallest q of the media with q <= 0.05" for v in vals):
            vals.append({"label": "smallest q of the media with q <= 0.05", "artifact": PH,
                         "selector": ["contrasts", "prfba_exact_tie__vs__picrust2", "controlled", "mcnemar_per_medium",
                                      {"each": "*"}, "q_bh"], "reduce": "min", "skip_null": True,
                         "printed": "0.002", "fmt": "f3",
                         "note": "the low end of the printed q range (0.029--0.046 since the re-score)"})
    b = [v for v in items["panels:panels.condF1_crude_bounds"]["values"]
         if str(v.get("label", "")).startswith("random vs original MAG crude cond F1")][0]
    if b.get("fmt") == "le":
        b["fmt"] = "f3"
        b["note"] = ("was a bound ('within 0.016', fmt le); the re-scored maximum 0.0167 no longer holds it, so the "
                     "sentence now prints the value rounded ('within 0.017')")
    # values the tex deliberately prints otherwise
    exc = {
        ("panels:panels.surviving_uncovered", "FL strata uncovered (low)"):
            "printed in words: 'down to none of the original panel's 274 records in the 50--99 stratum'",
        ("picrust:AS-12", "wins per loss (max)"): "printed in words: 'between three and seven amplicons'",
        ("front:abs-phenotype-controls", "ATP-pass rows (A-C13 configuration, re-run)"):
            ("the abstract keeps the pre-registered pass's 43,140 on purpose (the row-for-row reproduction on the "
             "original sample); the re-run's 43,174 is printed in sec:sharedref"),
    }
    for (iid, label), why in exc.items():
        for v in items[iid]["values"]:
            if v.get("label") == label:
                v["tex_exception"] = why
    # new prints since R5B (no earlier print: printed null, checked against the current tex)
    new_vals = {
        "picrust:PH-01": [
            {"label": "all-genome hamming delta (n 1,000)", "artifact": PH,
             "selector": ["contrasts", "prfba_exact_tie__vs__picrust2", "controlled", "delta_atp_hamming_all_genomes"],
             "printed": None, "fmt": "f2", "check_new_print": True},
            {"label": "all-genome Wilcoxon p", "artifact": PH,
             "selector": ["contrasts", "prfba_exact_tie__vs__picrust2", "controlled", "wilcoxon_p_all_genomes"],
             "printed": None, "fmt": "f3", "check_new_print": True}],
        "panels:panels.surviving_uncovered": [
            {"label": "stratum-standardized uncovered rate %s" % p, "artifact": EMP,
             "selector": ["tests", "standardised_non_ok_rate", p], "printed": None, "fmt": "pct1",
             "check_new_print": True} for p in ("original", "kmeans", "random")] + [
            {"label": "unstratified Fisher p original vs random", "artifact": EMP,
             "selector": ["tests", "survivors_pairwise", "original_vs_random", "fisher_p"], "printed": None,
             "fmt": "f3", "check_new_print": True}],
    }
    # the legacy reducer's as-shipped population figures, printed since R5B (RUN5B Flag 1)
    S100 = "region_validation_d100/picrust2_comparison/data/comparison_summary.json"
    S20 = "region_validation/picrust2_comparison/data/comparison_summary.json"
    C20 = "region_validation/picrust2_comparison/data/comparison_cells.csv.gz"

    def pop(summary, method):
        return {"artifact": summary, "selector": ["by_region", "ec", "all", "all", method, "population", "f1", "value"]}

    def pop_cells(cells, method):
        return {"artifact": cells, "recipe": "csv_mean",
                "args": {"usecols": CMP_COLS, "derive": {"x": ["f1_counts", "tp", "n_truth", "n_pred"]},
                         "filter": {"space": "ec", "method": method, "outcome": {"in": ["ok", "abstained"]}},
                         "column": "x", "fillna": 0.0}}
    gap = {"label": "legacy reducer minus PICRUSt2, as-shipped population EC F1 (depth 100)",
           "expr": "a - b", "vars": {"a": pop(S100, "prfba_exclude_legacy"), "b": pop(S100, "picrust2")},
           "precise": {"expr": "a - b", "vars": {"a": pop_cells(CMP, "prfba_exclude_legacy"),
                                                 "b": pop_cells(CMP, "picrust2")}},
           "precise_tol": 1e-4,          # the stored side is a difference of two 4-dp values
           "printed": None, "fmt": "+f3", "check_new_print": True}
    new_vals["picrust:AS-09"] = [
        dict(pop(S100, "prfba_exclude_legacy"), label="legacy reducer population EC F1, as shipped (depth 100)",
             precise=pop_cells(CMP, "prfba_exclude_legacy"), printed=None, fmt="f3", check_new_print=True),
        dict(pop(S20, "prfba_exclude_legacy"), label="legacy reducer population EC F1, as shipped (depth 20)",
             precise=pop_cells(C20, "prfba_exclude_legacy"), printed=None, fmt="f3", check_new_print=True),
        dict(gap)]
    new_vals["panels:disc.margin_reversal"] = [dict(gap)]
    for iid, vals in new_vals.items():
        cur = items[iid]["values"]
        for nv in vals:
            for i, v in enumerate(cur):
                if v.get("label") == nv["label"]:
                    cur[i] = nv
                    break
            else:
                cur.append(nv)
    spec.setdefault("_meta", {})["r5c"] = ("2026-11-17: precise sources, print formats, tex exceptions and new "
                                          "prints ($W/tools/r5c_spec_precise.py)")
    with open(SPEC + ".tmp", "w") as fh:
        json.dump(spec, fh, indent=1)
    import os
    os.replace(SPEC + ".tmp", SPEC)
    print("precise sources added:", n_prec, "(+4 shared-reference cell means)")


if __name__ == "__main__":
    main()
