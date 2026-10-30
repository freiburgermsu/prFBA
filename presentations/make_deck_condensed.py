#!/usr/bin/env python
"""Build the condensed prFBA deck: the method in four figures, then the validation results.

A short companion to the full update deck (make_deck.py). It shares that deck's loader,
palette, helpers and six of its result slides, so the two cannot disagree; the methods
slides and two result slides are new figures drawn from the pipeline's own files by
presentations/condensed/make_condensed_figures.py (run that first).

    python presentations/condensed/make_condensed_figures.py
    python presentations/make_deck_condensed.py      # -> presentations/prFBA_condensed_2026-10.pptx
    python presentations/make_deck_condensed.py --brief   # -> presentations/prFBA_brief_2026-10.pptx

--brief halves it again (9 slides): the workflow, the validation design and PICRUSt2's runs,
then where the calls land, why they fail, the controlled comparison, the phenotypes and the
summary; no section dividers.

Depth-100 records and the exact-tie union throughout, as in the full deck.
"""
import argparse
import json
import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt

import make_deck as MD
from make_deck import (M, W, H, BODY_TOP, DEEP, INK, TEAL, MINT, AMBER, AMBER_TXT, WHITE,
                       MUTED, HEAD, blank, background, title, footer, notes, para, rect,
                       textbox, stat, rd, sg)

FIG = os.path.join(MD.REPO, "presentations", "condensed")


def figure(s, name, *, y=BODY_TOP + 0.05, w=W - 2 * M):
    path = os.path.join(FIG, name)
    pic = s.shapes.add_picture(path, Inches(M), Inches(y), width=Inches(w))
    pic.name = name.rsplit(".", 1)[0]
    return pic


def section(prs, head, sub):
    s = blank(prs)
    background(s, DEEP)
    tf = textbox(s, M, 2.6, W - 2 * M, 2.0)
    para(tf, head, size=44, bold=True, color=WHITE, font=HEAD, first=True, space_after=4)
    para(tf, sub, size=20, color=MINT, font=HEAD, space_after=0)
    return s


# ------------------------------------------------------------------------ slides
def c_title(prs, D, F):
    s = blank(prs)
    background(s, DEEP)
    tf = textbox(s, M, 1.35, W - 2 * M, 2.2)
    para(tf, "prFBA", size=60, bold=True, color=WHITE, font=HEAD, first=True, space_after=2)
    para(tf, "Probabilistic flux balance analysis for amplicon-resolved metabolic modeling",
         size=22, color=MINT, font=HEAD, space_after=0)
    chips = ("One amplicon, end to end", "Five validation arms", "Like for like, PICRUSt2 leads")
    gap = 0.3
    cw = (W - 2 * M - gap * (len(chips) - 1)) / len(chips)
    for i, t in enumerate(chips):
        c = rect(s, M + i * (cw + gap), 4.0, cw, 0.7, RGBColor(0x0A, 0x4A, 0x55),
                 line=MD.SEAFOAM)
        ctf = c.text_frame
        ctf.margin_top = ctf.margin_bottom = 0
        ctf.vertical_anchor = MSO_ANCHOR.MIDDLE
        para(ctf, t, size=16, bold=True, color=WHITE, align=PP_ALIGN.CENTER, first=True,
             space_after=0)
    tf2 = textbox(s, M, 5.9, W - 2 * M, 0.4)
    para(tf2, "Andrew Freiburger  ·  Tyo & Henry labs  ·  Condensed update, October 2026",
         size=13, color=MINT, first=True, space_after=0)
    notes(s, (
        "A condensed version of the October update deck: the method, each slide built "
        "on one figure, then the validation results. Every benchmark number is at the "
        f"standard record depth of {D['aln_ms']['params']['k2']} hits under the shipped "
        "exact-tie union selection rule, and every value is read from the pipeline's result "
        "files by the same loader the full deck uses. The full deck carries the alignment "
        "engineering, the C12 diagnosis, the depth and rule decisions, the MAG build details "
        "and the audit trail."))
    return s


def c_workflow(prs, D, F):
    w = F["workflow"]
    ex = json.load(open(os.path.join(FIG, "worked_example.json")))
    s = blank(prs)
    title(s, "One amplicon, from 16S read to weighted genome",
          sub="A real V4 benchmark amplicon, source removed, through every stage")
    figure(s, "F_workflow.png")
    footer(s, f"{ex['amplicon']} · {MD.DEPTH}-hit record · exact-tie union · self-exclusion "
              "at τ = 0.987 · ECs from BV-BRC PGFams")
    n = len(w["reps"])
    pp = w["present_by_p"]
    notes(s, (
        "The workflow on one real benchmark amplicon, chosen by a fixed rule rather than by "
        "hand (extract_worked_example.py: the first V4 amplicon whose source is its own single "
        f"selection, and which falls to a {n}-genome genus-tier tie once the source is "
        f"removed). The amplicon is a V4 product of a {ex['truth_lineage']['Genus']} genome "
        f"({ex['truth_lineage']['Family']}, {ex['truth_lineage']['Phylum']}). Step 1: an edlib "
        f"prefilter keeps the {ex['n_edlib_candidates']} closest of 459,301 unique BV-BRC 16S "
        f"sequences, the GPU Smith-Waterman kernel ranks them, and Biopython aligns the top "
        f"{ex["n_hits"]} to give identity and coordinates "
        "(panel A). Step 2: self-exclusion removes the source genome and every hit at identity "
        f"0.987 or above ({w['n_removed']} hits here), then the exact-tie rule keeps every "
        f"reference tied at the best surviving identity, {w['best_identity']:.4f}, one genome "
        f"per species: {n} {ex['truth_lineage']['Genus']} genomes. A Bacillus 16S sits in the "
        "record just below the tie at rank 8, an off-phylum sequence the identity tie leaves "
        "out. Step 3: the selected genomes' functions are merged "
        f"into one synthetic genome of {w['n_union']} ECs (panel B), each weighted by p = k/n, "
        f"the share of the {n} genomes carrying it: "
        + ", ".join(f"{c} at p = {k}" for k, c in w["n_by_p"].items()) + ". "
        "Step 4: scored against the source genome's own ECs (panel C), the union reaches EC F1 "
        f"{w['f1']:.2f} (precision {w['precision']:.2f}, recall {w['recall']:.2f}, against "
        f"{w['n_truth']} true ECs). Within this one amplicon the share of ECs truly present "
        "falls with p: " + ", ".join(f"{v:.2f} at p = {k}" for k, v in pp.items()) + ". "
        "That is the ranking the calibration slide finds across the benchmark. In prFBA the "
        "weight then enters the flux objective, where a low-p reaction costs more flux; the "
        "weight never gates or removes a reaction. This example is typical in kind, not in "
        "outcome: most answered benchmark amplicons select one genome, so p = 1 throughout."))
    return s


def c_arms(prs, D, F):
    s = blank(prs)
    title(s, "Five validation arms, one unaltered pipeline",
          sub="Each changes what is searched, what is held out, and what counts as truth")
    figure(s, "F_arms.png")
    footer(s, "Every arm runs the shipped code; only the inputs and the scoring differ")
    g = D["kg"]
    notes(s, (
        "The validations share one pipeline and differ in four things: the query, the "
        "reference it is searched against, what is withheld, and the truth it is scored "
        "against. In-silico amplicons: "
        f"{D['tax']['_meta']['n_total']:,} amplicons cut by in-silico PCR from 10,000 "
        "taxonomy-stratified BV-BRC genomes in 8 regions (chart: amplicons per region), "
        "searched against the full reference; each is scored with its source present and with "
        "the source and every hit at 0.987 identity or above removed, against the source's "
        "NCBI lineage and gene families. As shipped against PICRUSt2: the same amplicons, each "
        f"system with its own reference ({D['cmp']['_meta']['picrust2']['n_reference_genomes']:,} "
        "genomes for PICRUSt2 against about 315,000 for prFBA), scored in BV-BRC's EC "
        "vocabulary, which is prFBA's own; PICRUSt2 has no way to exclude the source. "
        f"Controlled: {g['meta']['n_eligible']:,} dereplicated panel genomes split into "
        f"{g['meta']['n_folds']} folds with whole genera held out (chart: held-out genomes per "
        "fold, each stripe a genus group; one connected group fills fold 1), both methods "
        f"confined to one shared reference ({g['headline']['reference_density_genomes']:,} "
        "genomes, median per fold) and one answer key, plus ATP and gapfill phenotypes. "
        f"Synthetic MAGs: {F['arms']['n_drafts']:,} drafts made by fragmenting each panel "
        "genome and dropping whole pieces down to 1% completeness (schematic), matched by "
        "skani whole-genome ANI. Three panels: the same pipeline on three 10,000-organism "
        "panels sampled on different principles (chart: share of each panel's picks per "
        "density decile of the 16S embedding), to ask whether coverage is a property of the "
        "method or of the panel."))
    return s


def c_sankey(prs, D, F):
    k = F["sankey"]
    ex, inc, n, sc = k["exclude"], k["include"], k["n"], k["scoreable"]
    s = blank(prs)
    title(s, "Where every amplicon's genus call lands",
          sub=f"{n:,} amplicons by region · source present (left) vs removed (right)")
    figure(s, "F_sankey.png")
    footer(s, "Anchor predictor · exact-tie union · 100-hit records · shares of the "
              f"{sc:,} genus-scoreable amplicons")
    pc = lambda v: f"{100 * v / sc:.0f}%"  # noqa: E731
    notes(s, (
        f"Every one of the {n:,} benchmark amplicons flows from its region (middle) to its "
        "genus outcome, with the source organism present (left) and removed (right, self-"
        "exclusion at τ = 0.987). "
        f"{ex['unscoreable']:,} amplicons ({100 * ex['unscoreable'] / n:.0f}% of all) come "
        "from a source whose own truth has no named genus; they are unscoreable in both "
        f"modes and leave every denominator. Of the remaining {sc:,}: with the source present "
        f"{inc['correct']:,} ({pc(inc['correct'])}) get the right genus, {inc['wrong']:,} "
        f"({pc(inc['wrong'])}) a wrong one and {inc['unnamed']:,} ({pc(inc['unnamed'])}) a "
        "reference with no named genus, and none abstains. With the source removed, "
        f"{ex['correct']:,} ({pc(ex['correct'])}) are still right, {ex['wrong']:,} "
        f"({pc(ex['wrong'])}) wrong, {ex['unnamed']:,} ({pc(ex['unnamed'])}) land on an "
        f"unnamed genus and {ex['abstained']:,} ({pc(ex['abstained'])}) abstain: for "
        f"{k['abstain_emptied']:,} of them self-exclusion removed all {MD.DEPTH} hits of the "
        "record, and the rest had nothing at the family floor or past the reliability gate. "
        "Ribbon colour is the outcome, so the change in "
        "colour from left to right is what novelty costs, region by region. These counts "
        "reproduce the full deck's denominator slide (asv_region_concordance.csv and "
        "benchmark_truth.json, classified identically)."))
    return s


def c_picrust2(prs, D, F):
    m, p2 = D["cmp"]["_meta"], D["cmp"]["_meta"]["picrust2"]
    g = D["kg"]
    fz = g["meta"]["frozen"]["picrust2"]
    s = blank(prs)
    title(s, "PICRUSt2 run twice: as shipped, then like for like",
          sub="The comparison removes each built-in advantage in turn")
    figure(s, "F_picrust2.png")
    footer(s, f"PICRUSt2 {p2['version']} · installed without conda (uv venv, castor, EPA-ng, "
              "gappa, HMMER from source)")
    notes(s, (
        f"A, as shipped: PICRUSt2 {p2['version']} aligns each amplicon to its reference with "
        "HMMER profiles, places it in a fixed reference phylogeny with EPA-ng and gappa, and "
        f"infers gene content by {p2['hsp'].upper()} (maximum-parsimony) ancestral-state "
        f"reconstruction with castor over {p2['n_reference_genomes']:,} reference genomes; it "
        "always answers. prFBA aligns the same amplicon to 459,301 BV-BRC 16S sequences, keeps "
        "the exact ties and unions their gene families. prFBA predicts PGFam families and "
        "PICRUSt2 EC numbers, so the truth and both predictions pass through one map, family "
        "to BV-BRC product string to EC, and are scored on the "
        f"{m['ec_universe_shared']:,} ECs both can express ({m['rxn_universe_shared']:,} "
        "ModelSEED reactions). Two asymmetries are built in: the references differ about "
        "tenfold, and the truth is BV-BRC annotation, prFBA's own vocabulary. "
        f"B, controlled: {g['meta']['n_eligible']:,} panel genomes in "
        f"{g['meta']['n_folds']} genus-blocked folds; for each fold, the training genera "
        "become one species-dereplicated reference "
        f"({g['headline']['reference_density_genomes']:,} genomes, median), and PICRUSt2's "
        "alignment, tree and trait table are rebuilt from it with its own tooling (FastTree "
        f"{' '.join(g['meta']['frozen']['fasttree_args'])}, HSP {fz['hsp_method']}, edge "
        f"exponent {fz['edge_exponent']}, minimum alignment {fz['min_align']}), while prFBA "
        "searches exactly the same genomes, with its selection rule chosen by nested cross-"
        "validation inside each training fold. Both are scored on the held-out genera against "
        "one answer key, and the gapfilled models are compared on ATP production and growth. "
        f"The effect size that counts was fixed in advance at {g['meta']['frozen']['sesoi_delta_f1']} "
        "in population EC F1."))
    return s


def c_heatmap(prs, D, F):
    s = blank(prs)
    title(s, "Without the source, family transfers best",
          sub="Every region × measure · source present vs removed")
    figure(s, "F_heatmap.png", y=BODY_TOP + 0.1)
    footer(s, "Anchor predictor · exact-tie union · 100-hit records · ÷ answered = among "
              "amplicons that answer; ÷ all = every non-answer a miss")
    tax = D["tax"]["regions"]
    regs = list(tax)
    pr = lambda r, mo, rank, q: tax[r]["anchor"][mo]["all"]["all"]["per_rank"][rank][q]  # noqa: E731
    rng = lambda v: f"{rd(min(v), 2)}–{rd(max(v), 2)}"  # noqa: E731
    gc = [pr(r, "exclude_0.987", "Genus", "correct") for r in regs]
    fc = [pr(r, "exclude_0.987", "Family", "correct") for r in regs]
    gp = [pr(r, "exclude_0.987", "Genus", "correct_population") for r in regs]
    cap = D["cap"]["_meta"]
    notes(s, (
        "Each row is a 16S region and each cell a measure, with the source organism present, "
        "removed and scored among amplicons that still answer, and removed and scored over "
        "every scoreable amplicon (a non-answer counts as a miss; bold). With the source "
        f"present, almost everything is recovered. Removed, conditional genus accuracy falls "
        f"to {rng(gc)}, family accuracy holds at {rng(fc)}, and gene capture to "
        + rng([v / 100 for v in cap['headline_macro_recall_pct']['exclude_0.987'].values()])
        + "; over the whole population genus accuracy is " + rng(gp) + ". Full-length 16S "
        "and the two V1-anchored regions lose least; the regions starting at V3 or later lose "
        "a few points more genus accuracy, while gene capture barely depends on the region. "
        "Source: taxacc_region_summary.json and gene_capture_by_region.json (manuscript "
        "tab:main)."))
    return s


def c_space(prs, D, F):
    sp = F["space"]
    s = blank(prs)
    title(s, "Failures map to where the 16S space is crowded",
          sub="Original panel's sources in the BV-BRC 16S embedding · source removed")
    figure(s, "F_space.png", y=BODY_TOP + 0.1)
    footer(s, "UMAP of 1024-d 16S embeddings; density = 21st-neighbour distance in the full "
              "space · exact-tie union · 100-hit records")
    a_fl, a_v4 = sp["FullLength16S_abst_by_decile"], sp["V4_abst_by_decile"]
    c_fl, c_v4 = sp["FullLength16S_acc_by_decile"], sp["V4_acc_by_decile"]
    n_fl = sp["FullLength16S_n_by_decile"]
    notes(s, (
        f"Each dot is one of the {sp['n_located']:,} original-panel source genomes that sit on "
        "the embedding frame, placed at its representative 16S in the UMAP projection of all "
        "clusterable BV-BRC 16S sequences, and coloured by its amplicon's genus call once the "
        "source is removed (A full-length 16S, B V4); grey shading is local density measured in "
        "the full 1024-dimensional embedding. Panel C bins the sources by that density. In the "
        f"sparsest decile full-length 16S still calls the right genus {c_fl[0]:.0%} of the "
        f"time and abstains {a_fl[0]:.0%}; in the densest it abstains {a_fl[-1]:.0%} "
        f"(V4: {a_v4[-1]:.0%}) and is right only {c_fl[-1]:.0%}. That is not a failure to "
        "find relatives: in densely sequenced species every one of the record's 100 hits is "
        "the source's own taxon or within 0.987 identity of it, so self-exclusion empties the "
        f"record and nothing is left to select ({F['sankey']['abstain_emptied']:,} of the "
        f"{F['sankey']['exclude']['abstained']:,} abstentions over all regions; the full "
        "deck's emptied-records slide). The densest deciles hold few "
        f"panel sources ({n_fl[-2]} and {n_fl[-1]} for full-length 16S), because the "
        "taxonomy-stratified panel under-samples crowded space. V4 is right less often than "
        "full-length 16S at every density. Coverage figures under self-exclusion are therefore "
        "lower bounds set by the record depth, not by the organisms."))
    return s


def c_close(prs, D, F):
    s = blank(prs)
    background(s, DEEP)
    title(s, "What the validations say", dark=True)
    tax = D["tax"]["regions"]
    regs = list(tax)
    fc = [tax[r]["anchor"]["exclude_0.987"]["all"]["all"]["per_rank"]["Family"]["correct"]
          for r in regs]
    cg = D["kg"]["contrasts"]["prfba_nestedcv__vs__picrust2"]["ec"]
    k = F["sankey"]
    c12 = D["ctl"]["C12"]["status"]
    rows = [
        (f"{rd(min(fc), 2)}–{rd(max(fc), 2)}", "family accuracy for novel organisms, among answered"),
        (f"{100 * k['exclude']['abstained'] / k['scoreable']:.0f}%",
         "abstain once the source is removed: no false metabolism"),
        (sg(cg["delta_f1_population"]),
         f"prFBA − PICRUSt2, like for like · CI {sg(cg['ci95'][0])} to {sg(cg['ci95'][1])}"),
        ("PASS" if c12 == "pass" else "FAIL", "gapfill validity (C12), after the template fixes"),
    ]
    card = RGBColor(0x0A, 0x4A, 0x55)
    for i, (v, lab) in enumerate(rows):
        y = BODY_TOP + 0.1 + i * 1.22
        b = rect(s, M, y, W - 2 * M, 1.0, card)
        tf = textbox(s, M + 0.3, y + 0.08, 3.3, 0.85, anchor=MSO_ANCHOR.MIDDLE)
        para(tf, v, size=32, bold=True, color=AMBER if i == 2 else MINT, font=HEAD,
             first=True, space_after=0)
        tf = textbox(s, M + 3.8, y + 0.08, W - 2 * M - 4.1, 0.85, anchor=MSO_ANCHOR.MIDDLE)
        para(tf, lab, size=19, color=WHITE, first=True, space_after=0)
    footer(s, "github.com/freiburgermsu/prFBA  ·  full update deck: "
              "presentations/prFBA_updates_2026-10.pptx", dark=True)
    notes(s, (
        "Four things to take away. First, the selection transfers at family level even for "
        f"organisms absent from the database ({rd(min(fc), 2)}–{rd(max(fc), 2)} family "
        "accuracy across regions, among amplicons that answer), and genus is where short "
        "amplicons pay. Second, it abstains rather than guess: "
        f"{100 * k['exclude']['abstained'] / k['scoreable']:.0f}% of scoreable amplicons get "
        f"no call once the source is removed, {k['abstain_emptied']:,} of "
        f"{k['exclude']['abstained']:,} because every hit in the record was the source's own "
        "close kin, removed with it. Third, as shipped prFBA is more precise when it answers, but on the "
        "controlled benchmark, with one shared reference and a neutral answer key, PICRUSt2 "
        f"predicts EC content better: population EC F1 differs by {sg(cg['delta_f1_population'])} "
        f"(95% CI {sg(cg['ci95'][0])} to {sg(cg['ci95'][1])}), beyond the pre-registered "
        "0.02. Fourth, the phenotype arm's gapfill validity gate "
        f"{'passes' if c12 == 'pass' else 'fails'} after the disclosed template amendments, "
        "so the leave-genus-out run is valid; no EC number depended on it."))
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=None)
    ap.add_argument("--brief", action="store_true", help="the 9-slide version")
    args = ap.parse_args()
    if args.out is None:
        args.out = os.path.join(MD.REPO, "presentations", "prFBA_%s_2026-10.pptx"
                                % ("brief" if args.brief else "condensed"))
    MD.configure_depth(100, "union", None)
    D = MD.load()
    cm = D["ctl"]["C12"].get("measured") or {}
    D["c12_rescored"] = bool(D["c12fix"] and ("phenotype_gf" + D["amend_tag"] + "/") in
                             (cm.get("binding") or {}).get("phenotype_summary", ""))
    F = json.load(open(os.path.join(FIG, "figure_stats.json")))
    MD.NEW_SLIDES = set()                  # "new since the last deck" means nothing here
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    plan = [c_title,
            lambda p, D, F: section(p, "Methods", "Align · select · merge · validate"),
            c_workflow, c_arms, c_sankey, c_picrust2,
            lambda p, D, F: section(p, "Results", "What survives when the source is removed"),
            c_heatmap, c_space, MD.s11_results, MD.s16_kfold_result, MD.s17_phenotype,
            MD.s06c_calibration, MD.s06d_wobj, MD.s17_crossover, MD.s22_three_panel, c_close]
    if args.brief:
        plan = [c_title, c_workflow, c_arms, c_picrust2, c_sankey, c_space,
                MD.s16_kfold_result, MD.s17_phenotype, c_close]
    for fn in plan:
        MD.CURRENT_SLIDE = getattr(fn, "__name__", "section")
        MD._STATS.clear()
        if fn.__module__ == "make_deck":
            fn(prs, D)
        else:
            fn(prs, D, F)
        if MD._STATS:
            size = min(v for _, v in MD._STATS)
            for run, _ in MD._STATS:
                run.font.size = Pt(size)
    prs.core_properties.title = "prFBA: method and validation (%s)" % (
        "brief" if args.brief else "condensed")
    prs.core_properties.author = "Andrew Freiburger"
    prs.save(args.out)
    print(f"[condensed] {len(prs.slides._sldIdLst)} slides -> {args.out}")


if __name__ == "__main__":
    main()
