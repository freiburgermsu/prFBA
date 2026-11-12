#!/usr/bin/env python
"""Build the condensed prFBA deck: the method in four figures, then the validation results.

A short companion to the full update deck (make_deck.py). It shares that deck's loader,
palette, helpers and nine of its result slides, so the two cannot disagree; the methods
slides and two result slides are new figures drawn from the pipeline's own files by
presentations/condensed/make_condensed_figures.py (run that first).

    python presentations/condensed/make_condensed_figures.py
    python presentations/make_deck_condensed.py      # -> presentations/prFBA_condensed_2026-10.pptx
    python presentations/make_deck_condensed.py --brief   # -> presentations/prFBA_brief_2026-10.pptx
    python presentations/make_deck_condensed.py --panels  # -> presentations/prFBA_brief_panels_2026-10.pptx

--brief halves it again (10 slides): the workflow, the validation design and PICRUSt2's runs,
then where the calls land, why they fail, the controlled comparison, the phenotypes, a table
of what is new since the deck of 2026-11-01 and the summary; no section dividers.

--panels is the brief deck with one PICRUSt2 slide (the controlled result) instead of three;
the freed slides detail the studies on the 10,000-organism panels: the 16S and MAG arms'
design, the three-panel comparison and the paired MAG-versus-16S result.

Since 2026-11-11 the condensed deck also carries the full deck's 0.995 self-exclusion and
PICRUSt2-on-MAGs slides and one slide joining Stage E (a better p) and Stage O (the tuned
objective), each with the NEW chip; every variant's closing slide lists what is open, with
the integration item and the open items drawn from make_deck (integration_status(),
integration_sentence(), OPEN_ITEMS) so one edit updates all decks.

Depth-100 records and the exact-tie union throughout, as in the full deck.
"""
import argparse
import json
import os

from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.dml import MSO_LINE_DASH_STYLE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt

import make_deck as MD
from make_deck import (M, W, H, BODY_TOP, DEEP, INK, TEAL, SEAFOAM, MINT, AMBER, AMBER_TXT,
                       WHITE, MUTED, HEAD, blank, background, title, footer, notes, para, rect,
                       textbox, stat, rd, sg, takeaway, table, arrow_value)

FIG = os.path.join(MD.REPO, "presentations", "condensed")
VARIANT = "condensed"                       # condensed | brief | panels; set by main()

# Slides new since the condensed decks of 2026-11-01 carry the full deck's amber NEW chip
# (make_deck.title() reads MD.NEW_SLIDES, and enforces its 55-character title rule on them).
# The brief variants' one new slide is titled "New since Nov 1", so it needs no chip.
NEW_BY_VARIANT = {"condensed": {"s19b_picrust2_mags", "s09b_tau995", "c_estobj"},
                  "brief": set(), "panels": set()}

# What is open, as on the full deck's s15_next: the integration item is
# make_deck.integration_status(D) (the placeholder the integration pass replaces there, drawn
# dashed amber and italic while it is one), then the short forms of make_deck.OPEN_ITEMS. The
# integration's notes sentence is make_deck.integration_sentence(D).


def retarget_notes(s, old, new):
    """Rewrite a cross-reference in a reused full-deck slide's notes whose neighbour this
    deck does not carry; stop if the wording has changed, so a stale pointer cannot ship."""
    tf = s.notes_slide.notes_text_frame
    if tf.text.count(old) != 1:
        raise SystemExit(f"{MD.CURRENT_SLIDE}: notes no longer contain {old!r} exactly once")
    tf.text = tf.text.replace(old, new)


def bracket_sentence(D):
    """The real-ASV bracket (the full deck's s09c_bracket), as one parenthesis for s09b's
    notes: shares of each deployment's ASVs per best-hit identity band."""
    br = D["bracket"]
    pc = {k: [rd(100 * v / br[k][0], 0) for v in br[k][1:4]] for k in ("emilykin", "codif")}
    e, c = pc["emilykin"], pc["codif"]
    ah = {k: "/".join(rd(100 * v / MD.BRACKET_ANYHIT[k][0], 0) for v in MD.BRACKET_ANYHIT[k][1:])
          for k in ("emilykin", "codif")}
    return (f"(the full deck's bracket slide: of the real ASVs of the two depth-100 deployments, "
            f"EmilyKin's {br['emilykin'][0]:,} and the codiffusion bioreactor's "
            f"{br['codif'][0]:,}, {e[0]}% and {c[0]}% have their top hit, the highest-scoring of "
            f"the record's 100, below 0.987, {e[1]}% and {c[1]}% in the 0.987–0.995 band, and "
            f"{e[2]}% and {c[2]}% at 0.995 or above; column best_identity of each deployment's "
            "bvbrc_alignment_hits_d100/asv_alignment_summary.csv. Binned by the highest identity "
            f"of any hit in the record, which is what the exclusion removes, EmilyKin gives "
            f"{ah['emilykin']}% and codiffusion {ah['codif']}%)")


# cross-references in reused full-deck notes that point at slides this deck does not carry
NOTE_FIXES = {
    "s09b_tau995": lambda D: [("(next slide)", bracket_sentence(D))],
    "s19b_picrust2_mags": lambda D: [("(previous slide)",
                                      "(the full update deck's MAG-input slide)")],
    # the C12 diagnosis and fix slides are only in the full deck
    "s16_kfold_result": lambda D: [("the next three slides show why",
                                    "the full update deck shows why")],
    "s17_phenotype": lambda D: [("the next two slides show the",
                                 "the full update deck shows the")],
    # the MAG arm's abstention slide is only in the full deck
    "s17_crossover": lambda D: [("(see the abstention slide)",
                                 "(see the full update deck's MAG abstention slide)")],
}


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
    chips = (("One amplicon, end to end", "16S and MAGs on 3 × 10,000 organisms",
              "Like for like, PICRUSt2 leads") if VARIANT == "panels" else
             ("One amplicon, end to end", "Five validation arms", "Like for like, PICRUSt2 leads"))
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
    para(tf2, "Andrew Freiburger  ·  Tyo & Henry labs  ·  Condensed update, October–November "
         "2026", size=13, color=MINT, first=True, space_after=0)
    notes(s, (
        "A condensed version of the October–November update deck: the method, each slide built "
        "on one figure, then the validation results. Every benchmark number is at the "
        f"standard record depth of {D['aln_ms']['params']['k2']} hits under the shipped "
        "exact-tie union selection rule, and every value is read from the pipeline's result "
        "files by the same loader the full deck uses. The full deck carries the alignment "
        "engineering, the C12 diagnosis, the depth and rule decisions, the MAG build details "
        "and the audit trail. "
        + ("The slides marked NEW were added on 2026-11-11 and 2026-11-12 and carry the "
           "results since the deck of 2026-11-01: self-exclusion at 0.995, the better "
           "p(function) and tuned objective, and PICRUSt2 on MAGs."
           if VARIANT == "condensed" else
           "The 'New since Nov 1' slide, added on 2026-11-12, lists the results since the "
           "deck of 2026-11-01.")))
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
        + ("That is the ranking the calibration slide finds across the benchmark. "
           if VARIANT == "condensed" else
           "That is the ranking the full deck's calibration study finds across the benchmark. ")
        + "In prFBA the "
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


def c_panels(prs, D, F):
    f = F["panels"]
    s = blank(prs)
    title(s, "16S and MAGs on the same 10,000-organism panels",
          sub="Three panels, two ways to query each organism, one truth: its own genome")
    figure(s, "F_panels.png")
    footer(s, "Panels in the BV-BRC 16S embedding (UMAP) · 16S at record depth "
              f"{MD.DEPTH}, exact-tie union · MAGs by skani ANI")
    n = f["n_amp"]
    lev = f["mag_panel_levels"]
    g = f["n_gtdb"]
    seed_og = json.load(open(os.path.join(MD.AMP_RV, "data", "selection_10k.json")))["meta"]["seed"]
    notes(s, (
        "Every validation in this part starts from a panel of 10,000 BV-BRC organisms whose "
        "genomes are known, so each prediction can be scored against the organism's own "
        "lineage and gene content. Three panels are sampled on different principles (left, "
        "each panel's picks in the UMAP projection of the 16S embedding): the original panel "
        f"is taxonomy-stratified (seed {seed_og}); the k-means panel takes "
        "one organism per cluster of the 16S embedding; and the random panel draws "
        f"uniformly (seed {D['random_meta']['seed']}) from the "
        f"{D['random_meta']['n_admissible_genomes']:,} admissible organisms, "
        "so it carries the database's own skew toward densely sequenced species. "
        "16S lane: in-silico PCR cuts 8 regions of each organism's 16S, from V1–V2 to V7–V9 "
        f"with full-length 16S as the control ({n['original']:,}, {n['kmeans']:,} and "
        f"{n['random']:,} amplicons on the stratified, k-means and random panels). Each "
        "amplicon is aligned to all 459,301 unique BV-BRC 16S sequences and keeps a "
        f"{MD.DEPTH}-hit record. Self-exclusion then removes the source genome, its taxon and "
        "every hit at 0.987 identity or above, as if the organism were new to the database, "
        "and the exact-tie union selects from what is left; the anchor's lineage is scored "
        "at genus and family, and the union's gene families against the source's. "
        f"MAG lane: each stratified-panel genome is cut into ~{MD.FRAG_MEAN_KB} kb fragments "
        f"(a {MD.FRAG_MIN_KB} kb floor) and whole fragments are dropped, giving "
        f"{f['n_drafts']:,} drafts at {f['levels']} completeness levels from 100% down to "
        f"{f['lo']:.0%}; every panel is also rerun at {len(lev)} levels ({lev[0]:.0%} to "
        f"{lev[-1]:.0%}). Drafts are matched by skani whole-genome ANI against a sketch of the "
        f"panel itself and, on {min(g.values()):,}–{max(g.values()):,} organisms per panel, "
        f"against the {MD.GTDB_REPS:,} GTDB r226 species representatives. Self-exclusion "
        "removes the source genome, its taxon and every hit at 95 ANI or above, the species "
        "boundary and the "
        "analogue of the 16S lane's 0.987. Both lanes use the same truth and the same "
        "conditional and whole-population metrics, so MAG and 16S can be compared organism by "
        "organism."))
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


def estobj_figures(D):
    """Stage E's E3 against k/n and Stage O's flux contrasts, with the claims c_estobj and
    c_new make checked as the full deck's s06f/s06g check them."""
    P = D["pfe"]
    A, N1 = P["all"], P["n1"]
    est = lambda blk, e, k: blk[e][k]["est"] if isinstance(blk[e][k], dict) else blk[e][k]  # noqa: E731
    fams = ("E0", "E1", "E2", "E3")
    rows = {k: [est(A, e, k) for e in fams] for k in ("union_brier", "union_auroc")}
    O = D["otune"]
    C = O["confirm"]
    ss = C["contrasts"]["selected-shipped"]["metrics"]["taf_face_mid"]
    dec = D["otune_dec"]
    f = {"A": A, "N1": N1, "est": est, "brier": rows["union_brier"],
         "auroc": rows["union_auroc"], "O": O, "C": C, "both": ss["pooled"]["A"],
         "n1": ss["per_stratum"]["n1"], "e3": dec["est@X100-shipped"]["taf_face_mid"]["A"],
         "x1": dec["kn X1-X100"]["taf_face_mid"]["A"],
         "xe": dec["est X1-X100"]["taf_face_mid"]["A"], "sesoi": O["meta"]["sesoi"],
         "hd": C["headroom_confirm"]["A"],
         "deg": C["degeneracy"]["selected-shipped"]["faces.opt.n_var_1e-3"]["pooled"]["A"],
         "shf": C["degeneracy"]["selected-shuffled"]["faces.opt.n_var_1e-3"]["pooled"]["A"]}
    dd = P["decision"]["E3"]
    for ok, claim in (
            (all(P["decision"][e]["beats_k_over_n"] for e in fams[1:]), "every estimator beats k/n"),
            (dd["d_union_brier"]["hi"] < 0, "E3's union Brier contrast against k/n lies below 0"),
            (min(f["brier"]) == f["brier"][3] and max(f["auroc"]) == f["auroc"][3],
             "E3 is best on union Brier and AUROC"),
            (C["selected"] == "est/lin/G0/X1", "the selected objective is est/lin/G0/X1"),
            (C["labels"]["selected-shipped|A"].startswith("lowers truth-absent flux by >= SESOI"),
             "selected - shipped meets the SESOI at n >= 2"),
            (abs(D["otune_extra"]["selected-shipped"] - f["both"]["mean"]) < 1e-9,
             "per-amplicon rows reproduce the pooled contrast"),
            (f["n1"]["ci"][0] < 0 < f["n1"]["ci"][1], "the n = 1 interval spans 0"),
            (abs(f["e3"]["mean"] + f["xe"]["mean"] - f["both"]["mean"]) < 0.0005,
             "E3 alone + exchange weight 1 given E3 = the tuned objective"),
            (f["shf"]["comparator_mean"] < f["deg"]["arm_mean"],
             "shuffled p leaves fewer free reactions")):
        if not ok:
            raise SystemExit(f"{MD.CURRENT_SLIDE}: no longer true: " + claim)
    return f


def c_estobj(prs, D, F):
    """Stage E and Stage O on one slide: E3's better p, and the flux the tuned objective
    removes (the full deck's s06f_estimators and s06g_objective, condensed)."""
    f = estobj_figures(D)
    A, N1, est, C, O = f["A"], f["N1"], f["est"], f["C"], f["O"]
    s = blank(prs)
    title(s, "A better p and a tuned objective cut wrong flux",
          sub="E3: more hits plus a family prior · objective tuned on half the amplicons, "
              "tested on the rest")
    cw = 7.6
    MD.objective_bars(s, M, BODY_TOP + 0.05, cw, 4.1, f["e3"], f["xe"], f["both"], f["n1"],
                      f["sesoi"])
    rx = M + cw + 0.4
    rw = W - M - rx
    b, a = f["brier"], f["auroc"]
    c = stat(s, rx, BODY_TOP + 0.35, rw, f"{rd(b[0])} → {rd(b[3])}",
             "Brier on the union's ECs, k/n → E3 (lower better)", vsize=28, h=1.5)
    arrow_value(c)
    c = stat(s, rx, BODY_TOP + 2.15, rw, f"{rd(a[0])} → {rd(a[3])}",
             "AUROC: E3 ranks the union's ECs, k/n barely", vsize=28, h=1.5)
    arrow_value(c)
    tf = takeaway(s, "Tuned objective: past the SESOI at n ≥ 2, not at n = 1", M, 5.75,
                  W - 2 * M)
    tf.paragraphs[0].alignment = PP_ALIGN.CENTER
    footer(s, f"Controlled genus-blocked benchmark · declared post hoc · flux: confirm half, "
              f"{C['n_rows']} amplicons held out from selection · pooled n ≥ 2 unless marked")
    pts = lambda v: sg(100 * v, 2)                                                       # noqa: E731
    ci = lambda d: f"[{pts(d['ci'][0])}, {pts(d['ci'][1])}]"                            # noqa: E731
    hd, deg, shf = f["hd"], f["deg"], f["shf"]
    wu, w2 = MD.PFE_WITHIN_UNION_F1["all"], MD.PFE_WITHIN_UNION_F1["n2+"]
    notes(s, (
        "Two post-hoc studies, each a slide in the full deck. Stage E asked whether the "
        "per-function weight p can be better than k/n, the share of the selected genomes that "
        "carry a function (the calibration slide's subject). E3 takes identity-weighted "
        "prevalence over the first 100 hits of the record, not only the exact-tie set, shrunk "
        "toward the anchor's family prior (τ = 0.02, α = 2). Scored by nested cross-validation "
        f"on the controlled genus-blocked benchmark ({A['n_cells']:,} cells), it lowers the "
        f"Brier score on the union's ECs from {rd(b[0])} for k/n to {rd(b[3])} and raises their "
        f"AUROC from {rd(a[0])} to {rd(a[3])}; its universe Brier over all ECs is "
        f"{rd(est(A, 'E3', 'universe_brier'), 4)} against {rd(est(A, 'E0', 'universe_brier'), 4)}. "
        "About half of the union Brier gain is recalibration k/n could have without new "
        f"information (isotonic k/n alone reaches {rd(est(A, 'E4(E0)', 'union_brier'), 4)}). "
        f"Where prFBA usually runs, one selected genome ({N1['n_cells']:,} of "
        f"{A['n_cells']:,} cells), k/n gives every EC p = 1 and cannot rank (AUROC "
        f"{rd(est(N1, 'E0', 'union_auroc'), 1)}); E3 ranks them at "
        f"{rd(est(N1, 'E3', 'union_auroc'), 3)}. Thresholded, the gain within the union's ECs is "
        f"small: macro F1 {sg(wu[0], 4)} [{rd(wu[1], 4)}, {rd(wu[2], 4)}], below the 0.02 "
        f"SESOI, and {sg(w2[0], 4)} at n ≥ 2. Stage O then tuned the objective on half of the "
        f"weighted-objective test's amplicons ({O['tune']['n_rows']}) and tested it on the other "
        f"half ({C['n_rows']}): it selected E3's p, the shipped linear cost, gapfilled "
        "reactions at p = 0 and exchange weight 1. On the confirm half, pooling n ≥ 2, it "
        f"changes the truth-absent flux at the face midpoint by {pts(f['both']['mean'])} points "
        f"of the two-tier objective's gene-backed flux {ci(f['both'])} against the shipped "
        f"objective, past the SESOI of {100 * f['sesoi']:.0f} point (chart, teal; the dashed "
        "lines mark whole points). Exploratory decomposition: E3's p alone, at the shipped "
        f"exchange weight, gives {pts(f['e3']['mean'])} {ci(f['e3'])}, and lowering the "
        f"exchange weight to 1, given E3's p, adds {pts(f['xe']['mean'])} {ci(f['xe'])}, so the "
        "two add up to the whole, as in the manuscript (on k/n, exchange weight 1 alone gives "
        f"{pts(f['x1']['mean'])} {ci(f['x1'])}). For one-genome selections it is not "
        f"established: {pts(f['n1']['mean'])} {ci(f['n1'])} on {f['n1']['n']} rows (amber, "
        "hollow). "
        f"The tuned objective captures {rd(100 * hd['selected_own_cell'], 0)}% of an oracle "
        f"weighting's headroom on the confirm half, the shipped one {rd(100 * hd['shipped'], 0)}% "
        "there; it leaves "
        f"{rd(deg['arm_mean'], 1)} gene-backed reactions free on the optimal face against "
        f"{rd(deg['comparator_mean'], 1)} shipped, but shuffled costs leave "
        f"{rd(shf['comparator_mean'], 1)}, so fewer ties is not evidence of information. "
        "Caveats: both studies are post hoc; the benchmark reference is about 7,000 genomes "
        "per fold against about 315,000 in production; the lower exchange weight should be "
        "tested on measured exchange data, so ModelSEEDpy's default stays at 100 and prFBA ships it "
        "only as its recommended call. " + MD.integration_sentence(D) +
        " Sources: "
        "region_validation/shared_reference/PFUNC_ESTIMATORS_RESULTS.md and "
        "data/pfunc_estimators/scores.json; OBJECTIVE_TUNING_RESULTS.md, "
        "data/objective_tuning/summary.json, verification/out_c8.json and "
        "per_amplicon.jsonl.gz; both re-derived in their verification/VERIFICATION.md."))
    return s


def tmpl7_figures(D):
    """The seven-reaction template fix's headline counts, as the full deck's s17e reads them
    (the ModelSEEDTemplates-draftloops validation when present, else TMPL7_FALLBACK)."""
    t7 = D["tmpl7"]
    if t7:
        sc = t7["scan"]["by_templates"]
        return {"drafts": sc["21a8084"]["drafts"], "leak_before": sc["21a8084"]["leaking"],
                "leak_after": sc["tmpl7"]["leaking"], "rows": t7["counts"]["rows"],
                "gated": t7["counts"]["gated_after"],
                "gapfilled": t7["counts"]["rows_gapfilled"],
                "growing": t7["counts"]["growing_after"]}
    return dict(MD.TMPL7_FALLBACK)


def c_new(prs, D, F):
    """The brief variants' one slide on what is new since the deck of 2026-11-01: a result
    and its headline number per row, each a slide of the full update deck."""
    T = D["taxd"]["U100"]["regions"]
    regs = list(MD.REG_LABEL)
    m7, m5 = "exclude_0.987", "exclude_0.995"
    ab = {m: [sum(T[r]["abstention"][m]["all"]["all"][k] for r in regs)
              for k in ("n_abstained", "n")] for m in (m7, m5)}
    pc = {m: 100 * ab[m][0] / ab[m][1] for m in (m7, m5)}
    pop = D["gctau"]["_meta"]["headline_macro_recall_pct_population"]
    dpop = [pop[m5][r] - pop[m7][r] for r in regs]
    f = estobj_figures(D)
    t7 = tmpl7_figures(D)
    kb = MD.KBASE_D100
    mag = D["ssu"]["by_kind"]["MAG"]
    ret1 = next(r for r in D["ret"]["clean"] if r["completeness"] == 1.0)
    hi = D["ssu"]["by_kind_completeness"]["MAG"]["95-100"]
    over = (ret1["usable"] / hi["has_ssu"], ret1["full"] / hi["has_ssu_ge1200"])
    for ok, claim in (
            (ab[m5][0] < ab[m7][0] and ab[m5][1] == ab[m7][1], "abstention falls at 0.995"),
            (all(x > 0 for x in dpop), "population gene capture rises in every region"),
            (t7["gated"] == 0 and t7["growing"] == t7["gapfilled"],
             "the seven-reaction fix gates no row and every gapfilled model grows"),
            (t7["leak_after"] < t7["leak_before"], "fewer drafts leak"),
            (kb["failed"] == 0, "no depth-100 KBase upload failed"),
            (mag["has_ssu"] < 0.3, "under 30% of real MAGs carry any 16S")):
        if not ok:
            raise SystemExit("c_new: no longer true: " + claim)
    s = blank(prs)
    title(s, "New since Nov 1",
          sub="Five results since the last update; details in the full deck")
    rows = [
        ["Result", "Headline number"],
        ["Self-exclusion at 0.995 beside 0.987",
         f"No call (all amplicons) {rd(pc[m7], 1)}% → {rd(pc[m5], 1)}% · gene capture "
         f"+{rd(min(dpop), 1)}–{rd(max(dpop), 1)} pts"],
        ["E3 + tuned objective, shipped",
         f"Brier {rd(f['brier'][0])} → {rd(f['brier'][3])} · wrong flux "
         f"{sg(100 * f['both']['mean'], 2)} pts (n ≥ 2)"],
        ["Seven-reaction template fix",
         f"ATP-leaking drafts {t7['leak_before']:,} → {t7['leak_after']:,} · {t7['gated']} gated"],
        ["Depth-100 builds in KBase",
         f"{kb['codif']['n']:,} + {kb['emilykin']['n']:,} genomes · {kb['failed']} failed"],
        ["PICRUSt2 on MAGs: only via a 16S",
         f"{rd(100 * mag['has_ssu'], 0)}% of real MAGs carry one"],
    ]
    table(s, M, BODY_TOP + 0.15, W - 2 * M, rows, col_w=[4.9, W - 2 * M - 4.9], size=19,
          row_h=0.74, emph_col=1)
    footer(s, "Full update deck: presentations/prFBA_updates_2026-10.pptx · depth-100 records, "
              "exact-tie union")
    kc, ke = kb["codif"], kb["emilykin"]
    hd = f["hd"]
    notes(s, (
        "What has been finished since the deck of 2026-11-01, one row per result; the full "
        "update deck covers each, the first, second, third and fifth on slides of their own "
        "and the KBase upload in its notes and closing slide. First, self-exclusion re-run at 0.995 beside the "
        "headline 0.987 (0.987 is the full-length species boundary; for a short amplicon the "
        "identity that best matches species is nearer 100%, Edgar 2018): "
        f"{ab[m7][0]:,} of {ab[m7][1]:,} amplicons get no taxonomic call at 0.987 and "
        f"{ab[m5][0]:,} at 0.995 (shares of every amplicon; the "
        f"{100 * F['sankey']['exclude']['abstained'] / F['sankey']['scoreable']:.0f}% "
        "abstention elsewhere in this deck counts only the genus-scoreable amplicons of the "
        "genus-call slide), and whole-population gene capture rises by "
        f"{rd(min(dpop), 1)}–{rd(max(dpop), 1)} points in every region "
        "(region_validation_d100/data/taxacc_region_summary.json and "
        "gene_capture_tau_compare.json). Second, Stage E's estimator E3 (the first 100 hits "
        "shrunk toward a family prior) lowers the Brier score on the union's ECs from "
        f"{rd(f['brier'][0])} for k/n to {rd(f['brier'][3])} (AUROC {rd(f['auroc'][0])} to "
        f"{rd(f['auroc'][3])}); about half of that Brier gain is recalibration k/n could have "
        "without new information, since the out-of-fold isotonic map of k/n alone reaches "
        f"{rd(f['est'](f['A'], 'E4(E0)', 'union_brier'))}, so E3 against recalibrated k/n is "
        f"{sg(f['est'](f['A'], 'E3', 'union_brier') - f['est'](f['A'], 'E4(E0)', 'union_brier'))}. "
        "Stage O's objective tuned on it changes truth-absent flux by "
        f"{sg(100 * f['both']['mean'], 2)} points against the shipped objective on held-out "
        f"amplicons with two or more genomes selected [{sg(100 * f['both']['ci'][0], 2)}, "
        f"{sg(100 * f['both']['ci'][1], 2)}], past the SESOI of {100 * f['sesoi']:.0f} point, "
        f"{rd(100 * hd['selected_own_cell'], 0)}% of an oracle's headroom; it is not "
        f"established for one-genome selections ({sg(100 * f['n1']['mean'], 2)}, interval "
        "spanning 0). Both are post hoc (PFUNC_ESTIMATORS_RESULTS.md, "
        "OBJECTIVE_TUNING_RESULTS.md). " + MD.integration_sentence(D) + " Third, seven more "
        "reversible template reactions made forward-only (ModelSEEDTemplates branch "
        "fix-draft-atp-loops) cut the raw drafts that make ATP from nothing from "
        f"{t7['leak_before']:,} to {t7['leak_after']:,} of {t7['drafts']:,}; run through ATP "
        f"correction and gapfilling, {t7['gated']} of {t7['rows']:,} rows fail the C12 gate and "
        f"all {t7['gapfilled']:,} gapfilled models still grow (C12_SCOPE.md section 2b). "
        f"Fourth, both depth-100 deployment builds are in KBase since {kb['date']}: codiffusion "
        f"({kc['n']:,} genomes) as GenomeSet {kc['ref']}, beside the 20-hit build, and EmilyKin "
        f"({ke['n']:,} genomes) as GenomeSet {ke['ref']} in a new narrative, with "
        f"{kb['failed']} failed uploads (deployment_d100/logs/kbase_upload_d100.log). Fifth, "
        "PICRUSt2 cannot take a MAG except through the MAG's own 16S, and across GTDB r232 "
        f"only {rd(100 * mag['has_ssu'], 1)}% of {mag['n']:,} real MAGs carry any 16S "
        f"({rd(100 * mag['has_ssu_ge1200'], 1)}% a full-length one), so on real MAGs PICRUSt2 "
        "could answer for about one in four; complete synthetic MAGs overstate its coverage "
        f"{rd(over[0], 1)}–{rd(over[1], 1)}-fold against real MAGs of 95–100% completeness "
        "(region_validation/picrust2_comparison/README.md, 'PICRUSt2 and MAGs'; "
        "data/mag_gtdb_ssu_summary.json and mag_retention_summary.json)."))
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
        (f"{rd(min(fc), 2)}–{rd(max(fc), 2)}", "family accuracy, novel organisms"),
        (f"{100 * k['exclude']['abstained'] / k['scoreable']:.0f}%",
         "abstain, source removed"),
        (sg(cg["delta_f1_population"]), "prFBA − PICRUSt2, like for like"),
        ("PASS" if c12 == "pass" else "FAIL", "gapfill validity (C12)"),
    ]
    if VARIANT == "panels":
        tregs = list(MD.REG_LABEL)
        cell = lambda r, ep: MD.tp3_cell(D, r, "exclude_0.987", ep)  # noqa: E731
        rmo = "delta_random_minus_original"
        cov = [cell(r, ep)[rmo]["est"] for r in tregs
               for ep in ("genus_coverage", "genecap_coverage")]
        q = max(abs(cell(r, "genus_cond_acc")[k]["est"]) for r in tregs
                for k in (rmo, "delta_kmeans_minus_original"))
        P = D["magp"]
        dlt = [P[f"exclude|c1_x0|{r}|f1"]["conditional"]["mag"]
               - P[f"exclude|c1_x0|{r}|f1"]["conditional"]["amp"] for r in tregs]
        rows = [rows[0],
                (f"{sg(max(cov), 2)} to {sg(min(cov), 2)}",
                 f"random-panel coverage; quality within {rd(q, 2)}"),
                (f"{sg(min(dlt))} to {sg(max(dlt))}", "MAG − 16S F1, both answering"),
                rows[2]]
    card = RGBColor(0x0A, 0x4A, 0x55)
    rh, step = 0.8, 0.94                    # four result rows, then the "Open" strip
    for i, (v, lab) in enumerate(rows):
        y = BODY_TOP + 0.05 + i * step
        rect(s, M, y, W - 2 * M, rh, card)
        tf = textbox(s, M + 0.3, y + 0.04, 3.3, rh - 0.08, anchor=MSO_ANCHOR.MIDDLE)
        para(tf, v, size=30, bold=True, color=AMBER if "PICRUSt2" in lab else MINT, font=HEAD,
             first=True, space_after=0)
        tf = textbox(s, M + 3.8, y + 0.04, W - 2 * M - 4.1, rh - 0.08, anchor=MSO_ANCHOR.MIDDLE)
        para(tf, lab, size=19, color=WHITE, first=True, space_after=0)
    # what is open, as on the full deck's s15_next; the integration item is the shared
    # placeholder (make_deck.INTEGRATION_STATUS), drawn dashed amber and italic
    oy, oh = BODY_TOP + 0.05 + 4 * step + 0.2, 0.78
    chip = rect(s, M, oy, 1.15, oh, AMBER)
    ctf = chip.text_frame
    ctf.margin_left = ctf.margin_right = ctf.margin_top = ctf.margin_bottom = 0
    ctf.vertical_anchor = MSO_ANCHOR.MIDDLE
    para(ctf, "Open", size=18, bold=True, color=INK, font=HEAD, align=PP_ALIGN.CENTER,
         first=True, space_after=0)
    integ, integ_pending = MD.integration_status(D)
    items = [integ] + [short for _, short in MD.OPEN_ITEMS]
    gx, gap = M + 1.15 + 0.2, 0.18
    bw = (W - M - gx - gap * (len(items) - 1)) / len(items)
    for i, it in enumerate(items):
        pending = integ_pending and it == integ
        box = rect(s, gx + i * (bw + gap), oy, bw, oh, card, line=AMBER if pending else None)
        if pending:
            box.line.dash_style = MSO_LINE_DASH_STYLE.DASH
            box.line.width = Pt(2)
        btf = box.text_frame
        btf.margin_left = btf.margin_right = Inches(0.12)
        btf.margin_top = btf.margin_bottom = 0
        btf.vertical_anchor = MSO_ANCHOR.MIDDLE
        para(btf, it, size=14, color=WHITE, align=PP_ALIGN.CENTER, first=True, space_after=0,
             italic=pending, line=0.95)
    footer(s, "github.com/freiburgermsu/prFBA  ·  full update deck: "
              "presentations/prFBA_updates_2026-10.pptx", dark=True)
    kc, ke = MD.KBASE_D100["codif"], MD.KBASE_D100["emilykin"]
    n_real = D["accmap"]["coverage"]["mag_selections"]["real_mags_gtdb_r226"]["n_mags"]
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
        "so the leave-genus-out run is valid; no EC number depended on it."
        + (" In this version the second and third rows replace the abstention and gate rows: "
           "across the three 10,000-organism panels, self-exclusion coverage changes with how "
           "the panel was sampled while the quality of the answers given does not; and where "
           "both answer, a synthetic MAG matched by ANI predicts gene content better than any "
           "16S region, though it answers for far fewer organisms." if VARIANT == "panels"
           else "")
        + " Done since the deck of 2026-11-01: the revised manuscript is pushed to the "
        "prFBA_paper main branch, with the 0.995 results and the Stage E and Stage O sections "
        f"written in; both depth-100 deployment builds are in KBase since "
        f"{MD.KBASE_D100['date']} (codiffusion, {kc['n']:,} genomes, GenomeSet {kc['ref']}, "
        f"beside the 20-hit build; EmilyKin, {ke['n']:,} genomes, GenomeSet {ke['ref']} in a new "
        f"narrative; {MD.KBASE_D100['failed']} failed uploads); self-exclusion is scored at "
        "0.995 beside 0.987; and the collaborators gave permission to use the "
        f"{n_real} real MAGs. Open, the strip along the bottom. "
        # PLACEHOLDER: make_deck.integration_sentence() is what the integration pass updates
        # (as in s15_next)
        "The tuned objective is E3's p, the linear cost, gapfilled reactions at p = 0 and "
        "exchange weight 1. " + MD.integration_sentence(D) + " The ModelSEEDpy fix for "
        "MSProbability.prFBA is pushed on branch fix-prfba-single-model, but its pull request "
        "is not opened yet. The five-reaction template fix and the two curated direction "
        "branches are open upstream as ModelSEED/ModelSEEDTemplates pull requests "
        f"{MD.TEMPLATE_PRS_OPEN} (filed 2026-10-27, awaiting review); the seven-reaction branch "
        "fix-draft-atp-loops is pushed to the fork, and its upstream pull request is not opened "
        "yet. The GTDB accession map stays deferred to future work, because wiring it in "
        "changes which genomes a GTDB hit selects."))
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=None)
    ap.add_argument("--brief", action="store_true", help="the 10-slide version")
    ap.add_argument("--panels", action="store_true",
                    help="the 10-slide version with one PICRUSt2 slide and the panel studies")
    args = ap.parse_args()
    global VARIANT
    VARIANT = "panels" if args.panels else "brief" if args.brief else "condensed"
    if args.out is None:
        args.out = os.path.join(MD.REPO, "presentations", "prFBA_%s_2026-10.pptx" % {
            "panels": "brief_panels", "brief": "brief", "condensed": "condensed"}[VARIANT])
    MD.configure_depth(100, "union", None)
    D = MD.load()
    cm = D["ctl"]["C12"].get("measured") or {}
    D["c12_rescored"] = bool(D["c12fix"] and ("phenotype_gf" + D["amend_tag"] + "/") in
                             (cm.get("binding") or {}).get("phenotype_summary", ""))
    F = json.load(open(os.path.join(FIG, "figure_stats.json")))
    MD.NEW_SLIDES = NEW_BY_VARIANT[VARIANT]     # new since the condensed decks of 2026-11-01
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    plan = [c_title,
            lambda p, D, F: section(p, "Methods", "Align · select · merge · validate"),
            c_workflow, c_arms, c_sankey, c_picrust2,
            lambda p, D, F: section(p, "Results", "What survives when the source is removed"),
            c_heatmap, c_space, MD.s09b_tau995, MD.s11_results, MD.s16_kfold_result,
            MD.s17_phenotype, MD.s06c_calibration, MD.s06d_wobj, c_estobj, MD.s17_crossover,
            MD.s19b_picrust2_mags, MD.s22_three_panel, c_close]
    if VARIANT == "brief":
        plan = [c_title, c_workflow, c_arms, c_picrust2, c_sankey, c_space,
                MD.s16_kfold_result, MD.s17_phenotype, c_new, c_close]
    elif VARIANT == "panels":
        plan = [c_title, c_workflow, c_panels, c_sankey, c_space, MD.s22_three_panel,
                MD.s17_crossover, MD.s16_kfold_result, c_new, c_close]
    for fn in plan:
        MD.CURRENT_SLIDE = getattr(fn, "__name__", "section")
        MD._STATS.clear()
        if fn.__module__ == "make_deck":
            s = fn(prs, D)
            for old, new in NOTE_FIXES.get(fn.__name__, lambda D: [])(D):
                retarget_notes(s, old, new)
        else:
            fn(prs, D, F)
        if MD._STATS:
            size = min(v for _, v in MD._STATS)
            for run, _ in MD._STATS:
                run.font.size = Pt(size)
    prs.core_properties.title = "prFBA: method and validation (%s)" % VARIANT
    prs.core_properties.author = "Andrew Freiburger"
    prs.save(args.out)
    print(f"[condensed] {len(prs.slides._sldIdLst)} slides -> {args.out}")


if __name__ == "__main__":
    main()
