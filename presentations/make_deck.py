#!/usr/bin/env python
"""Build the prFBA update deck from the repositories' own result files.

Every number on a slide is read from the JSON/CSV the pipeline produced --
``taxacc_region_summary.json``, ``gene_capture_by_region.json``,
``comparison_summary.json``, ``annotation_ceiling.json``,
``depth_probe_results.json``, the MAG arm's ``mag_summary_full.json`` /
``crossover.json`` / ``gtdb_summary_2000.json``, and the two-panel comparison's
``panel_comparison.json`` / ``mag_panel_comparison.json`` / ``space_deciles.json``
-- so the deck cannot drift from the data it reports.

The one exception is the block of constants below: a handful of MAG-arm figures
that the pipeline reports only in prose. Each names the file it came from.

    python presentations/make_deck.py --out presentations/prFBA_updates_2026-10.pptx

Record-depth switch: ``--depth 100 [--rule legacy|union] [--constants JSON]`` reads the
amplicon-arm summaries (taxacc/gene-capture, PICRUSt2 comparison, three-panel comparison,
alignment wall time) from the depth-100 rerun trees instead (see
region_validation_d100/RUNBOOK_d100.md) and writes ``prFBA_updates_2026-10_d100.pptx``.
The hand-entered self-recovery / emptied constants are then read from ``--constants``
(default presentations/deck_constants_d100.json) and must all be filled in. Slide PROSE
still describes the depth-20 result; it is rewritten in the deck-regeneration step.
The default (--depth 20) is byte-for-byte the cited deck: every path and constant unchanged.
"""
import argparse
import json
import os

from pptx import Presentation
from pptx.chart.data import CategoryChartData
from pptx.dml.color import RGBColor
from pptx.enum.chart import XL_CHART_TYPE, XL_LEGEND_POSITION, XL_LABEL_POSITION
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.util import Inches, Pt, Emu

REPO = os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
RV = os.path.join(REPO, "region_validation")
CMP = os.path.join(RV, "picrust2_comparison", "data")
MAG = os.path.join(REPO, "mag_validation")
TP = os.path.join(REPO, "two_panel")
SR = os.path.join(RV, "shared_reference")

# --- record-depth switch (see module docstring; --depth 20 keeps every value below) ----
# Only the amplicon-arm inputs that depend on the record depth move; the MAG arm, the
# k-fold benchmark, the depth probe, the annotation ceiling / dialect and the two-panel
# files are depth-independent and keep their sources.
DEPTH = 20
AMP_RV = RV                                                    # "tax" / "cap"
CMP_SUMMARY = os.path.join(CMP, "comparison_summary.json")      # "cmp"
TP3_JSON = os.path.join(REPO, "three_panel", "three_panel_comparison.json")  # "tp3"
ALIGN_STAT = ("4.5 h", "to align all 60,461 unique benchmark amplicons against the full "
              "reference set")                                  # s04 stat (depth-20 run)
DEPTH_CONSTANTS = ("SELF_TOP20", "SELF_SAT", "SELF_CRED", "EMPTIED")
DEPTH20_PROSE = ("s04 'Two stages' text", "s07/s11 'record depth 20' (line ~373)",
                 "s14_depth (probe framing, '20-hit record')", "s16/s18 '20-hit record' text",
                 "s21/s22 deciles text", "s23_reversal 'The flip' / 'Self-inclusion too' / "
                 "'What it licenses' bullets")


def configure_depth(depth, rule="legacy", constants=None):
    """Point the depth-dependent inputs at the depth-100 rerun trees (no-op at depth 20)."""
    global DEPTH, AMP_RV, CMP_SUMMARY, TP3_JSON, ALIGN_STAT
    if depth == 20:
        return
    if depth != 100:
        raise SystemExit("--depth must be 20 or 100")
    DEPTH = depth
    sfx = "_d100" + ("_legacy" if rule == "legacy" else "")
    AMP_RV = os.path.join(REPO, "region_validation" + sfx)
    CMP_SUMMARY = os.path.join(REPO, "region_validation_d100", "picrust2_comparison", "data",
                               "comparison_summary.json")
    TP3_JSON = os.path.join(REPO, "three_panel", "three_panel_comparison%s.json" % sfx)
    ms_path = os.path.join(REPO, "region_validation_d100", "data", "hits", "raw",
                           "method_selection.json")
    if not os.path.exists(ms_path):
        raise SystemExit("--depth %d: no depth-%d alignment yet (%s)" % (depth, depth, ms_path))
    ms = json.load(open(ms_path))
    if ms["params"]["k2"] != depth:
        raise SystemExit("region_validation_d100 alignment has k2=%s, not %d"
                         % (ms["params"]["k2"], depth))
    ALIGN_STAT = ("%.1f h" % (ms["wall_seconds_total"] / 3600.0),
                  "to align all {:,} unique benchmark amplicons at depth {} against the full "
                  "reference set".format(ms["n_asvs"], depth))
    path = constants or os.path.join(REPO, "presentations", "deck_constants_d%d.json" % depth)
    vals = json.load(open(path))
    missing = [k for k in DEPTH_CONSTANTS
               if not isinstance(vals.get(k), dict)
               or any(vals[k].get(p) is None for p in ("original", "kmeans", "random"))]
    if missing:
        raise SystemExit("%s: fill in %s (each {original, kmeans, random}) before building the "
                         "depth-%d deck" % (path, ", ".join(missing), depth))
    for k in DEPTH_CONSTANTS:
        globals()[k] = {p: vals[k][p] for p in ("original", "kmeans", "random")}
    print("[make_deck] DEPTH=%d rule=%s: tax/cap <- %s, cmp <- %s, tp3 <- %s, constants <- %s"
          % (depth, rule, AMP_RV, CMP_SUMMARY, TP3_JSON, path))
    print("[make_deck] WARNING: slide prose still describes the depth-20 result; revise: "
          + "; ".join(DEPTH20_PROSE))

# --- palette: teal/mint for the biology, amber reserved for what is new ------
INK = RGBColor(0x10, 0x29, 0x2E)
DEEP = RGBColor(0x05, 0x34, 0x3D)
TEAL = RGBColor(0x02, 0x80, 0x90)
SEAFOAM = RGBColor(0x00, 0xA8, 0x96)
MINT = RGBColor(0x02, 0xC3, 0x9A)
AMBER = RGBColor(0xE8, 0x87, 0x1A)      # chip fills only
AMBER_TXT = RGBColor(0xC9, 0x6E, 0x0A)  # amber readable on the light tint
SEA_TXT = RGBColor(0x00, 0x7A, 0x6C)    # seafoam readable on the light tint
MUTED = RGBColor(0x5B, 0x71, 0x79)
LIGHT = RGBColor(0xEF, 0xF5, 0xF5)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
HEAD, BODY = "Cambria", "Calibri"

W, H = 13.333, 7.5
M = 0.6                      # left/right margin
TITLE_Y, TITLE_H = 0.42, 0.85
BODY_TOP = 1.5               # every content slide starts its body here


def load():
    j = lambda p: json.load(open(p))  # noqa: E731
    return {
        "tax": j(os.path.join(AMP_RV, "data", "taxacc_region_summary.json")),
        "cap": j(os.path.join(AMP_RV, "data", "gene_capture_by_region.json")),
        "cmp": j(CMP_SUMMARY),
        "ceil": j(os.path.join(CMP, "annotation_ceiling.json")),
        "probe": j(os.path.join(CMP, "depth_probe", "depth_probe_results.json")),
        "dialect": j(os.path.join(CMP, "annotation_dialect.json")),
        # the MAG arm: -n 100 is the audited RESULTS_v2 basis, used everywhere
        # except the cross-panel slide, which is depth-matched at -n 500
        "mag": j(os.path.join(MAG, "out", "mag_summary_full.json")),
        "magx": j(os.path.join(MAG, "out", "crossover.json")),
        "gtdb": j(os.path.join(MAG, "out_gtdb", "gtdb_summary_2000.json")),
        # the two-panel comparison
        "pan": j(os.path.join(TP, "panel_comparison.json")),
        "panmag": j(os.path.join(TP, "mag_panel_comparison.json")),
        "dec": j(os.path.join(TP, "space_deciles.json")),
        "knn": j(os.path.join(TP, "knn_density.json")),
        # the shared-reference k-fold benchmark and its phenotype arm
        "kg": j(os.path.join(SR, "data", "summary_genus.json")),
        "kf": j(os.path.join(SR, "data_family", "summary_family.json")),
        # the phenotype arm, re-run 2026-10-18 with one output directory per pass, so each
        # summary covers exactly its pre-specified selection (1,000 ATP / 200 gapfilled)
        "phatp": j(os.path.join(SR, "data", "phenotype_atp", "phenotype_summary.json")),
        "phgf": j(os.path.join(SR, "data", "phenotype_gf", "phenotype_summary.json")),
        # truth genomes whose gapfilled model fails the sanity gate (C12 leg b), counted
        # as genomes: meta.n_gated_out counts rows, one per (genome, pipeline)
        "phgf_gated": gated_truth_genomes(os.path.join(SR, "data", "phenotype_gf",
                                                       "per_genome.jsonl")),
        # the three-panel comparison: taxonomy-stratified / k-means / uniform random
        "tp3": j(TP3_JSON),
        "tp3g": j(os.path.join(REPO, "three_panel", "three_panel_mag_gtdb.json")),
        "tp3p": j(os.path.join(REPO, "three_panel", "three_panel_mag_panelarm.json")),
        "dec3": j(os.path.join(REPO, "three_panel", "space_deciles_three.json")),
        "sel3": {k: j(p)["genomes"] for k, p in (
            ("original", os.path.join(RV, "data", "selection_10k.json")),
            ("kmeans", os.path.join(REPO, "kmeans_panel", "selection_kmeans.json")),
            ("random", os.path.join(REPO, "random_panel", "selection_random.json")))},
        # why the MAG arm has to exist, and what it costs to build one
        "ssu": j(os.path.join(CMP, "mag_gtdb_ssu_summary.json")),
        "ret": j(os.path.join(CMP, "mag_retention_summary.json")),
        "feed": j(os.path.join(CMP, "mag_direct_feed_summary.json")),
        "hmm": j(os.path.join(CMP, "src16s_hmm_check_summary.json")),
    }


# --- numbers the pipeline reports only in prose; source named on each ---------
# mag_validation/out/RESULTS_v2.md, "the AF gate is almost everything"
# RESULTS_v2.md's table counts qualifying hits over all 9,997 MAGs but divides by the 9,837
# scoreable ones; recounted on one basis (scoreable organisms, c100_x0 hits after the
# scorer's exclusion) on 2026-10-18. It is the share with a qualifying surviving hit, before
# selection, not the resolve rate (0.308).
AF_SWEEP = [
    ("as shipped — ANI \u2265 83, AF \u2265 15%", 0.311),
    ("ANI floor 83 \u2192 0", 0.316),
    ("AF gate 15% \u2192 5%", 0.379),
    ("AF gate 15% \u2192 2%", 0.452),
    ("both relaxed — ANI \u2265 80, AF \u2265 2%", 0.567),
    ("any surviving hit at all (ceiling)", 0.655),
]
# GTDB r226's skani sketch holds 143,614 species representatives; 732,475 is the accession
# count of the bac120 + ar53 TAXONOMY files, which RESULTS_v2.md's "73x" divided by mistake.
GTDB_FOLD = 14                 # 143,614 / 10,000-genome panel (14.4x), verified 2026-10-18
PANEL_KM_SPECIES = 4935        # two_panel/RESULTS.md section 1
PANEL_KM_PHYLA = 142           # two_panel/RESULTS.md section 1
PANEL_OG_PHYLA = 206           # two_panel/RESULTS.md section 1
PANEL_OVERLAP = 804            # TWO_PANEL_RUNBOOK.md part 2
# mag_validation/synth_mags.py: module constants and the docstring's justification
FRAG_MEAN_KB, FRAG_MIN_KB = 50, 2.5
PANEL_CLOSED_PCT, PANEL_LT10_PCT = 19.8, 40.4
# --- the three-panel analysis (prFBA_paper sec:twopanel; verified 2026-10-18) ---------
# self-recovery, amplicon level, from each panel's data/hits/asv_top20_alignment_hits.json
# (676 MB per panel, so not loaded here): source in its own top 20; misses where the top
# 20 are all exact ties; misses where the source's sequence is credited to another genome
SELF_TOP20 = {"original": 96.7, "kmeans": 69.0, "random": 61.2}
SELF_SAT = {"original": 3.3, "kmeans": 15.6, "random": 24.6}
SELF_CRED = {"original": 0.0, "kmeans": 15.4, "random": 14.2}
EMPTIED = {"original": 14.3, "kmeans": 40.4, "random": 49.9}   # full-length, tau = 0.987
GTDB_BEST_ANI = {"original": 89.5, "kmeans": 91.4, "random": 92.1}
# --- the synthetic-genome re-evaluation (prFBA_paper sec:recon) ----------------------
BENCH_N1, DEPLOY_N1 = 86.8, 62.5     # % of answered queries whose selection is one genome
XPHYLUM_DEPLOY = (31.7, 59.0)        # % of n=2 / n>=3 deployment merges spanning phyla
XPHYLUM_BENCH = 6.4                  # same, controlled benchmark (all multi-genome sets)
# mean per-genome EC-F1 difference, prFBA exact-tie minus PICRUSt2, over the 755 paired
# ATP-informative genomes (phenotype_atp rows, deduplicated; recomputed 2026-10-18)
PHENO_EC_GAP = -0.0284
# --- the claim-level audits of the manuscript (prFBA_paper f7fdc38, d2d7516) ----------
# the second audit's record holds 48 flags, one overturned on recheck ("no correction needed")
AUDIT_CLAIMS, AUDIT_DEFECTS = (417, 178), (114, 47)
# the verification of this deck and the manuscript patch, 2026-10-18: five auditors, each
# flag re-derived by a skeptic (claims checked, defects confirmed)
DECK_AUDIT = (598, 54)
# the second pass, over the fixes themselves (claims checked, distinct defects confirmed)
FIX_AUDIT = (377, 14)
# faults our own runs and reproduction checks caught (the phenotype gating and duplicate-row
# faults were claim-audit catches, so they are listed with those, not here)
PIPELINE_CAUGHT = [
    "A gene-family fetcher with no User-Agent: every request refused, fds leaked",
    "A three-panel scorer that counted any non-abstention as a call at every rank "
    "(genus coverage overstated by up to 0.19)",
    "A density-decile comparison that read its stored array in reverse",
    "A reproduction check that passed while comparing zero cells (a wrong key path)"]


def sg(x, d=3):
    """Signed number with a typographic minus, so a line never breaks inside it.

    Rounds half away from zero on the stored decimal (-0.1235 -> -0.124), as the
    manuscript does; binary float formatting would give -0.123."""
    from decimal import Decimal, ROUND_HALF_UP
    v = Decimal(repr(x)).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP)
    return f"{v:+.{d}f}".replace("-", "\u2212")


def rd(x, d=3):
    """Unsigned counterpart of sg(): half away from zero on the stored decimal."""
    from decimal import Decimal, ROUND_HALF_UP
    return f"{Decimal(repr(x)).quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP):.{d}f}"


def gated_truth_genomes(path):
    out = set()
    with open(path) as fh:
        for line in fh:
            if '"truth"' in line and '"gated"' in line:
                r = json.loads(line)
                if r["arm"] == "truth" and r["status"] == "gated":
                    out.add(r["src_genome_id"])
    return sorted(out)


def composition(gs):
    """Panel composition, counting unclassified placeholders as their own genus bins."""
    import collections
    named = lambda v: bool(v) and not str(v).startswith("unclassified")   # noqa: E731
    gen = collections.Counter(e["genus"] for e in gs)
    v = sorted(gen.values()); n, tot = len(v), sum(v)
    gini = (2 * sum((i + 1) * x for i, x in enumerate(v))) / (n * tot) - (n + 1) / n
    arc = sum(1 for e in gs if e["domain"] == "Archaea")
    return dict(n=len(gs), archaea=arc, arc_pct=100 * arc / len(gs),
                phyla=len({e["phylum"] for e in gs}),
                families=len({e["family"] for e in gs if named(e["family"])}),
                genera=len({e["genus"] for e in gs if named(e["genus"])}),
                species=len({e["species"] for e in gs}), gini=gini, largest=max(v))


def tp3_cell(D, region, regime, endpoint):
    return next(c for c in D["tp3"]["cells"] if c["region"] == region
                and c["regime"] == regime and c["endpoint"] == endpoint)


REG_LABEL = {"FullLength16S": "Full-length", "V1-V2": "V1–V2", "V1-V3": "V1–V3",
             "V3-V4": "V3–V4", "V4": "V4", "V4-V5": "V4–V5", "V4-V5_944R": "V4–V5 944R",
             "V7-V9": "V7–V9"}
PANEL_COLORS = None   # set after the palette, see below


# --------------------------------------------------------------------- helpers
def textbox(slide, x, y, w, h, *, anchor=MSO_ANCHOR.TOP):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    tf.vertical_anchor = anchor
    return tf


def para(tf, text, *, size=15, bold=False, color=INK, font=BODY, space_after=6,
         align=PP_ALIGN.LEFT, italic=False, first=False, line=None):
    p = tf.paragraphs[0] if first else tf.add_paragraph()
    p.alignment = align
    p.space_after = Pt(space_after)
    if line:
        p.line_spacing = line
    r = p.add_run()
    r.text = text
    f = r.font
    f.size, f.bold, f.italic, f.name = Pt(size), bold, italic, font
    f.color.rgb = color
    return p


def rect(slide, x, y, w, h, fill, *, shape=MSO_SHAPE.ROUNDED_RECTANGLE, line=None,
         adj=0.06):
    s = slide.shapes.add_shape(shape, Inches(x), Inches(y), Inches(w), Inches(h))
    s.fill.solid()
    s.fill.fore_color.rgb = fill
    if line is None:
        s.line.fill.background()
    else:
        s.line.color.rgb = line
        s.line.width = Pt(1)
    s.shadow.inherit = False
    if shape == MSO_SHAPE.ROUNDED_RECTANGLE:
        try:
            s.adjustments[0] = adj
        except (IndexError, AttributeError):
            pass
    s.text_frame.word_wrap = True
    s.text_frame.vertical_anchor = MSO_ANCHOR.TOP
    return s


def blank(prs):
    return prs.slides.add_slide(prs.slide_layouts[6])


def background(slide, color):
    slide.background.fill.solid()
    slide.background.fill.fore_color.rgb = color


def title(slide, text, *, sub=None, dark=False, tag=None):
    """Content-slide title: identical position, size and alignment on every slide."""
    if tag:
        chip = rect(slide, W - M - 0.72, TITLE_Y + 0.05, 0.72, 0.32, AMBER)
        tf = chip.text_frame
        tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
        tf.vertical_anchor = MSO_ANCHOR.MIDDLE
        para(tf, tag, size=11.5, bold=True, color=INK, align=PP_ALIGN.CENTER, first=True)
    tf = textbox(slide, M, TITLE_Y, W - 2 * M - (0.95 if tag else 0), TITLE_H)
    para(tf, text, size=30, bold=True, color=WHITE if dark else INK, font=HEAD, first=True,
         space_after=2)
    if sub:
        para(tf, sub, size=13.5, color=MINT if dark else MUTED, space_after=0)


def stat(slide, x, y, w, value, label, *, color=TEAL, vsize=None, h=1.45, tint=LIGHT):
    # one rule for every stat card: 34pt, dropping to 28pt for long values
    if vsize is None:
        vsize = 34 if len(value) <= 7 else 28
    card = rect(slide, x, y, w, h, tint)
    tf = card.text_frame
    tf.margin_left = tf.margin_right = Inches(0.14)
    tf.margin_top = Inches(0.12)
    tf.margin_bottom = Inches(0.12)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    para(tf, value, size=vsize, bold=True, color=color, font=HEAD, first=True, space_after=1)
    para(tf, label, size=11.5, color=MUTED, space_after=0, line=0.92)
    return card


def bullets(tf, items, *, size=14.5, space=9, color=INK, head_color=TEAL):
    for i, it in enumerate(items):
        if isinstance(it, tuple):
            head, rest = it
            p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
            p.space_after = Pt(space)
            r = p.add_run()
            r.text = head
            r.font.size, r.font.bold, r.font.name = Pt(size), True, BODY
            r.font.color.rgb = head_color
            r2 = p.add_run()
            r2.text = rest
            r2.font.size, r2.font.name = Pt(size), BODY
            r2.font.color.rgb = color
        else:
            para(tf, it, size=size, color=color, space_after=space, first=(i == 0))


def hanging(p, indent=0.17):
    """Wrap wrapped lines under the text, not under the bullet glyph."""
    pPr = p._p.get_or_add_pPr()
    pPr.set("marL", str(int(indent * 914400)))
    pPr.set("indent", str(-int(indent * 914400)))
    return p


def footer(slide, text, *, dark=False):
    tf = textbox(slide, M, H - 0.52, W - 2 * M, 0.3)
    para(tf, text, size=10, color=MINT if dark else MUTED, first=True, space_after=0)


def style_chart(chart, *, legend=True, size=11):
    chart.has_title = False          # otherwise a single series auto-titles itself
    chart.font.size = Pt(size)
    chart.font.name = BODY
    chart.font.color.rgb = MUTED
    if legend:
        chart.has_legend = True
        chart.legend.position = XL_LEGEND_POSITION.BOTTOM
        chart.legend.include_in_layout = False
        chart.legend.font.size = Pt(11)
        chart.legend.font.color.rgb = INK
    else:
        chart.has_legend = False


def series_color(series, color):
    series.format.fill.solid()
    series.format.fill.fore_color.rgb = color


# ------------------------------------------------------------------- the slides
def s01_title(prs, D):
    s = blank(prs)
    background(s, DEEP)
    tf = textbox(s, M, 2.0, 11.4, 2.4)
    para(tf, "prFBA", size=58, bold=True, color=WHITE, font=HEAD, first=True, space_after=4)
    para(tf, "Probabilistic flux balance analysis for amplicon-resolved metabolic modeling",
         size=22, color=MINT, font=HEAD, space_after=16)
    para(tf, "Project update — a controlled comparison that reverses our own ranking at "
             "the gene level, a phenotype assay that barely separates the methods, and a "
             "third panel that shows which results replicate",
         size=15, color=WHITE)
    tf2 = textbox(s, M, 6.3, 11.4, 0.7)
    para(tf2, "Andrew Freiburger  ·  Tyo & Henry labs  ·  October 2026",
         size=13, color=MINT, first=True, space_after=0)
    s.notes_slide.notes_text_frame.text = (
        "The headline since the last update: the controlled benchmark the dialect slide "
        "called for reverses our earlier result. PICRUSt2 is more accurate at EC level, "
        "beyond the effect size we fixed in advance, at both genus and family rank; "
        "prFBA forced to answer only ties copying the nearest reference. At the "
        "phenotype level the mean ATP score shows no difference at n = 755 — a null, not "
        "an equivalence — though PICRUSt2 is ahead on two of 22 media. A third, random "
        "panel shows the headline self-exclusion coverage is a property of the stratified "
        "panel at record depth 20, and two claim-level audits found and corrected 161 "
        "defects in the manuscript.")
    return s


def s02_problem(prs, D):
    s = blank(prs)
    title(s, "Amplicons name the members, not their metabolism",
          sub="16S surveys capture more species per unit of DNA than shotgun, but carry no genes")
    tf = textbox(s, M, BODY_TOP, 6.5, 4.4)
    bullets(tf, [
        # wording below is the author's own edit of the deck (2026-10-18), kept verbatim
        ("The gap.  ", "16S tells you who is present but GEM modeling needs gene content."),
        ("The usual fixes.  ", "MAGs are incomplete, expensive, and biased to abundant "
         "members. Amplicons are cheap but functionally mute."),
        ("Our approach.  ", "Map each ASV to reference genome(s) via alignment scoring, keep "
         "co-optimal candidates, and preserve ambiguity as probabilities that are sorted out "
         "in optimization."),
    ], size=15, space=14)
    stat(s, 7.42, BODY_TOP, 2.5, "1.43M", "BV-BRC 16S records,\n460,117 unique sequences")
    stat(s, 10.23, BODY_TOP, 2.5, "67,991", "benchmark amplicons\nfrom 10,000 genomes")
    stat(s, 7.42, BODY_TOP + 1.75, 2.5, "3,224", "synthetic genomes\ndeposited in KBase")
    stat(s, 10.23, BODY_TOP + 1.75, 2.5, "8", "16S regions compared\nhead to head",
         color=SEA_TXT)
    footer(s, "Probability-weighted genomes feed a probabilistic FBA objective")
    s.notes_slide.notes_text_frame.text = (
        "The uncertainty in mapping an amplicon to a genome is real information. prFBA "
        "represents it instead of discarding it.")
    return s


def s03_pipeline(prs, D):
    s = blank(prs)
    title(s, "The pipeline, end to end",
          sub="Every selection is auditable against a published identity threshold")
    steps = [
        ("1", "Amplicon", "ASV or in-silico\namplicon sequence"),
        ("2", "Alignment", "edlib top-500 prefilter,\nGPU Smith-Waterman"),
        ("3", "Selection", "exact-tie union,\nspecies de-duplication"),
        ("4", "Synthetic genome", "union of gene content,\nper-function probability"),
        ("5", "prFBA", "probabilities constrain\nthe flux objective"),
    ]
    x, w, gap = M, 2.186, 0.3
    for i, (n, head, desc) in enumerate(steps):
        card = rect(s, x, BODY_TOP + 0.35, w, 2.5, LIGHT)
        tf = card.text_frame
        tf.margin_left = tf.margin_right = Inches(0.16)
        tf.margin_top = Inches(0.78)
        para(tf, head, size=16, bold=True, color=INK, font=HEAD, first=True, space_after=7)
        para(tf, desc, size=12.5, color=MUTED, space_after=0, line=1.05)
        dot = rect(s, x + 0.16, BODY_TOP + 0.08, 0.54, 0.54, TEAL, shape=MSO_SHAPE.OVAL)
        dtf = dot.text_frame
        dtf.margin_left = dtf.margin_right = dtf.margin_top = dtf.margin_bottom = 0
        dtf.vertical_anchor = MSO_ANCHOR.MIDDLE
        para(dtf, n, size=16, bold=True, color=WHITE, font=HEAD,
             align=PP_ALIGN.CENTER, first=True, space_after=0)
        x += w + gap
    tf = textbox(s, M, BODY_TOP + 3.25, W - 2 * M, 1.5)
    bullets(tf, [
        ("Hardware-adaptive.  ", "The prefilter is hardware-independent; only the "
         "scoring backend changes (CUDA kernel, or Biopython on CPU), so outputs are "
         "comparable across machines."),
        ("Decision-tracked.  ", "Every hit carries its disposition and reason, so a "
         "selection can be audited rather than trusted."),
    ], size=14, space=9)
    footer(s, "align_hits.py → select_references.py → build_synthetic_genomes.py → KBase")
    return s


def s04_alignment(prs, D):
    s = blank(prs)
    title(s, "Exact alignment at database scale",
          sub="A classical k-mer retriever beat dense embeddings, so retrieval is exact, not approximate")
    for i, (v, lab, col) in enumerate([
            ("57.3", "GCUPS on one RTX 5070 Ti\n(GPU Smith-Waterman)", TEAL),
            ("7.8×", "faster than a 60-process\nBiopython CPU baseline", TEAL),
            ("400/400", "control pairs bit-exact\nto Biopython", SEA_TXT),
            ("459,301", "unique references\nscanned per amplicon", SEA_TXT)]):
        stat(s, M + i * 3.108, BODY_TOP, 2.808, v, lab, color=col)
    tf = textbox(s, M, BODY_TOP + 1.9, 7.6, 3.0)
    bullets(tf, [
        ("Two stages.  ", "An edlib edit-distance scan keeps the 500 closest references "
         "per amplicon; only those are scored by Smith-Waterman. The 500-deep shortlist "
         "is hardware-independent, so every machine sees the same candidates."),
        ("Why not embeddings.  ", "A nucleotide-transformer space put the true best hit "
         "first only 50.5% of the time; a k-mer retriever reached 75.5%. The 16S cosine "
         "space is too concentrated — random pairs already average 0.90."),
        ("Packaged.  ", "The CUDA kernel is now the standalone gpusw package, so the "
         "aligner is reusable outside this project."),
    ], size=14, space=11)
    stat(s, 8.5, BODY_TOP + 2.0, 4.23, ALIGN_STAT[0], ALIGN_STAT[1], h=1.55)
    footer(s, "Scoring scheme identical on every backend: match +2 / mismatch −3 / gap −5, −2")
    return s


def s05_selection(prs, D):
    s = blank(prs)
    title(s, "Selection keeps the ambiguity it finds",
          sub="Near-identical 16S genes cannot distinguish genomes — so do not pretend otherwise")
    tf = textbox(s, M, BODY_TOP, 7.3, 3.6)
    bullets(tf, [
        ("Exact-tie union.  ", "Ties for the best identity are kept and collapsed to one "
         "genome per BV-BRC species."),
        ("Gene-richest representative.  ", "Among co-optimal conspecifics the most "
         "complete genome represents the species — without letting gene count change "
         "the taxonomic call."),
        ("Thresholds, not vibes.  ", "Species ≥ 0.987, genus ≥ 0.945, family "
         "≥ 0.865; below the family floor the pipeline abstains rather than guess."),
    ], size=14.5, space=12)
    card = rect(s, 8.2, BODY_TOP, 4.53, 2.45, LIGHT)
    tf2 = card.text_frame
    tf2.margin_left = tf2.margin_right = Inches(0.2)
    tf2.margin_top = Inches(0.18)
    para(tf2, "Self-recovery", size=13, bold=True, color=INK, font=HEAD, first=True, space_after=10)
    para(tf2, "Does the selector keep the source organism when it is among the 20 retained hits?",
         size=12, color=MUTED, space_after=12, line=1.05)
    p = tf2.add_paragraph()
    p.alignment = PP_ALIGN.LEFT
    for text, col, sz in (("82.9%", MUTED, 28), ("  →  ", MUTED, 18), ("99.2%", TEAL, 28)):
        r = p.add_run()
        r.text = text
        r.font.size, r.font.bold, r.font.name = Pt(sz), True, HEAD
        r.font.color.rgb = col
    p.space_after = Pt(4)
    para(tf2, "ordered reducer → exact-tie union (99.9% on full-length 16S)",
         size=11, color=MUTED, space_after=0)
    footer(s, "Tie clusters reach ~200 genomes at identical score; the union records them as equivalents")
    return s


def s06_synthetic(prs, D):
    s = blank(prs)
    title(s, "From selected genomes to a probabilistic model",
          sub="Each function carries the fraction of selected references that encode it")
    card = rect(s, M, BODY_TOP, 6.2, 1.5, DEEP)
    tf = card.text_frame
    tf.margin_left = tf.margin_right = Inches(0.25)
    tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    para(tf, "p(function) = genomes encoding it ÷ genomes selected",
         size=15, bold=True, color=WHITE, font=HEAD, first=True, space_after=6)
    para(tf, "core → p = 1.0     accessory → p < 1.0     one reference → p = 1 for all",
         size=13, color=MINT, space_after=0)
    tf2 = textbox(s, M, BODY_TOP + 1.85, 6.2, 2.6)
    bullets(tf2, [
        ("Into the objective.  ", "Each reaction takes the largest weight among its genes, "
         "and prFBA penalizes flux through low-weight reactions. The weight is an "
         "annotation prevalence, not a calibrated probability — next slide."),
        ("Reproducible by construction.  ", "Each genome ships with provenance: git "
         "commit, input hashes, per-stage parameters, and the abstention reason when no "
         "model was built."),
    ], size=14, space=11)
    stat(s, 7.3, BODY_TOP, 2.6, "3,224", "genomes + GenomeSet\nin KBase workspace 261948")
    stat(s, 10.2, BODY_TOP, 2.53, "52", "of 3,276 ASVs abstained\nrather than guess",
         color=AMBER_TXT)
    card2 = rect(s, 7.3, BODY_TOP + 1.85, 5.43, 2.55, LIGHT)
    tf3 = card2.text_frame
    tf3.margin_left = tf3.margin_right = Inches(0.22)
    tf3.margin_top = Inches(0.2)
    para(tf3, "Also extended to whole genomes", size=14, bold=True, color=INK, font=HEAD,
         first=True, space_after=9)
    para(tf3, "Real MAGs enter the same selection path by skani average nucleotide "
              "identity instead of 16S identity — 291 dereplicated MAGs against GTDB R226 "
              "and 2,160 BV-BRC sludge genomes, with an alignment-fraction guard so a "
              "high-ANI, near-zero-coverage hit cannot win. Only BV-BRC matches can go on "
              "to reconstruction. That path now has a validation arm of its own.",
         size=12.5, color=MUTED, space_after=0, line=1.08)
    footer(s, "ANI classes use thresholds recalibrated for genomes: ≥95% same species, 85–95% putatively species-distinct")
    return s


def s07_validation(prs, D):
    s = blank(prs)
    title(s, "A 10,000-organism in-silico validation",
          sub="Run the unaltered pipeline on amplicons whose true genome and gene set we know")
    tf = textbox(s, M, BODY_TOP, 6.1, 3.4)
    bullets(tf, [
        ("The design.  ", "10,000 taxonomy-stratified BV-BRC genomes, deterministic "
         "(seed 1729), in-silico PCR for eight regions — V1–V2 through V7–V9 "
         "with full-length 16S as positive control — giving 67,991 ASVs."),
        ("Scored two ways.  ", "Self-inclusion leaves the source organism in the database "
         "(the realistic upper bound). Self-exclusion removes it and every ≥0.987 "
         "neighbour, then re-selects — performance as if the organism were novel."),
        ("Two measures.  ", "Taxonomic accuracy against the NCBI lineage, and PGFam gene "
         "capture against the source genome's actual gene set."),
    ], size=14, space=12)
    stat(s, 7.1, BODY_TOP, 2.7, "0.94–0.99", "genus accuracy when the\nsource is present")
    stat(s, 10.03, BODY_TOP, 2.7, "95–99%", "of its metabolism\nrecovered")
    stat(s, 7.1, BODY_TOP + 1.75, 2.7, "0.60–0.73", "genus accuracy when removed,\namong amplicons that answer",
         color=AMBER_TXT)
    stat(s, 10.03, BODY_TOP + 1.75, 2.7, "68–72%", "gene capture when removed,\namong amplicons that answer", color=AMBER_TXT)
    tf2 = textbox(s, 7.1, BODY_TOP + 3.55, 5.63, 0.9)
    para(tf2, "Family-level assignment is far more robust than genus: 0.79–0.86 even "
              "for novel organisms — family-level metabolism transfers, genus-level is "
              "where short amplicons pay.", size=12.5, color=MUTED, first=True, space_after=0,
         line=1.08)
    footer(s, "All eight regions differ significantly at genus under self-exclusion (paired McNemar, Holm p < 0.05)")
    return s


def s08_regions(prs, D):
    s = blank(prs)
    title(s, "Longer reads buy taxonomic resolution, not raw gene recovery",
          sub="Genus accuracy by region, source organism present vs removed")
    regions = ["FullLength16S", "V1-V2", "V1-V3", "V7-V9", "V3-V4", "V4-V5_944R", "V4-V5", "V4"]
    tax = D["tax"]["regions"]
    inc = [tax[r]["anchor"]["include"]["all"]["all"]["per_rank"]["Genus"]["correct"] for r in regions]
    exc = [tax[r]["anchor"]["exclude_0.987"]["all"]["all"]["per_rank"]["Genus"]["correct"] for r in regions]
    cd = CategoryChartData()
    cd.categories = [r.replace("FullLength16S", "Full-length").replace("_944R", " 944R")
                     .replace("V1-V2", "V1–V2").replace("V1-V3", "V1–V3")
                     .replace("V3-V4", "V3–V4").replace("V4-V5", "V4–V5")
                     .replace("V7-V9", "V7–V9") for r in regions]
    cd.add_series("Source present", inc, number_format="0.00")
    cd.add_series("Source removed", exc, number_format="0.00")
    gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(BODY_TOP),
                            Inches(8.1), Inches(4.4), cd)
    ch = gf.chart
    style_chart(ch)
    series_color(ch.series[0], TEAL)
    series_color(ch.series[1], AMBER)
    ch.value_axis.maximum_scale = 1.08
    ch.value_axis.minimum_scale = 0.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    for sr in ch.series:
        sr.has_data_labels = True
        dl = sr.data_labels
        dl.show_value = True         # has_data_labels alone never sets showVal
        dl.number_format, dl.number_format_is_linked = "0.00", False
        dl.position = XL_LABEL_POSITION.OUTSIDE_END
        dl.font.size = Pt(9)
        dl.font.color.rgb = MUTED
        dl.font.name = BODY
    tf = textbox(s, 9.0, BODY_TOP + 0.15, 3.73, 4.2)
    bullets(tf, [
        ("The gradient.  ", "Full-length 16S and the long V1–V3 generalize best; "
         "V4 — the most commonly sequenced region — generalizes worst."),
        ("Asymmetric cost.  ", "The gene-capture gap is nearly uniform at ~27 points "
         "across every region, while the genus-accuracy gap widens for short amplicons."),
        ("Reading it.  ", "Read length costs taxonomic resolution more than it costs "
         "functional recovery."),
    ], size=13, space=11)
    footer(s, "Anchor predictor, Bacteria + Archaea pooled; conditional on the pipeline making a call")
    return s


def s09_denominators(prs, D):
    s = blank(prs)
    title(s, "Correcting how we scored ourselves",
          sub="Two scoring bugs made the published generalization numbers wrong in both directions")
    tf = textbox(s, M, BODY_TOP, 6.6, 3.6)
    bullets(tf, [
        ("The bug.  ", "The scorer could not tell “the pipeline made no call” from "
         "“the truth has no name at this rank.” 22–25% of source organisms "
         "have no named NCBI genus, and every one was counted as a prediction failure."),
        ("The fix.  ", "Unscoreable ranks now leave the denominator entirely, as the "
         "design spec always specified. Conditional accuracy, correct counts and called "
         "counts are unchanged — verified cell by cell."),
        ("Gene capture.  ", "The mirror-image bug was fixed earlier: abstentions had been "
         "dropped from the mean instead of scoring zero, which inflated capture."),
    ], size=14, space=12)
    card = rect(s, 7.5, BODY_TOP, 5.23, 2.55, LIGHT)
    tf2 = card.text_frame
    tf2.margin_left = tf2.margin_right = Inches(0.22)
    tf2.margin_top = Inches(0.18)
    para(tf2, "Whole-population genus accuracy", size=13.5, bold=True, color=INK, font=HEAD,
         first=True, space_after=12)
    for label, old, new, col in (("novel organisms", "0.23–0.37", "0.30–0.48", TEAL),
                                 ("source present", "0.66–0.76", "0.88–0.98", SEA_TXT)):
        p = tf2.add_paragraph()
        p.alignment = PP_ALIGN.LEFT
        p.space_after = Pt(9)
        for text, c, sz, b in ((old, MUTED, 19, False), ("  →  ", MUTED, 13, False),
                               (new, col, 19, True), (f"   {label}", MUTED, 12, False)):
            r = p.add_run()
            r.text = text
            r.font.size, r.font.bold, r.font.name = Pt(sz), b, HEAD if sz > 14 else BODY
            r.font.color.rgb = c
    para(tf2, "Reported, not previously reported — the published figures were "
              "pessimistic by 7–11 points for novel organisms, 20–23 with the source present.", size=11.5, color=MUTED, space_after=0, line=1.05)
    tf3 = textbox(s, 7.5, BODY_TOP + 2.75, 5.23, 1.6)
    para(tf3, "Direction of travel matters more than the digits: both corrections were "
              "found by re-deriving the metric from its own definition, and both are now "
              "unit-checked in the scorer.", size=12.5, color=MUTED, first=True,
         space_after=0, line=1.08)
    footer(s, "Manuscript tables, figures and the method paper were regenerated from the corrected scorer")
    return s


def s10_benchmark_design(prs, D):
    s = blank(prs)
    title(s, "Head to head against PICRUSt2",
          sub="The established route from amplicon to function, on exactly the same benchmark")
    cols = [
        ("PICRUSt2 2.6.3", SEA_TXT, [
            "Places each amplicon in a fixed reference phylogeny (EPA-ng, gappa)",
            "Infers gene content by maximum-parsimony ancestral-state reconstruction",
            "27,870 reference genomes; always returns a gene complement",
        ]),
        ("prFBA", TEAL, [
            "Aligns each amplicon to named reference genomes",
            "Unions the gene content those genomes actually carry",
            "~315,000 BV-BRC genomes behind 459,301 unique 16S sequences; abstains below the family floor",
        ]),
        ("Shared namespace", AMBER_TXT, [
            "Truth and both predictions pass through one identical projection",
            "Gene family → product string → EC number → ModelSEED reaction",
            "Scored on the 2,132 ECs both annotation systems can express",
        ]),
    ]
    x, w = M, 3.82
    for head, col, items in cols:
        card = rect(s, x, BODY_TOP, w, 3.3, LIGHT)
        tf = card.text_frame
        tf.margin_left = tf.margin_right = Inches(0.2)
        tf.margin_top = Inches(0.22)
        para(tf, head, size=16, bold=True, color=col, font=HEAD, first=True, space_after=11)
        for it in items:
            hanging(para(tf, "•  " + it, size=12.5, color=INK, space_after=9, line=1.06))
        x += w + 0.335
    tf = textbox(s, M, BODY_TOP + 3.55, W - 2 * M, 1.1)
    bullets(tf, [
        ("Why a shared projection.  ", "Family-level products lose some EC numbers, so "
         "truth and prediction go through the same map and the loss cancels. 953 "
         "amplicons whose source genome yields no EC are unscoreable and excluded."),
    ], size=13.5, space=6)
    footer(s, "Installed without conda: uv venv, user-space R + castor, epa-ng/gappa/HMMER from source (60/61 tests pass)")
    return s


def s11_results(prs, D):
    s = blank(prs)
    title(s, "Coverage against precision",
          sub="Scored against our own annotation — a choice the controlled comparison later removes")
    c = D["cmp"]["by_region"]["ec"]["all"]["all"]
    rows = [
        ("PICRUSt2", c["picrust2"], SEA_TXT),
        ("prFBA, source removed", c["prfba_exclude"], TEAL),
        ("prFBA, source present", c["prfba_include"], MUTED),
    ]
    hdr = ["", "Coverage", "Recall", "Precision", "F1", "Population F1"]
    xs = [M + 0.16, 4.35, 5.95, 7.55, 9.3, 10.75]
    widths = [3.6, 1.5, 1.5, 1.65, 1.5, 1.82]
    y = BODY_TOP + 0.1
    for i, h in enumerate(hdr):
        tf = textbox(s, xs[i], y, widths[i], 0.3)
        para(tf, h, size=11.5, bold=True, color=MUTED,
             align=PP_ALIGN.LEFT if i == 0 else PP_ALIGN.CENTER, first=True, space_after=0)
    y += 0.42
    for name, blk, col in rows:
        rect(s, M, y - 0.1, W - 2 * M, 0.72, LIGHT)
        vals = [f"{blk['coverage']:.3f}",
                f"{blk['conditional']['recall']['value']:.3f}",
                f"{blk['conditional']['precision']['value']:.3f}",
                f"{blk['conditional']['f1']['value']:.3f}",
                f"{blk['population']['f1']['value']:.3f}"]
        tf = textbox(s, xs[0], y, widths[0], 0.5, anchor=MSO_ANCHOR.MIDDLE)
        para(tf, name, size=14, bold=True, color=col, first=True, space_after=0)
        for i, v in enumerate(vals, start=1):
            tf = textbox(s, xs[i], y, widths[i], 0.5, anchor=MSO_ANCHOR.MIDDLE)
            bold = (i == 4 and "removed" in name) or (i == 5 and name == "PICRUSt2")
            para(tf, v, size=14, bold=bold, color=INK if not bold else col,
                 align=PP_ALIGN.CENTER, first=True, space_after=0)
        y += 0.82
    tf = textbox(s, M, y + 0.15, 7.0, 1.8)
    bullets(tf, [
        ("When it answers.  ", "prFBA leads in all eight regions — median per-amplicon "
         "F1 advantage +0.08 to +0.10, Holm p ≈ 0 — but the ranking flips under "
         "PICRUSt2's own annotation, and stays reversed under a neutral key. Read this "
         "slide as a measurement of our annotation, not of the method."),
        ("When it does not.  ", "prFBA abstains on a quarter of novel amplicons, so over "
         "the whole population PICRUSt2's mean F1 is higher. The depth slide shows most "
         "of that abstention is our benchmark's doing, not the method's."),
    ], size=13.5, space=10)
    stat(s, 7.9, y + 0.15, 4.83, "+0.083", "mean EC F1 advantage on the 50,332 amplicons "
         "both answer (0.798 vs 0.715)", h=1.5)
    footer(s, "EC space, n = 67,038 scoreable amplicons; ModelSEED reaction space tracks it within 0.01")
    return s


def s12_novelty(prs, D):
    s = blank(prs)
    title(s, "The failure modes differ",
          sub="Precision as the query moves away from the reference database")
    bins = ["<0.01", "0.01-0.05", "0.05-0.15", "0.15-0.5", ">=0.5"]
    st = D["cmp"]["by_stratum"]["ec"]["nsti_bin"]
    p2 = [st[b]["picrust2"]["conditional"]["precision"]["value"] for b in bins]
    pf = [st[b]["prfba_exclude"]["conditional"]["precision"]["value"] for b in bins]
    cd = CategoryChartData()
    cd.categories = ["< 0.01", "0.01–0.05", "0.05–0.15", "0.15–0.5", "≥ 0.5"]
    cd.add_series("prFBA", pf, number_format="0.00")
    cd.add_series("PICRUSt2", p2, number_format="0.00")
    gf = s.shapes.add_chart(XL_CHART_TYPE.LINE_MARKERS, Inches(M), Inches(BODY_TOP),
                            Inches(7.6), Inches(4.3), cd)
    ch = gf.chart
    style_chart(ch)
    for sr, col in zip(ch.series, (TEAL, SEAFOAM)):
        sr.format.line.color.rgb = col
        sr.format.line.width = Pt(2.5)
        sr.smooth = False
        sr.has_data_labels = True
        dl = sr.data_labels
        dl.show_value = True
        dl.number_format, dl.number_format_is_linked = "0.00", False
        dl.position = XL_LABEL_POSITION.ABOVE
        dl.font.size = Pt(9)
        dl.font.color.rgb = MUTED
        dl.font.name = BODY
    ch.value_axis.maximum_scale = 1.0
    ch.value_axis.minimum_scale = 0.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    tf = textbox(s, 8.5, BODY_TOP + 0.1, 4.23, 4.3)
    bullets(tf, [
        ("PICRUSt2.  ", "Precision falls 0.73 → 0.23 while recall barely moves "
         "(0.87 → 0.77). Distant placements keep asserting gene content, and most of "
         "the extra assertions are wrong."),
        ("prFBA.  ", "Precision falls only 0.83 → 0.61 — not because it abstains (it "
         "still answers 81.5% at NSTI ≥ 0.5, against 76.4% in the nearest bin) but because "
         "it asserts far smaller gene sets: ~135 ECs per answer there, against ~406."),
        ("Why it matters.  ", "A model built from an abstention is one you know not to "
         "build. A model built from a confidently wrong genome is a plausible artifact "
         "that propagates into every flux prediction."),
    ], size=13, space=11)
    footer(s, "NSTI = nearest-sequenced-taxon index, PICRUSt2's own distance measure")
    return s


def s13_ceiling(prs, D):
    """What the control measures: two annotations of one genome, no prediction."""
    s = blank(prs)
    background(s, DEEP)
    title(s, "A control: how much do the two annotations even agree?", dark=True,
          sub="Before reading any score, ask what the scoring truth is made of")
    sh = D["ceil"]["ceiling"]["shared"]

    tf = textbox(s, M, BODY_TOP, 6.5, 3.6)
    bullets(tf, [
        ("The setup.  ", "2,689 scoreable benchmark genomes sit in both databases. BV-BRC annotated "
         "that genome and so did the PICRUSt2 project — independently, with different "
         "gene callers and EC rules."),
        ("The thought experiment.  ", "Suppose PICRUSt2 places an amplicon perfectly, on "
         "its own source genome's tip. What it returns is that genome's stored row. So "
         "compare the two annotations of the one genome directly."),
        ("No prediction is involved.  ", "Any disagreement here is two projects reading "
         "the same DNA differently. It sets the reference level for what PICRUSt2 can "
         "score against our truth — not a strict bound: for 43% of in-reference amplicons "
         "PICRUSt2 is more precise than it."),
    ], size=14, space=13, color=WHITE, head_color=MINT)

    gx, gy, gw = 7.5, BODY_TOP - 0.05, 5.23
    chip = rect(s, gx + 1.35, gy, 2.5, 0.52, RGBColor(0x0A, 0x4A, 0x55))
    ctf = chip.text_frame
    ctf.vertical_anchor = MSO_ANCHOR.MIDDLE
    ctf.margin_left = ctf.margin_right = 0
    para(ctf, "one genome, e.g. GCF_002009845", size=11.5, bold=True, color=WHITE,
         align=PP_ALIGN.CENTER, first=True, space_after=0)
    for i, (who, what) in enumerate((("BV-BRC says", "our truth: gene set A"),
                                     ("PICRUSt2 says", "its reference row: gene set B"))):
        card = rect(s, gx + i * 2.73, gy + 0.72, 2.5, 0.95, RGBColor(0x0A, 0x4A, 0x55))
        ctf = card.text_frame
        ctf.margin_left = ctf.margin_right = Inches(0.14)
        ctf.margin_top = Inches(0.12)
        para(ctf, who, size=12.5, bold=True, color=MINT, font=HEAD, first=True, space_after=4)
        para(ctf, what, size=11.5, color=WHITE, space_after=0, line=1.0)
    tf2 = textbox(s, gx, gy + 1.82, gw, 0.3)
    para(tf2, "A compared with B — for the same organism:", size=12, bold=True,
         color=WHITE, first=True, space_after=0)
    for i, (v, lab) in enumerate((
            (f"{sh['recall']['mean']*100:.1f}%", "of A's EC numbers also appear in B"),
            (f"{sh['precision']['mean']*100:.1f}%", "of B's EC numbers also appear in A"),
            (f"{sh['jaccard']['mean']:.3f}", "Jaccard overlap of the two sets"))):
        card = rect(s, gx, gy + 2.2 + i * 0.78, gw, 0.66, RGBColor(0x0A, 0x4A, 0x55))
        ctf = card.text_frame
        ctf.margin_left = Inches(0.2)
        ctf.margin_right = Inches(0.2)
        ctf.vertical_anchor = MSO_ANCHOR.MIDDLE
        pp = ctf.paragraphs[0]
        pp.alignment = PP_ALIGN.LEFT
        r = pp.add_run()
        r.text = v
        r.font.size, r.font.bold, r.font.name = Pt(22), True, HEAD
        r.font.color.rgb = MINT
        r2 = pp.add_run()
        r2.text = "   " + lab
        r2.font.size, r2.font.name = Pt(11), BODY
        r2.font.color.rgb = WHITE
    footer(s, "A quarter of PICRUSt2's enzyme calls, and a tenth of BV-BRC's, have no "
              "counterpart in the other — neither is wrong, they are different annotation "
              "pipelines", dark=True)
    s.notes_slide.notes_text_frame.text = (
        "This is a control, not a result. It says the measuring stick itself is only "
        "accurate to about 0.69 Jaccard when comparing across the two annotation systems.")
    return s


def s13b_dialect(prs, D):
    """What the control means: each method scores best in its own dialect."""
    s = blank(prs)
    title(s, "Each method wins at home",
          sub="Identical predictions, scored against two different answer keys for the "
              "same genome")
    sc = D["dialect"]["by_scope"]["conditional"]
    n = sc["picrust2"]["bvbrc"]["f1"]["n"]

    hdr_y = BODY_TOP + 0.34
    tfk = textbox(s, M, BODY_TOP - 0.06, W - 2 * M, 0.34)
    para(tfk, "Same 16,282 amplicons, same predictions, same EC vocabulary — only the "
              "answer key changes between the two columns", size=12.5, color=MUTED,
         italic=True, first=True, space_after=0)
    cols = [("", M + 0.16, 4.3), ("Answer key:\nBV-BRC's annotation", 5.1, 3.4),
            ("Answer key:\nPICRUSt2's annotation", 8.8, 3.6)]
    for text, x, w in cols:
        tf = textbox(s, x, hdr_y, w, 0.62)
        for i, line in enumerate(text.split("\n")):
            para(tf, line, size=12, bold=True, color=MUTED,
                 align=PP_ALIGN.LEFT if x == M + 0.16 else PP_ALIGN.CENTER,
                 first=(i == 0), space_after=0)
    y = hdr_y + 0.72
    for name, key, col in (("PICRUSt2", "picrust2", SEA_TXT),
                           ("prFBA (source removed)", "prfba_exclude", TEAL),
                           ("prFBA (source present)", "prfba_include", MUTED)):
        rect(s, M, y - 0.1, W - 2 * M, 0.74, LIGHT)
        tf = textbox(s, M + 0.16, y, 4.3, 0.54, anchor=MSO_ANCHOR.MIDDLE)
        para(tf, name, size=14, bold=True, color=col, first=True, space_after=0)
        for x, w, t in ((5.1, 3.4, "bvbrc"), (8.8, 3.6, "picrust2")):
            v = sc[key][t]["f1"]["mean"]
            best = (t == "bvbrc" and key == "prfba_exclude") or \
                   (t == "picrust2" and key == "picrust2")
            tf = textbox(s, x, y, w, 0.54, anchor=MSO_ANCHOR.MIDDLE)
            para(tf, f"{v:.3f}", size=17 if best else 15, bold=best,
                 color=col if best else INK, align=PP_ALIGN.CENTER, first=True, space_after=0)
        y += 0.80
    tf = textbox(s, M, y + 0.10, 7.3, 1.9)
    bullets(tf, [
        ("The flip.  ", "Under our answer key prFBA leads by 0.050 F1; under PICRUSt2's "
         "it trails by 0.204. The diagonal is near-tautological — each system is handed "
         "back its own annotation of a genome it already contains — and each off-diagonal "
         "cell is roughly 0.82–0.88 of the home score."),
        ("So do not claim.  ", "“prFBA predicts function more accurately than "
         "PICRUSt2.” This benchmark cannot support it in either direction — so we built "
         "one that can. It is the next three slides, and it does not favour us."),
    ], size=13.5, space=11)
    card = rect(s, 8.1, y + 0.10, 4.63, 1.8, LIGHT)
    ctf = card.text_frame
    ctf.margin_left = ctf.margin_right = Inches(0.2)
    ctf.margin_top = Inches(0.16)
    para(ctf, "What does survive", size=13.5, bold=True, color=INK, font=HEAD, first=True,
         space_after=8)
    for t in ("Coverage and abstention — a design difference, not a measurement",
              "The shape of precision decay with novelty, within each method",
              "prFBA's own validation — its output is a BV-BRC genome, so that "
              "truth is the operative one"):
        hanging(para(ctf, "•  " + t, size=11.5, color=MUTED, space_after=5, line=1.0))
    footer(s, f"Conditional on prFBA answering, n = {n:,} amplicons whose source genome is in both databases")
    s.notes_slide.notes_text_frame.text = (
        "A fair cross-method ranking needs a truth neither system authors. That is "
        "exactly what the shared-reference k-fold benchmark does, and it is why its "
        "result is credible against our own interest.")
    return s


def s14_depth(prs, D):
    s = blank(prs)
    title(s, "Most of our reported abstentions were our own doing",
          sub="Self-exclusion deletes hits from a record already truncated to 20 candidates")
    pr = D["probe"]
    dn = pr["depth_needed"]
    order = ["21-50", "51-100", "101-200", "201-500"]
    total = sum(dn.values())
    cum, ys = 0, []
    for k in order:
        cum += dn.get(k, 0)
        ys.append(cum / total)
    cd = CategoryChartData()
    cd.categories = ["50", "100", "200", "500"]
    cd.add_series("Recovered", ys, number_format="0%")
    gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(BODY_TOP + 1.92),
                            Inches(5.7), Inches(2.75), cd)
    ch = gf.chart
    style_chart(ch, legend=False)
    series_color(ch.series[0], TEAL)
    ch.value_axis.maximum_scale = 1.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    sr = ch.series[0]
    sr.has_data_labels = True
    dl = sr.data_labels
    dl.show_value = True
    dl.number_format, dl.number_format_is_linked = "0%", False
    dl.position = XL_LABEL_POSITION.OUTSIDE_END
    dl.font.size = Pt(9)
    dl.font.color.rgb = MUTED
    dl.font.name = BODY
    tf0 = textbox(s, M, BODY_TOP + 1.62, 5.7, 0.3)
    para(tf0, "Amplicons recovered, by how many hits we keep", size=11.5, bold=True,
         color=MUTED, first=True, space_after=0)
    stat(s, M, BODY_TOP, 2.7, "94.6%", "of self-exclusion abstentions\nwere emptied records",
         color=AMBER_TXT)
    stat(s, M + 3.0, BODY_TOP, 2.7, "1.000", "include-self coverage\nin every identity bin")
    tf = textbox(s, 6.5, BODY_TOP, 6.23, 4.3)
    bullets(tf, [
        ("What happens.  ", "For an organism with many near-identical relatives, all 20 "
         "retained hits exceed the exclusion threshold. The filter deletes every one and "
         "the selector abstains for want of input — not for want of a qualifying "
         "reference, which the database still holds."),
        ("Not a production problem.  ", "Nothing deletes hits in real use, and a genuinely "
         "novel organism has no near-identical relatives to remove. Include-self coverage "
         "is always 1.000."),
        ("Measured, not assumed.  ", "320 affected amplicons re-aligned keeping 500 hits "
         "instead of 20: 88.4% recover an answer, at EC F1 0.783 against 0.798 for "
         "amplicons that already answered — ordinary calls, not dregs."),
        ("So.  ", "Every self-exclusion coverage and population number we report is a lower "
         "bound; the conditional ones are not. Coverage would "
         "rise 0.77 → ~0.96 and population F1 0.60 → ~0.75 — under our own annotation as "
         "the answer key. The controlled comparison puts the neutral-key ceiling for a "
         "fully-answering prFBA at −0.03 against PICRUSt2, not ahead of it."),
    ], size=13, space=10)
    footer(s, "The prefilter already scores 500 candidates, so keeping 100 costs essentially nothing")
    return s


def s15_next(prs, D):
    s = blank(prs)
    background(s, DEEP)
    title(s, "Where this leaves us", dark=True)
    items = [
        ("Report the controlled result as the headline", "PICRUSt2 is more accurate at EC "
         "level at genus and family rank, and prFBA forced to answer only ties copying the "
         "nearest reference. The manuscript now says so.", MINT),
        ("Lead with what replicates", "Conditional answer quality is robust to panel "
         "composition. The ATP phenotype shows no detectable mean difference at n = 755 — "
         "worth reporting, and worth not overselling.", MINT),
        ("Test the weighting, which nothing yet does", "Calibrate p(function) against the "
         "truth genomes the panels already provide, and run prFBA's objective on the "
         "probability-weighted genomes it was built for.", MINT),
        ("Raise retention depth to 100", "The self-exclusion coverage figures are a "
         "property of the 20-hit record; the prefilter already scores 500 candidates, so "
         "the fix costs almost nothing.", MINT),
        ("Close the loop to flux", "Propagate probability-weighted genomes through "
         "community prFBA and validate against paired metagenomes.", MINT),
    ]
    y = BODY_TOP + 0.05
    for i, (head, body, col) in enumerate(items):
        dot = rect(s, M, y + 0.04, 0.42, 0.42, RGBColor(0x0A, 0x4A, 0x55), shape=MSO_SHAPE.OVAL)
        dtf = dot.text_frame
        dtf.margin_left = dtf.margin_right = dtf.margin_top = dtf.margin_bottom = 0
        dtf.vertical_anchor = MSO_ANCHOR.MIDDLE
        para(dtf, str(i + 1), size=14, bold=True, color=MINT, font=HEAD,
             align=PP_ALIGN.CENTER, first=True, space_after=0)
        tf = textbox(s, M + 0.62, y, 11.4, 1.0)
        para(tf, head, size=16, bold=True, color=WHITE, font=HEAD, first=True, space_after=4)
        para(tf, body, size=12, color=MINT, space_after=0, line=1.04)
        y += 1.06
    footer(s, "github.com/freiburgermsu/prFBA  ·  github.com/freiburgermsu/prFBA_paper", dark=True)
    s.notes_slide.notes_text_frame.text = (
        "One follow-up is nearly free and changes headline numbers: the depth re-run. "
        "Calibrating the weighting is the other gating experiment, because nothing "
        "reported yet tests it.")
    return s


# ------------------------------------------------ the MAG arm (synthetic MAGs)
def s16_mag_arm(prs, D):
    s = blank(prs)
    title(s, "A second arm: whole genomes instead of amplicons",
          sub="The same 10,000 organisms, degraded into synthetic MAGs and matched by whole-genome ANI")
    m = D["mag"]["summary"]
    tiers = ["c1_x0_exclude", "c0.9_x0_exclude", "c0.7_x0_exclude",
             "c0.5_x0_exclude", "c0.3_x0_exclude"]
    cov = [m[t]["coverage_answered"] for t in tiers]
    f1 = [m[t]["genecap_cond_f1"]["mean"] for t in tiers]
    gen = [m[t]["Genus"]["correct"] for t in tiers]
    # the 0%-contamination grid only (the contaminated half was withdrawn); each MAG is
    # scored twice, once per self-inclusion regime
    n_cells_clean = sum(v["n"] for k, v in m.items() if "_x0_" in k)
    n_mags = n_cells_clean // 2
    for i, (v, lab, col) in enumerate([
            (f"{n_mags:,}", "synthetic MAGs from\n9,997 panel genomes", TEAL),
            (f"{n_cells_clean:,}", "scored cells across five\ncompleteness tiers", TEAL),
            (f"{f1[0]:.3f}", f"conditional gene-capture F1\n— {f1[-1]:.3f} at 30% complete", SEA_TXT),
            (f"{gen[0]:.3f}", f"genus accuracy\n— {gen[-1]:.3f} at 30% complete", SEA_TXT)]):
        stat(s, M + i * 3.108, BODY_TOP, 2.808, v, lab, color=col)

    tf = textbox(s, M, BODY_TOP + 1.62, 6.5, 3.5)
    bullets(tf, [
        ("The design.  ", "Each panel genome is cut into ~50 kb synthetic contigs — "
         "removing a contiguity confound across genomes spanning 1 to 5,265 deposited "
         "contigs — then fragments are dropped at random to 100 / 90 / 70 / 50 / 30% "
         "completeness. The draft is matched by skani ANI instead of 16S identity."),
        ("Same everything else.  ", "The amplicon arm's truth, metric definitions, "
         "three-state cells and conditional-vs-population denominators — not "
         "reimplementations — so the two arms are directly comparable."),
        ("Completeness is a null axis.  ", "Nothing the pipeline decides moves with it: "
         "not which reference is selected, not whether one is selected, not how much "
         "metabolism is recovered."),
        ("Why.  ", "With the source present, AF_query stays flat at ~99.8% — a clean "
         "contig subset still aligns end to end — while AF_ref ≈ 0.99 × completeness; "
         "under exclusion AF_query stays flat too and AF_ref still scales. The gate reads "
         "the first and ignores the second, a free completeness signal the pipeline "
         "already computes and discards."),
    ], size=13.5, space=10)

    cd = CategoryChartData()
    cd.categories = ["100%", "90%", "70%", "50%", "30%"]
    cd.add_series("Conditional gene-capture F1", f1, number_format="0.000")
    cd.add_series("Resolve rate (coverage)", cov, number_format="0.000")
    gf = s.shapes.add_chart(XL_CHART_TYPE.LINE_MARKERS, Inches(7.4), Inches(BODY_TOP + 1.72),
                            Inches(5.33), Inches(3.2), cd)
    ch = gf.chart
    style_chart(ch)
    for sr, col in zip(ch.series, (TEAL, AMBER)):
        sr.format.line.color.rgb = col
        sr.format.line.width = Pt(2.5)
        sr.smooth = False
        sr.has_data_labels = True
        dl = sr.data_labels
        dl.show_value = True
        dl.number_format, dl.number_format_is_linked = "0.000", False
        dl.position = XL_LABEL_POSITION.ABOVE
        dl.font.size = Pt(9)
        dl.font.color.rgb = MUTED
        dl.font.name = BODY
    ch.value_axis.maximum_scale = 1.0
    ch.value_axis.minimum_scale = 0.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    tfc = textbox(s, 7.4, BODY_TOP + 1.42, 5.33, 0.3)
    para(tfc, "Self-exclusion, by MAG completeness", size=11.5, bold=True, color=MUTED,
         first=True, space_after=0)
    footer(s, "0%-contamination grid only: the synthesizer that built this grid overshot "
              "its dose at low completeness, so every contamination claim was withdrawn")
    s.notes_slide.notes_text_frame.text = (
        "Two reasons this arm exists: most real MAGs carry no usable 16S at all, and the "
        "proposal asked at what completeness a MAG stops being worth more than an "
        "amplicon. The answer to the second is that nothing moves between 100% and 30%.")
    return s


def s17_crossover(prs, D):
    s = blank(prs)
    title(s, "There is no completeness at which the amplicon catches up",
          sub="Paired within organism, down to a 1%-complete draft — a 37 kb fragment")
    rows = D["magx"]
    cats = [f"{int(r['completeness'] * 100)}%" for r in rows]
    mag = [r["FullLength16S_mag"] for r in rows]
    amp = [r["FullLength16S_amp"] for r in rows]
    top, bot = rows[0], rows[-1]

    cd = CategoryChartData()
    cd.categories = cats
    cd.add_series("MAG, whole-genome ANI", mag, number_format="0.000")
    cd.add_series("Full-length 16S amplicon", amp, number_format="0.000")
    gf = s.shapes.add_chart(XL_CHART_TYPE.LINE_MARKERS, Inches(M), Inches(BODY_TOP + 0.3),
                            Inches(7.3), Inches(4.0), cd)
    ch = gf.chart
    style_chart(ch)
    for sr, col in zip(ch.series, (TEAL, AMBER)):
        sr.format.line.color.rgb = col
        sr.format.line.width = Pt(2.5)
        sr.smooth = False
    ch.value_axis.maximum_scale = 0.88
    ch.value_axis.minimum_scale = 0.60
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    tfc = textbox(s, M, BODY_TOP, 7.3, 0.3)
    para(tfc, "Conditional gene-capture F1 against MAG completeness  (note the axis "
              "starts at 0.60)", size=11.5, bold=True, color=MUTED, first=True, space_after=0)

    stat(s, 8.1, BODY_TOP, 4.63,
         f"+{top['FullLength16S_delta']:.3f}  →  +{bot['FullLength16S_delta']:.3f}",
         "MAG advantage in gene-capture F1 over a full-length amplicon, from a complete "
         "genome down to a 1% draft", vsize=24, h=1.4)
    tf = textbox(s, 8.1, BODY_TOP + 1.6, 4.63, 3.4)
    bullets(tf, [
        ("Why it holds.  ", "The prediction is the gene content of the selected "
         "reference, not of the MAG. Completeness matters only if it changes which "
         "reference is chosen, and ANI is robust enough that it rarely does."),
        ("Against V4.  ", f"The margin is larger still — {sg(top['V4_delta'])} at 100% "
         f"and {sg(bot['V4_delta'])} at 1%."),
        ("What does degrade.  ", f"Coverage, not quality: the resolve rate falls "
         f"{rd(top['resolve_rate'])} → {rd(bot['resolve_rate'])}, flat to ~15% and then "
         "declining as the remaining sequence stops clearing skani's screen."),
        ("Do not over-read it.  ", "This is conditional on both arms answering. Over the "
         "whole population the amplicon arm wins, 0.535 to 0.248, because the MAG arm "
         "abstains on 69% of organisms."),
    ], size=13, space=10)
    footer(s, f"Means of paired per-organism F1, n = 1,260–1,580 organisms per tier; "
              f"median paired delta +{top['FullLength16S_median_delta']:.3f} → "
              f"+{bot['FullLength16S_median_delta']:.3f}; Wilcoxon p < 1e-80 throughout")
    s.notes_slide.notes_text_frame.text = (
        "At 1% completeness — 50x below MIMAG's medium-quality floor, not a MAG any "
        "study would keep — whole-genome ANI still beats a full-length amplicon.")
    return s


def s18_abstain(prs, D):
    s = blank(prs)
    res = D["mag"]["summary"]["c1_x0_exclude"]["coverage_answered"]
    title(s, f"The MAG arm abstains on {100 * (1 - res):.0f}% of organisms — and should", sub="One gate, not the ANI floor, accounts for nearly all of it")
    labels = [t for t, _ in AF_SWEEP][::-1]
    vals = [v for _, v in AF_SWEEP][::-1]
    cd = CategoryChartData()
    cd.categories = labels
    cd.add_series("Share with a qualifying hit", vals, number_format="0.000")
    gf = s.shapes.add_chart(XL_CHART_TYPE.BAR_CLUSTERED, Inches(M), Inches(BODY_TOP + 0.32),
                            Inches(6.3), Inches(3.3), cd)
    ch = gf.chart
    style_chart(ch, legend=False, size=10.5)
    series_color(ch.series[0], TEAL)
    ch.value_axis.maximum_scale = 0.72
    ch.value_axis.minimum_scale = 0.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    sr = ch.series[0]
    sr.has_data_labels = True
    dl = sr.data_labels
    dl.show_value = True
    dl.number_format, dl.number_format_is_linked = "0.000", False
    dl.position = XL_LABEL_POSITION.OUTSIDE_END
    dl.font.size = Pt(9)
    dl.font.color.rgb = MUTED
    dl.font.name = BODY
    tfc = textbox(s, M, BODY_TOP, 6.3, 0.3)
    para(tfc, f"Qualifying-hit share as the gates are relaxed (resolve rate {rd(res)})",
         size=11.5, bold=True,
         color=MUTED, first=True, space_after=0)

    tf = textbox(s, 7.1, BODY_TOP, 5.63, 3.6)
    bullets(tf, [
        ("Marker reach is the cause.  ", "At genus distance two genomes keep high local "
         "identity over only a small part of their length. Among the 6,011 organisms "
         "that do have a same-genus congener in the panel, the best surviving hit "
         "aligns over a median 20.3% of the query — and 45.1% fall under the 15% gate."),
        ("An amplicon has no such problem.  ", "It aligns end to end by construction. In "
         "the 0.945–0.987 identity band the amplicon answers 98–99% of the time; "
         "whole-genome ANI answers 8–20%."),
        ("Reference density, not the marker.  ", f"Against GTDB r226 — {GTDB_FOLD}× denser "
         f"— the same pipeline resolves {rd(D['gtdb']['summary']['c1_x0_exclude']['resolve_rate'])} "
         f"of {D['gtdb']['n_genomes']:,} organisms instead of {rd(res)}, with "
         "completeness-invariance intact."),
    ], size=13, space=11)

    card = rect(s, M, 5.08, 6.3, 1.72, LIGHT)
    ctf = card.text_frame
    ctf.margin_left = ctf.margin_right = Inches(0.2)
    ctf.margin_top = Inches(0.16)
    para(ctf, "The gate is not a defect — it is the method working", size=13.5, bold=True,
         color=INK, font=HEAD, first=True, space_after=8)
    para(ctf, "The 1,367 organisms a 2% gate would newly resolve come back at mean "
              "gene-capture F1 0.599 against 0.807 for those already resolving, and only "
              "27.0% of them exceed the full-length 16S amplicon's mean of 0.709. prFBA abstains "
              "where answering would be worse than the marker it replaces.",
         size=12, color=MUTED, space_after=0, line=1.06)

    card2 = rect(s, 7.1, 5.08, 5.63, 1.72, LIGHT)
    c2 = card2.text_frame
    c2.margin_left = c2.margin_right = Inches(0.2)
    c2.margin_top = Inches(0.16)
    para(c2, "Neither is it the amplicon arm's truncation artifact", size=13.5, bold=True,
         color=INK, font=HEAD, first=True, space_after=8)
    para(c2, "Re-searching 647 probe queries at depth 2,000 leaves every uncapped "
             "selection identical (200/200), and depth 200 is already fully converged. "
             "Truncation censors 0.16% of MAG queries, against 94.6% of the amplicon "
             "arm's self-exclusion abstentions.",
         size=12, color=MUTED, space_after=0, line=1.06)
    footer(s, "Self-exclusion at ANI 95, 100% completeness; all 11 selections that depth "
              "did change were abstention → resolution, never one reference for another")
    return s


def s19_two_panel(prs, D):
    s = blank(prs)
    title(s, "Does the benchmark panel decide the answer?",
          sub="Three 10,000-organism panels built on different principles — taxonomy-stratified, embedding-uniform, and a random null")
    regs = ["FullLength16S", "V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V4-V5_944R", "V7-V9"]
    cd = CategoryChartData()
    cd.categories = [REG_LABEL[r] for r in regs]
    for key, lab in (("original", "Taxonomy-stratified"), ("kmeans", "k-means"),
                     ("random", "Uniform random")):
        cd.add_series(lab, [tp3_cell(D, r, "exclude_0.987", "genecap_cond_f1")[key]["crude"]
                            for r in regs], number_format="0.00")
    gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(BODY_TOP + 0.3),
                            Inches(7.5), Inches(3.05), cd)
    ch = gf.chart
    style_chart(ch, size=10.5)
    for sr, col in zip(ch.series, (TEAL, AMBER, SEAFOAM)):
        series_color(sr, col)
    ch.value_axis.maximum_scale = 0.8
    ch.value_axis.minimum_scale = 0.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    tfc = textbox(s, M, BODY_TOP, 7.5, 0.3)
    para(tfc, "Conditional gene-capture F1 by region, self-exclusion — the three panels nearly agree",
         size=11.5, bold=True, color=MUTED, first=True, space_after=0)

    card = rect(s, M, 4.98, 7.5, 1.82, LIGHT)
    ctf = card.text_frame
    ctf.margin_left = ctf.margin_right = Inches(0.2)
    ctf.margin_top = Inches(0.16)
    para(ctf, "Coverage is not panel-invariant — because retention depth is not",
         size=13.5, bold=True, color=INK, font=HEAD, first=True, space_after=8)
    gcov = tp3_cell(D, "FullLength16S", "exclude_0.987", "genecap_coverage")
    para(ctf, f"On full-length 16S, self-exclusion empties the whole 20-hit record for "
              f"{EMPTIED['random']}% of random-panel queries, {EMPTIED['kmeans']}% of "
              f"k-means and {EMPTIED['original']}% of the original's. Where the record "
              "survives, only 2.5–3.3% go uncovered on any panel, so emptied records are "
              f"nearly all of the coverage gap ({gcov['original']['crude']:.2f} against "
              f"{gcov['random']['crude']:.2f}).",
         size=12, color=MUTED, space_after=0, line=1.06)

    tf = textbox(s, 8.3, BODY_TOP, 4.43, 5.3)
    bullets(tf, [
        ("Three compositions.  ", f"Stratified to over-weight rare taxa; k-means++ over "
         f"{D['knn']['n']:,} 16S embeddings; and a uniform random draw over the same "
         "admissible organisms — the null, which inherits the database's own skew."),
        ("What replicates.  ", "Conditional gene capture, everywhere: within 0.026 between "
         "any two panels once phylum mix is standardized, in both arms."),
        ("What does not.  ", "Self-exclusion coverage, which falls 0.27–0.31 on the random "
         "panel in every region and survives standardization — a property of record "
         "depth, not of the scheme's taxa."),
        ("The defensible claim.  ", "Gene capture given a reference does not depend on how "
         "the panel was built. Coverage at fixed depth measures panel redundancy."),
    ], size=12.5, space=10)
    footer(s, "Crude values; the phylum-standardized comparison follows. All three panels "
              "scored against one PGFam cache state, verified intact against a snapshot")
    return s


def s20_panel_space(prs, D):
    s = blank(prs)
    title(s, "Where each panel samples the 16S space",
          sub="Local density measured in the full 1024-d embedding, not in the projection")
    # three_panel/plot_space_three.py: A k-means, B original, C random, D decile shares
    img = os.path.join(REPO, "three_panel", "F_panel_space_three.png")
    s.shapes.add_picture(img, Inches(M), Inches(1.45), width=Inches(W - 2 * M))
    d = D["dec3"]
    sp = lambda k: 100 * sum(d[k][:2])        # noqa: E731  two sparsest deciles
    tf = textbox(s, M, 5.32, W - 2 * M, 1.2)
    para(tf, f"The stratified panel over-weights the sparse space hardest — {sp('og'):.1f}% of "
             f"its picks from the two sparsest deciles, against {sp('km'):.1f}% for k-means — "
             f"while the random panel (C) tracks the population: {sp('rnd'):.1f}% "
             f"against {sp('population'):.1f}%, and {100 * d['rnd'][9]:.1f}% of picks in the "
             f"densest decile against the population's {100 * d['population'][9]:.1f}% "
             f"(stratified {100 * d['og'][9]:.1f}%, k-means {100 * d['km'][9]:.1f}%). The dense "
             "deciles are where self-exclusion empties a 20-hit record.",
         size=12, color=MUTED, first=True, space_after=6, line=1.06)
    para(tf, "Black ticks in D are each decile's share of all sequences: the density measure "
             "is heavily tied (166 distinct values), so the deciles hold 7.7–12.6% of "
             "sequences rather than 10% each.", size=11, color=MUTED, space_after=0, line=1.04)
    footer(s, f"{d['n_og']:,} of the original panel's 10,000 organisms have a position in "
              f"this space (the rest have short or atypical 16S, or a sequence the index credits to "
              f"no admissible taxon); all 10,000 k-means and random picks do. Shares are of "
              f"each panel's on-frame picks")
    s.notes_slide.notes_text_frame.text = (
        "Density ranking is stable in k — Spearman "
        f"{D['knn']['spearman_21_50']:.3f} between the 21st and 50th nearest neighbour. "
        "Panel C is the uniform-random panel; all 10,000 of its organisms sit on the "
        "clusterable set. The manuscript's figure (F16) is still the two-panel version.")
    return s


# ----------------------------------------- the controlled, shared-reference arm
def s15_kfold_design(prs, D):
    s = blank(prs)
    title(s, "A controlled comparison: one reference set, one answer key", tag="NEW",
          sub="Neither system authors the truth it is scored against")
    g, f = D["kg"], D["kf"]
    mg = g["meta"]
    for i, (v, lab, col) in enumerate([
            (f"{mg['n_eligible']:,}", "eligible panel genomes,\nspecies-dereplicated", TEAL),
            (f"{g['headline']['reference_density_genomes']:,}", "reference genomes per fold\n(median, leave-genus-out)", TEAL),
            (f"{mg['n_cells_scored']:,}", "scored holdout cells\nat genus; "
             f"{f['meta']['n_cells_scored']:,} at family", SEA_TXT),
            ("0.02", "smallest effect worth caring\nabout, fixed in advance", AMBER_TXT)]):
        stat(s, M + i * 3.108, BODY_TOP, 2.808, v, lab, color=col)

    tf = textbox(s, M, BODY_TOP + 1.62, 7.2, 3.5)
    bullets(tf, [
        ("The fix.  ", "Both methods predict the EC content of held-out genomes from the "
         "same reference set, and one annotation — BV-BRC's EC projection — is both "
         "methods' gene content and the answer key. The dialect problem disappears "
         "because there is only one dialect."),
        ("Held out by clade, not at random.  ", "Leave-genus-out and leave-family-out, "
         "k = 5 folds (k = 10 was impossible: one leakage component holds 23.5% of the panel), "
         "so a holdout genome's relatives leave the reference with it. "
         "Significance is a cluster bootstrap over holdout groups, not over amplicons."),
        ("Chosen in advance.  ", "The effect size worth caring about, the selection mode "
         "(leave-one-fold-out selection, which picked the same mode in all five folds), the "
         "controls and the exit criteria were fixed before the run."),
    ], size=13.5, space=11)

    card = rect(s, 8.0, BODY_TOP + 1.62, 4.73, 3.5, LIGHT)
    ctf = card.text_frame
    ctf.margin_left = ctf.margin_right = Inches(0.2)
    ctf.margin_top = Inches(0.18)
    para(ctf, "Five baselines, because a benchmark without a floor is a press release",
         size=13, bold=True, color=INK, font=HEAD, first=True, space_after=9)
    for t in ("constant — the most common ECs, every time",
              "nearest neighbour — copy the single closest reference",
              "family donor — copy a random same-family genome",
              "clade-max — best nearby genome, chosen against the truth: an oracle",
              "random abstain — prFBA's own forced answers, declined at random to its coverage"):
        hanging(para(ctf, "•  " + t, size=12, color=MUTED, space_after=7, line=1.04))
    footer(s, f"Panel caveat: an HMM check finds {D['hmm']['n_not_16s']} of "
              f"{D['hmm']['n_records']:,} panel records are not 16S at all and "
              f"{D['hmm']['n_low_cover']} cover less than 80% of the model")
    s.notes_slide.notes_text_frame.text = (
        "This is the experiment the dialect slide said was needed. It was run, and the "
        "result does not favour us — which is why it is worth showing.")
    return s


def s16_kfold_result(prs, D):
    s = blank(prs)
    title(s, "Under a neutral answer key, the ranking reverses", tag="NEW",
          sub="Population EC-F1, leave-genus-out — and every method sits near a trivial baseline")
    g, f = D["kg"], D["kf"]
    arms = [("base_constant", "Constant (floor)"),
            ("base_random_abstain", "Random, prFBA's coverage"),
            ("prfba_nestedcv", "prFBA"),
            ("base_nn1", "Nearest neighbour"),
            ("prfba_nestedcv_forced", "prFBA, forced to answer"),
            ("picrust2", "PICRUSt2"),
            ("base_clade_max", "Clade-max oracle (ceiling)")]
    vals = [g["arms"][k]["ec"]["population"]["f1"]["value"] for k, _ in arms]
    cd = CategoryChartData()
    cd.categories = [lab for _, lab in arms][::-1]
    cd.add_series("Population EC-F1", vals[::-1], number_format="0.000")
    gf = s.shapes.add_chart(XL_CHART_TYPE.BAR_CLUSTERED, Inches(M), Inches(BODY_TOP + 0.3),
                            Inches(6.6), Inches(4.0), cd)
    ch = gf.chart
    style_chart(ch, legend=False, size=10.5)
    series_color(ch.series[0], TEAL)
    ch.value_axis.maximum_scale = 0.85
    ch.value_axis.minimum_scale = 0.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    sr = ch.series[0]
    sr.has_data_labels = True
    dl = sr.data_labels
    dl.show_value = True
    dl.number_format, dl.number_format_is_linked = "0.000", False
    dl.position = XL_LABEL_POSITION.OUTSIDE_END
    dl.font.size = Pt(9)
    dl.font.color.rgb = MUTED
    dl.font.name = BODY
    tfc = textbox(s, M, BODY_TOP, 6.6, 0.3)
    para(tfc, "Every arm scored on the same cells, same answer key", size=11.5,
         bold=True, color=MUTED, first=True, space_after=0)

    cg = g["contrasts"]["prfba_nestedcv__vs__picrust2"]["ec"]
    cf = f["contrasts"]["prfba_nestedcv__vs__picrust2"]["ec"]
    fg = g["contrasts"]["prfba_nestedcv_forced__vs__picrust2"]["ec"]
    stat(s, 7.5, BODY_TOP, 5.23, f"{sg(cg['delta_f1_population'], 3)}",
         f"prFBA minus PICRUSt2 at genus, 95% CI [{cg['ci95'][0]:.3f}, {cg['ci95'][1]:.3f}] "
         f"— {sg(cf['delta_f1_population'], 3)} at family", color=AMBER_TXT, h=1.4)
    tf = textbox(s, 7.5, BODY_TOP + 1.6, 5.23, 3.4)
    cnn = g["contrasts"]["prfba_nestedcv__vs__base_nn1"]["ec"]
    p2nn = g["contrasts"]["picrust2__vs__base_nn1"]["ec"]
    bullets(tf, [
        ("Two thirds of the gap is abstention.  ", f"Forced to answer, prFBA trails PICRUSt2 by "
         f"{abs(fg['delta_f1_population']):.3f} — closer, not level. Refusing costs most of the "
         "deficit; the rest is in the answers themselves."),
        ("The bar the design set.  ", "Copying the single nearest surviving reference. "
         "prFBA forced only ties it, a margin with no paired test; the shipped, abstaining "
         f"arm falls {abs(cnn['delta_f1_population']):.3f} below it, while PICRUSt2 clears it by "
         f"{p2nn['delta_f1_population']:.3f}."),
        ("Harder at family rank.  ", f"The gap widens to {sg(cf['delta_f1_population'], 3)} "
         f"and prFBA's skill falls to {sg(f['arms']['prfba_nestedcv']['ec']['skill_score'], 2)}, "
         "a point estimate whose interval against the prevalence floor spans zero. "
         "Clade-max is an oracle chosen against the truth, not a method to beat."),
    ], size=12.5, space=9)
    # the second leave-genus-out pass records C11/C12 and stamps the artefact invalid on C12
    bad = ", ".join(r.split(":")[-1] for r in g["headline"].get("run_invalid_reasons") or [])
    footer(s, "Cluster bootstrap over holdout groups. Every fold-scoring control passed"
              + (f"; the phenotype gate {bad} fails (next slide), so the design's rule marks "
                 "this run invalid — no number here changes" if bad else ""))
    return s


def s17_phenotype(prs, D):
    s = blank(prs)
    background(s, DEEP)
    title(s, "Where prFBA answers, ATP barely registers the gap", dark=True, tag="NEW",
          sub="Every arm's gene set rebuilt into a metabolic model, then asked the same question")
    pa, pg = D["phatp"], D["phgf"]
    A = {k: v["controlled"] for k, v in pa["arms"].items()}
    G = {k: v["controlled"] for k, v in pg["arms"].items()}
    c = pa["contrasts"]["prfba_exact_tie__vs__picrust2"]["controlled"]
    cnn = pa["contrasts"]["prfba_exact_tie__vs__base_nn1"]["controlled"]
    med = lambda blk, k: blk[k]["median"]                                  # noqa: E731
    gated = len(D["phgf_gated"])
    mcn = {m: v for m, v in pa["contrasts"]["prfba_exact_tie__vs__picrust2"]["controlled"]
           ["mcnemar_per_medium"].items() if v.get("q_bh") is not None and v["q_bh"] < 0.05}
    mcn_pic = sum(1 for v in mcn.values() if v["c"] > v["b"])      # c: PICRUSt2 right, prFBA wrong
    n_tested = sum(1 for v in pa["contrasts"]["prfba_exact_tie__vs__picrust2"]["controlled"]
                   ["mcnemar_per_medium"].values() if v.get("p") is not None)

    tf = textbox(s, M, BODY_TOP, 6.6, 4.4)
    bullets(tf, [
        ("The readout.  ", "Each arm's EC set goes through one fixed ModelSEED "
         "reconstruction, and the model is asked which of 54 media it can make ATP on — "
         "22 ever carry signal. The score is media called wrong against the truth model."),
        ("The result.  ", f"Over the {c['n_pairs']} ATP-informative genomes (prFBA answers "
         f"all {pa['meta']['n_candidates_eligible']:,} assayed), where its EC-F1 gap is only "
         f"{sg(PHENO_EC_GAP, 3)}, it calls "
         f"{abs(c['delta_atp_hamming']):.2f} "
         f"{'more' if c['delta_atp_hamming'] > 0 else 'fewer'} media wrong than PICRUSt2, "
         f"p = {c['wilcoxon_p']:.2f}: no mean difference — a null, not an equivalence. Per "
         f"medium, over all {pa['meta']['n_candidates_eligible']:,}, PICRUSt2 is ahead on "
         f"{mcn_pic} of {n_tested} "
         f"(q = {rd(max(v['q_bh'] for v in mcn.values()))})."),
        ("Against the bar.  ", f"Here prFBA beats copying the nearest reference "
         f"({abs(cnn['delta_atp_hamming']):.2f} fewer media wrong, q = {cnn['q_bh']:.3f}) — "
         "the baseline it "
         "ties on ECs over these same genomes."),
        ("Structure and essentiality.  ", "PICRUSt2's drafts share more of the truth's "
         f"reactions (median Jaccard {med(A['picrust2'], 'rxn_jaccard_draft'):.3f} against "
         f"{med(A['prfba_exact_tie'], 'rxn_jaccard_draft'):.3f}). On reaction Jaccard, "
         "essentiality and growth prFBA is indistinguishable from the nearest-reference "
         "baseline (unpaired medians; no paired test specified)."),
    ], size=12.5, space=9, color=WHITE, head_color=MINT)

    gx, gw = 7.55, 5.18
    tf2 = textbox(s, gx, BODY_TOP - 0.02, gw, 0.3)
    para(tf2, "Media called wrong, per genome (lower is better)", size=12, bold=True,
         color=MINT, first=True, space_after=0)
    rows = [("Truth, rebuilt", A["truth"]["atp_hamming"]["mean"], MINT),
            ("PICRUSt2", A["picrust2"]["atp_hamming"]["mean"], WHITE),
            ("prFBA", A["prfba_exact_tie"]["atp_hamming"]["mean"], WHITE),
            ("Nearest neighbour", A["base_nn1"]["atp_hamming"]["mean"], WHITE),
            ("Constant (floor)", A["base_constant"]["atp_hamming"]["mean"], WHITE)]
    y = BODY_TOP + 0.36
    for lab, v, col in rows:
        card = rect(s, gx, y, gw, 0.52, RGBColor(0x0A, 0x4A, 0x55))
        ctf = card.text_frame
        ctf.margin_left = ctf.margin_right = Inches(0.2)
        ctf.vertical_anchor = MSO_ANCHOR.MIDDLE
        pp = ctf.paragraphs[0]
        pp.alignment = PP_ALIGN.LEFT
        r = pp.add_run()
        r.text = f"{v:.2f}"
        r.font.size, r.font.bold, r.font.name = Pt(18), True, HEAD
        r.font.color.rgb = col
        r2 = pp.add_run()
        r2.text = "   " + lab
        r2.font.size, r2.font.name = Pt(12), BODY
        r2.font.color.rgb = WHITE
        y += 0.6
    card = rect(s, gx, y + 0.08, gw, 1.32, RGBColor(0x0A, 0x4A, 0x55))
    ctf = card.text_frame
    ctf.margin_left = ctf.margin_right = Inches(0.18)
    ctf.margin_top = Inches(0.1)
    para(ctf, "The two phenotype controls", size=12.5, bold=True, color=AMBER, font=HEAD,
         first=True, space_after=4)
    ident = (pg["meta"].get("identity_control") or {}).get("pass")
    para(ctf, f"Identity (C11): {'passes' if ident else 'FAILS'} — one EC set built twice "
              "gives identical models. Validity gate (C12): "
              + (f"fails — {gated} of {pg['meta']['n_candidates_eligible']} gapfilled truth "
                 "models make ATP or biomass from water and protons alone. Those genomes are now "
                 "dropped from every arm's gapfilled comparison (a scorer bug had kept "
                 "them)." if gated else "passes."),
         size=11, color=WHITE, space_after=0, line=1.04)
    footer(s, f"ATP: {pa['meta']['n_candidates_eligible']:,} pre-specified genomes, "
              f"Wilcoxon signed-rank, BH q = {c['q_bh']:.2f}; gapfilled readouts: "
              f"{pg['meta']['n_candidates_eligible']} genomes, unpaired per-arm summaries",
           dark=True)
    s.notes_slide.notes_text_frame.text = (
        "The EC-level loss is real and the manuscript reports it. The phenotype assay can "
        "only score genomes prFBA answers, where the EC gap is small. There, over the 755 "
        "ATP-informative genomes, we cannot detect a mean difference — a null result, on an "
        "assay that resolves 22 "
        "media — though the per-medium tests favour PICRUSt2 on two of them.")
    return s


# ------------------------------------------------ why the MAG arm exists at all
def s19_mag_why(prs, D):
    s = blank(prs)
    title(s, "Most real MAGs have no 16S to give", tag="NEW",
          sub="The amplicon route is not available for the genomes people actually assemble")
    ssu = D["ssu"]["by_kind"]
    cd = CategoryChartData()
    cd.categories = ["Carries any SSU gene", "Carries SSU ≥ 1,200 bp"]
    cd.add_series("Isolate genomes", [ssu["isolate"]["has_ssu"],
                                      ssu["isolate"]["has_ssu_ge1200"]], number_format="0%")
    cd.add_series("MAGs", [ssu["MAG"]["has_ssu"], ssu["MAG"]["has_ssu_ge1200"]],
                  number_format="0%")
    gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(BODY_TOP + 0.32),
                            Inches(6.1), Inches(3.3), cd)
    ch = gf.chart
    style_chart(ch, size=10.5)
    series_color(ch.series[0], TEAL)
    series_color(ch.series[1], AMBER)
    ch.value_axis.maximum_scale = 1.0
    ch.value_axis.minimum_scale = 0.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    for sr in ch.series:
        sr.has_data_labels = True
        dl = sr.data_labels
        dl.show_value = True
        dl.number_format, dl.number_format_is_linked = "0%", False
        dl.position = XL_LABEL_POSITION.OUTSIDE_END
        dl.font.size = Pt(9)
        dl.font.color.rgb = MUTED
        dl.font.name = BODY
    tfc = textbox(s, M, BODY_TOP, 6.1, 0.3)
    para(tfc, f"GTDB r232, {D['ssu']['n_genomes']:,} genomes", size=11.5, bold=True,
         color=MUTED, first=True, space_after=0)

    ret = {r["completeness"]: r for r in D["ret"]["clean"]}
    tf = textbox(s, 7.0, BODY_TOP, 5.73, 3.5)
    bullets(tf, [
        ("The gap.  ", f"{ssu['MAG']['has_ssu_ge1200'] * 100:.1f}% of "
         f"{ssu['MAG']['n']:,} real MAGs carry a full-length 16S, against "
         f"{ssu['isolate']['has_ssu_ge1200'] * 100:.1f}% of isolates. A 16S method "
         "simply cannot be run on most of them."),
        ("Our synthetic MAGs are kinder.  ", f"Full-length 16S survives in "
         f"{ret[1.0]['full'] * 100:.1f}% of complete drafts and "
         f"{ret[0.3]['full'] * 100:.1f}% at 30% completeness — about four times the rate in "
         f"real MAGs of similar completeness "
         f"({D['ssu']['by_kind_completeness']['MAG']['95-100']['has_ssu_ge1200'] * 100:.1f}% "
         "at 95–100%) — so they understate the problem."),
        ("And contigs are not a substitute.  ", "Fed a MAG directly, PICRUSt2 refuses "
         "the default run outright. Forced through, it treats each contig as its own "
         "organism."),
    ], size=13, space=11)

    fd = D["feed"]
    card = rect(s, M, 5.12, 12.13, 1.68, LIGHT)
    ctf = card.text_frame
    ctf.margin_left = ctf.margin_right = Inches(0.2)
    ctf.margin_top = Inches(0.16)
    para(ctf, "What “forced through” looks like, on one bacterial genome",
         size=13.5, bold=True, color=INK, font=HEAD, first=True, space_after=8)
    para(ctf, f"Its {fd['n_contigs']} contigs draw {fd['distinct_closest_refs']} different "
              f"closest references and {fd['domains']['arc']} of {fd['n_contigs']} are "
              f"called archaeal. The union asserts {fd['union_over_contigs']['n_ec']:,} EC "
              f"numbers against {fd['n_truth_ec']} real ones — precision "
              f"{fd['union_over_contigs']['precision']:.3f}. "
              f"{fd['no16s_contigs']['n']} of the {fd['n_contigs']} contigs carry no 16S, "
              f"yet {fd['no16s_contigs']['nsti_lt_0_15']} of them are reported at NSTI "
              "< 0.15, which reads as a confident placement.",
         size=12, color=MUTED, space_after=0, line=1.06)
    footer(s, "The one full-length 16S that genome does carry scores precision "
              f"{fd['fulllength16s_amplicon']['precision']:.3f} — the marker works when it "
              "is there; the problem is that usually it is not")
    return s


def s20_mag_build(prs, D):
    s = blank(prs)
    title(s, "How a synthetic MAG is built", tag="NEW",
          sub="Two-step degradation, because dropping deposited contigs would confound with assembly quality")
    steps = [
        ("1", "Fragment", f"cut each deposited\ncontig to ~{FRAG_MEAN_KB} kb mean,\n"
         f"{FRAG_MIN_KB} kb floor"),
        ("2", "Drop", "remove whole fragments\nto hit the target\ncompleteness"),
        ("3", "Contaminate", "splice in fragments\nfrom a different-phylum\ndonor"),
        ("4", "Label", "prefix every contig\nheader with the MAG id"),
        ("5", "Seed", "RNG from a CRC32 of\nthe genome id"),
    ]
    x, w, gap = M, 2.186, 0.3
    for i, (n, head, desc) in enumerate(steps):
        card = rect(s, x, BODY_TOP + 0.2, w, 2.4, LIGHT)
        tf = card.text_frame
        tf.margin_left = tf.margin_right = Inches(0.16)
        tf.margin_top = Inches(0.74)
        para(tf, head, size=15.5, bold=True, color=INK, font=HEAD, first=True, space_after=7)
        para(tf, desc, size=12, color=MUTED, space_after=0, line=1.05)
        dot = rect(s, x + 0.16, BODY_TOP - 0.07, 0.54, 0.54, TEAL, shape=MSO_SHAPE.OVAL)
        dtf = dot.text_frame
        dtf.margin_left = dtf.margin_right = dtf.margin_top = dtf.margin_bottom = 0
        dtf.vertical_anchor = MSO_ANCHOR.MIDDLE
        para(dtf, n, size=16, bold=True, color=WHITE, font=HEAD,
             align=PP_ALIGN.CENTER, first=True, space_after=0)
        x += w + gap

    tf = textbox(s, M, BODY_TOP + 2.82, 6.4, 2.3)
    bullets(tf, [
        ("Why fragment first.  ", f"{PANEL_CLOSED_PCT}% of panel genomes are closed and "
         f"{PANEL_LT10_PCT}% have fewer than ten contigs, so dropping deposited contigs "
         "cannot reach a 10%-granular grid and would make degradation severity depend on "
         "each reference's assembly quality — across a panel spanning 1 to 5,265 contigs."),
        ("Bases are never mutated.  ", "skani compares k-mer sketches, and what a MAG "
         "actually suffers is incompleteness, not divergence."),
    ], size=13, space=10)

    card = rect(s, 7.2, BODY_TOP + 2.82, 5.53, 1.08, LIGHT)
    ctf = card.text_frame
    ctf.margin_left = ctf.margin_right = Inches(0.2)
    ctf.margin_top = Inches(0.14)
    para(ctf, "Step 4 is not cosmetic", size=13, bold=True, color=INK, font=HEAD,
         first=True, space_after=6)
    para(ctf, "skani keys a query on its first contig header and silently collapses "
              "collisions: eight variants sharing one header returned five queries.",
         size=12, color=MUTED, space_after=0, line=1.05)

    cr = D["ret"]["c100_reasons"]
    card2 = rect(s, 7.2, BODY_TOP + 4.04, 5.53, 1.08, LIGHT)
    c2 = card2.text_frame
    c2.margin_left = c2.margin_right = Inches(0.2)
    c2.margin_top = Inches(0.14)
    para(c2, "What the 2.5 kb floor costs", size=13, bold=True, color=INK, font=HEAD,
         first=True, space_after=6)
    para(c2, f"At full completeness {cr['16s_only_on_contigs_dropped_by_min_frag']:,} MAGs "
             f"lose their only 16S to it, and {cr['genome_has_no_16s']} source genomes "
             "never had one.", size=12, color=MUTED, space_after=0, line=1.05)
    footer(s, "mag_validation/synth_mags.py — a given (genome, completeness, contamination) "
              "is byte-identical across runs and machines")
    return s



# ------------------------------------------- the synthetic-genome step, re-examined
def s06b_weighting(prs, D):
    s = blank(prs)
    title(s, "What the weighting can — and cannot — claim", tag="NEW",
          sub="A literature and code review of the synthetic-genome step; ten defects survived adversarial review")
    for i, (v, lab, col) in enumerate([
            (f"{BENCH_N1}%", "of answered benchmark queries\nselect one genome: p = 1 for all", TEAL),
            (f"{DEPLOY_N1}%", "of answered ASVs in the 3,950-ASV\ndeployment run do the same", TEAL),
            ("0", "experiments in this work that\nexercise the weighting", AMBER_TXT),
            (f"{XPHYLUM_DEPLOY[1]:.0f}%", "of 3+-genome merges in that run\nspan different phyla", AMBER_TXT)]):
        stat(s, M + i * 3.108, BODY_TOP, 2.808, v, lab, color=col)
    tf = textbox(s, M, BODY_TOP + 1.7, 7.4, 3.6)
    bullets(tf, [
        ("Degenerate where it matters.  ", "With one reference, 1 − p^exp vanishes and "
         "prFBA reduces to a two-level parsimony objective. Where the weight is active, "
         "its commonest fractional value is 0.5."),
        ("It validates the selector, not the weight.  ", "Every arm in this deck predicts "
         "the unweighted union of the selected references' content. No reported number "
         "tests the weighting."),
        ("A prevalence, not a probability.  ", "Genomes at the species boundary share only "
         "65–90% of their genes, so p = 1 overstates presence — and absence means no annotation, "
         "not no gene."),
        ("Allocated inversely to need.  ", "Selection size tracks database redundancy: "
         "largest where one reference suffices, one where the query is novel."),
    ], size=13, space=10)
    card = rect(s, 8.3, BODY_TOP + 1.7, 4.43, 3.2, LIGHT)
    ctf = card.text_frame
    ctf.margin_left = ctf.margin_right = Inches(0.2)
    ctf.margin_top = Inches(0.18)
    para(ctf, "Prior art it re-derives", size=14, bold=True, color=INK, font=HEAD,
         first=True, space_after=9)
    for it in ("Likelihood-based gapfilling on the same ModelSEED/KBase stack (Benedict 2014)",
               "The objective is GIMME's, with prevalence for expression (Becker & Palsson 2008)",
               "Pangenome profiles from a 16S call (PanFP); nearest-genome copying (Piphillin)",
               "New is the combination: exact-tie selection at database depth, plus abstention"):
        hanging(para(ctf, "•  " + it, size=12, color=MUTED, space_after=8, line=1.04))
    footer(s, f"Tie sets on the controlled benchmark span phyla {XPHYLUM_BENCH}% of the time; "
              f"in the 3,950-ASV deployment run {XPHYLUM_DEPLOY[0]}% (two genomes) to "
              f"{XPHYLUM_DEPLOY[1]}% (three or more)")
    return s


# ------------------------------------------------------------ the random panel
def s21_random_panel(prs, D):
    s = blank(prs)
    title(s, "A third panel, drawn at random", tag="NEW",
          sub="The null composition: uniform over the same 140,197 admissible organisms the k-means panel clusters")
    C = {k: composition(v) for k, v in D["sel3"].items()}
    rows = [("Archaea", lambda c: f"{c['archaea']:,} ({c['arc_pct']:.1f}%)"),
            ("Phyla", lambda c: f"{c['phyla']}"),
            ("Named families", lambda c: f"{c['families']:,}"),
            ("Named genera", lambda c: f"{c['genera']:,}"),
            ("Species", lambda c: f"{c['species']:,}"),
            ("Genus Gini", lambda c: f"{c['gini']:.3f}"),
            ("Largest genus", lambda c: f"{c['largest']}")]
    hdr = ["", "Stratified", "k-means", "Random"]
    tbl = s.shapes.add_table(len(rows) + 1, 4, Inches(M), Inches(BODY_TOP),
                             Inches(6.6), Inches(0.42 * (len(rows) + 1))).table
    for j, wid in enumerate((2.1, 1.5, 1.5, 1.5)):
        tbl.columns[j].width = Inches(wid)
    for i in range(len(rows) + 1):
        for j in range(4):
            cell = tbl.cell(i, j)
            cell.fill.solid()
            cell.fill.fore_color.rgb = DEEP if i == 0 else (LIGHT if i % 2 else WHITE)
            text = hdr[j] if i == 0 else (rows[i - 1][0] if j == 0
                                          else rows[i - 1][1](C[("original", "kmeans", "random")[j - 1]]))
            tf = cell.text_frame
            tf.margin_left = tf.margin_right = Inches(0.1)
            p = tf.paragraphs[0]
            p.alignment = PP_ALIGN.LEFT if j == 0 else PP_ALIGN.CENTER
            r = p.add_run()
            r.text = text
            r.font.size, r.font.name = Pt(13), BODY
            r.font.bold = (i == 0) or (j == 3)
            r.font.color.rgb = WHITE if i == 0 else (SEA_TXT if j == 3 else INK)
    tf = textbox(s, 7.6, BODY_TOP, 5.13, 4.9)
    bullets(tf, [
        ("The draw.  ", "Uniform over organisms, not sequences (seed 1729), with the "
         "original panel's representative-16S rule. Admissibility is identical to k-means, "
         "so any difference is the sampling scheme alone."),
        ("The null.  ", f"No correction in either direction: Archaea fall to "
         f"{C['random']['arc_pct']:.1f}% against 3.5% of the pool, and one genus "
         f"contributes {C['random']['largest']} organisms."),
        ("Run through everything.  ", "The full eight-region 16S arm, the panel-as-reference "
         "MAG arm at depth 500, and the GTDB arm on a 20% phylum-proportional subsample."),
        ("Checked before trusted.  ", "Each comparison's crude column reproduces the panel's "
         "own scorer exactly — 288 of 288 amplicon cells, 30 of 30 GTDB cells."),
    ], size=13, space=11)
    footer(s, "Genus Gini over genus-level bins, placeholders included; overlap 804 "
              "stratified–k-means, 524 stratified–random, 701 k-means–random, 49 in all three")
    return s


def s22_three_panel(prs, D):
    s = blank(prs)
    title(s, "Coverage falls on the random panel; quality holds", tag="NEW",
          sub="Phylum-standardized over the 113 phyla in all three panels, self-exclusion at τ = 0.987")
    regs = ["FullLength16S", "V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V4-V5_944R", "V7-V9"]
    cd = CategoryChartData()
    cd.categories = [REG_LABEL[r] for r in regs]
    for key, lab in (("original", "Taxonomy-stratified"), ("kmeans", "k-means"),
                     ("random", "Uniform random")):
        cd.add_series(lab, [tp3_cell(D, r, "exclude_0.987", "genus_coverage")[key]["phystd"]
                            for r in regs], number_format="0.00")
    gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(BODY_TOP + 0.3),
                            Inches(7.9), Inches(4.5), cd)
    ch = gf.chart
    style_chart(ch, size=10.5)
    for sr, col in zip(ch.series, (TEAL, AMBER, SEAFOAM)):
        series_color(sr, col)
    ch.value_axis.maximum_scale = 0.8
    ch.value_axis.minimum_scale = 0.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    tfc = textbox(s, M, BODY_TOP, 7.9, 0.3)
    para(tfc, "Genus coverage of novel organisms, by region", size=11.5, bold=True,
         color=MUTED, first=True, space_after=0)
    cov = [tp3_cell(D, r, "exclude_0.987", "genus_coverage")["delta_random_minus_original"]["est"] for r in regs]
    acc = [abs(tp3_cell(D, r, "exclude_0.987", "genus_cond_acc")["delta_random_minus_original"]["est"]) for r in regs]
    f1 = [abs(tp3_cell(D, r, "exclude_0.987", "genecap_cond_f1")["delta_random_minus_original"]["est"]) for r in regs]
    x = 8.85
    for i, (v, lab, col) in enumerate([
            (f"{sg(max(cov), 2)} to {sg(min(cov), 2)}", "random − stratified coverage,\nevery region, CI excludes 0", AMBER_TXT),
            (f"≤ {max(acc):.2f}", "difference in conditional genus\naccuracy, either direction", TEAL),
            (f"≤ {max(f1):.3f}", "difference in conditional\ngene-capture F1", TEAL)]):
        stat(s, x, BODY_TOP + i * 1.6, 3.88, v, lab, color=col, h=1.42)
    footer(s, "Standardization keeps 77–91% of the crude coverage gap, so composition is a minority of it; "
              "1,000 phylum-stratified bootstrap replicates")
    return s


def s23_reversal(prs, D):
    s = blank(prs)
    title(s, "Against a dense reference the sign reverses", tag="NEW",
          sub="The coverage gap belongs to the 20-hit record, not to the organisms")
    v4 = tp3_cell(D, "V4", "exclude_0.987", "genus_coverage")
    gt = next(c for c in D["tp3g"]["cells"] if c["regime"] == "exclude" and c["completeness"] == 1.0)
    cd = CategoryChartData()
    cd.categories = ["16S V4 amplicon, genus coverage", "Whole-genome MAG vs GTDB, resolve rate"]
    for key, lab in (("original", "Taxonomy-stratified"), ("kmeans", "k-means"),
                     ("random", "Uniform random")):
        cd.add_series(lab, [v4[key]["phystd"], gt[key]["phystd"]], number_format="0.00")
    gf = s.shapes.add_chart(XL_CHART_TYPE.COLUMN_CLUSTERED, Inches(M), Inches(BODY_TOP + 0.3),
                            Inches(6.4), Inches(4.5), cd)
    ch = gf.chart
    style_chart(ch, size=10.5)
    for sr, col in zip(ch.series, (TEAL, AMBER, SEAFOAM)):
        series_color(sr, col)
        sr.has_data_labels = True
        dl = sr.data_labels
        dl.show_value = True
        dl.number_format, dl.number_format_is_linked = "0.00", False
        dl.position = XL_LABEL_POSITION.OUTSIDE_END
        dl.font.size = Pt(9.5)
        dl.font.color.rgb = MUTED
        dl.font.name = BODY
    ch.value_axis.maximum_scale = 0.85
    ch.value_axis.minimum_scale = 0.0
    ch.value_axis.has_major_gridlines = True
    ch.value_axis.major_gridlines.format.line.color.rgb = RGBColor(0xDC, 0xE6, 0xE6)
    ch.value_axis.major_gridlines.format.line.width = Pt(0.75)
    ch.category_axis.has_major_gridlines = False
    tfc = textbox(s, M, BODY_TOP, 6.4, 0.3)
    para(tfc, "Novel organisms, phylum-standardized", size=11.5, bold=True, color=MUTED,
         first=True, space_after=0)
    dg = gt["delta_random_minus_original"]
    tf = textbox(s, 7.35, BODY_TOP, 5.38, 5.1)
    bullets(tf, [
        ("The flip.  ", f"On the amplicon arm the random panel answers less; against GTDB "
         f"it resolves more ({sg(dg['est'], 2)}, 95% CI {sg(dg['ci'][0], 2)} to {sg(dg['ci'][1], 2)}), "
         "holding down to 30% completeness."),
        ("Why.  ", f"On full-length 16S, exclusion empties a full 20-hit record for "
         f"{EMPTIED['original']}%, {EMPTIED['kmeans']}% and {EMPTIED['random']}% of "
         "queries — a truncation artefact; a full 100-hit GTDB record empties for under 1%. "
         "And representative organisms keep closer relatives: median "
         f"best surviving ANI {GTDB_BEST_ANI['original']} → {GTDB_BEST_ANI['random']}."),
        ("Self-inclusion too.  ", f"The source is in its own top 20 for "
         f"{SELF_TOP20['original']}%, {SELF_TOP20['kmeans']}% and {SELF_TOP20['random']}% "
         "of amplicons — misses split between exact-tie eviction and the index crediting "
         "the sequence to another organism."),
        ("What it licenses.  ", "Our headline self-exclusion coverage is a property of the "
         "stratified panel at depth 20. Conditional quality is the claim that replicates."),
    ], size=12.5, space=9)
    footer(s, "GTDB r226: 143,614 species-representative genomes; exclusion at ANI ≥ 95; "
              "~2,000-genome phylum-proportional subsample of each panel")
    return s


# --------------------------------------------------------------- the audit
def s24_audit(prs, D):
    s = blank(prs)
    title(s, "How the work was checked", tag="NEW",
          sub="Every falsifiable claim in the revised and new sections re-derived from the data; a second auditor re-checked each flag")
    for i, (v, lab, col) in enumerate([
            (f"{sum(AUDIT_CLAIMS)}", "manuscript claims audited\nagainst the artefacts", TEAL),
            (f"{sum(AUDIT_DEFECTS)}", "defects confirmed on\nrecheck, all corrected", AMBER_TXT),
            ("288/288", "amplicon cells reproduced\nbefore any comparison", TEAL),
            (f"{len(PIPELINE_CAUGHT)}", "pipeline and analysis\nfaults caught and fixed", AMBER_TXT)]):
        stat(s, M + i * 3.108, BODY_TOP, 2.808, v, lab, color=col)
    def col(x, head, items):
        card = rect(s, x, BODY_TOP + 1.7, 5.92, 3.4, LIGHT)
        ctf = card.text_frame
        ctf.margin_left = ctf.margin_right = Inches(0.2)
        ctf.margin_top = Inches(0.18)
        para(ctf, head, size=14, bold=True, color=INK, font=HEAD, first=True, space_after=9)
        for it in items:
            hanging(para(ctf, "•  " + it, size=12, color=MUTED, space_after=7, line=1.04))
    col(M, "What the claim audits caught", [
        "The searched database: ~315,000 genomes, not 112,000 — nor the 268,000 a first "
        "audit substituted",
        "GTDB's sketch holds 143,614 species representatives, not 732,475 accessions",
        "prFBA's objective is probabilistic; it was described as protein-constrained",
        "A tenfold pangenome finding cited to Roary belongs to Panaroo",
        "Cross-phylum merge rates that counted an unknown phylum as a phylum (40.6% → 31.7%)",
        "A phenotype scorer that compared arms against gate-failing truth models, and a "
        "summary that counted duplicated rows twice"])
    col(M + 6.21, "What our own runs and reproduction checks caught", PIPELINE_CAUGHT)
    footer(s, f"This deck and the manuscript update were audited the same way — {DECK_AUDIT[0]} "
              f"claims, {DECK_AUDIT[1]} confirmed defects; then {FIX_AUDIT[0]} claims over the "
              f"fixes, {FIX_AUDIT[1]} more — all fixed. Gene-family cache intact against its "
              "pre-run snapshot")
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=None,
                    help="default presentations/prFBA_updates_2026-10.pptx "
                         "(prFBA_updates_2026-10_d100.pptx at --depth 100)")
    ap.add_argument("--depth", type=int, default=20, choices=[20, 100],
                    help="amplicon record depth whose results the deck reports (default 20)")
    ap.add_argument("--rule", default="legacy", choices=["legacy", "union"],
                    help="--depth 100 only: selection rule of the original-panel and "
                         "three-panel inputs (default legacy, the rule tab:main uses)")
    ap.add_argument("--constants", default=None,
                    help="--depth 100 only: JSON with SELF_TOP20/SELF_SAT/SELF_CRED/EMPTIED")
    args = ap.parse_args()
    configure_depth(args.depth, args.rule, args.constants)
    if args.out is None:
        args.out = os.path.join(REPO, "presentations", "prFBA_updates_2026-10%s.pptx"
                                % ("" if args.depth == 20 else "_d%d" % args.depth))
    D = load()
    prs = Presentation()
    prs.slide_width, prs.slide_height = Inches(W), Inches(H)
    for fn in (s01_title, s02_problem, s03_pipeline, s04_alignment, s05_selection,
               s06_synthetic, s06b_weighting, s07_validation, s08_regions,
               s09_denominators, s10_benchmark_design, s11_results, s12_novelty,
               s13_ceiling, s13b_dialect, s15_kfold_design, s16_kfold_result,
               s17_phenotype, s14_depth, s19_mag_why, s20_mag_build, s16_mag_arm,
               s17_crossover, s18_abstain, s19_two_panel, s21_random_panel,
               s20_panel_space, s22_three_panel, s23_reversal, s24_audit, s15_next):
        fn(prs, D)
    prs.core_properties.title = "prFBA project update"
    prs.core_properties.author = "Andrew Freiburger"
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    prs.save(args.out)
    print(f"[make_deck] {len(prs.slides.__iter__.__self__._sldIdLst)} slides -> {args.out}")


if __name__ == "__main__":
    main()
