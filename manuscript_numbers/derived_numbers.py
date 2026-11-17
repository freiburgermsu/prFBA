#!/usr/bin/env python
"""Derived manuscript numbers that no committed artefact stores (rescore plan item W6).

Each number in this file was printed in prFBA_paper/prFBA.tex (or carried by hand in
presentations/make_deck.py) from an ad hoc recomputation over committed cell tables. This
script is that recomputation, written down once, so that a re-score can refresh the
numbers mechanically.

    python manuscript_numbers/derived_numbers.py --state pre    # control: pre-repair inputs
    python manuscript_numbers/derived_numbers.py --state post   # the live tree

--state pre resolves every repo-relative input as
    $W/pre_repair/<rel>   if it exists (old outputs moved aside by the re-score),
    $W/start/<rel>        if <rel> is tracked at the start tag (rescore-2026-11-start),
    $R/<rel>              otherwise (git-ignored and not yet moved aside), with a warning.
--state post reads the live tree. Every input is recorded with its sha256 in
_meta.inputs, and every number with the manuscript string it reproduces (`printed`), the
computed value, the formatted value, and whether the two match. On --state pre each
mismatch carries a `mismatch_reason`.

Writes manuscript_numbers/derived_numbers.<state>.json and prints a side-by-side table.

Recipes (sources in brackets):
  sr_*        controlled benchmark selection sizes and tie-set coherence [cells.csv.gz
              arm prfba_exact_tie/prfba_legacy, density 1.0, space ec, n_reps > 0, one row per
              (md5, fold); tie sets from prfba/f0{0..4}_d100/predictions.json exact_tie reps,
              lineages from panel_eligible.json; first computed at prFBA_paper 3f95f1c]
  recon_modal_weight   modal fractional k/n feature weight in multi-genome syntheses
  as10_*      matched both-answer EC F1, as shipped [D100_VS_D20.md section 2D]
  magvs16s_*  tab:magvs16s and the MAG-vs-16S paired claims [UNION_THROUGHOUT_CHANGES.md
              items 1.1, 1.48-1.61; compare_arms.py / crossover.py pairing logic]
  mag_band_*  identity-bin answer rates [item 1.56; RESULTS_v2.md]
  mag_af_*    ANI-floor / AF-gate sweep [make_deck.py AF_SWEEP; af_gate_counterfactual.py]
  mag_reach_* marker reach [RESULTS_v2.md "Why the MAG arm abstains"]
  panels_*    Welch tests, MAG singleton abstention and multiplicity, pooled surviving-
              uncovered rates [items 1.40, 1.72, 1.80; two_panel/compare_mag_panels.py]
"""
import argparse, collections, csv, datetime, fnmatch, glob, gzip, hashlib, json, math, os, re
import statistics, subprocess, sys

R_DEFAULT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
W_DEFAULT = os.environ.get("PRFBA_RESCORE_W",
                           os.path.expanduser("~/Documents/prFBA_rescore_2026-11"))
TAG = "rescore-2026-11-start"
PAPER = os.path.expanduser("~/Documents/prFBA_paper")
PAPER_PRE_COMMIT = "0214a14"          # manuscript state the printed values are taken from
EMILYKIN_SHIPPED_D100 = os.path.expanduser("~/Documents/EmilyKin/synthetic_genomes_d100")
KN_EMILYKIN_BVBRC = "deployment_d100/downstream_e3/kn_builds/emilykin"

REGIONS = ("FullLength16S", "V1-V2", "V1-V3", "V3-V4", "V4", "V4-V5", "V4-V5_944R", "V7-V9")
TEX_REGION = {"FullLength16S": "Full-length", "V1-V2": "V1--V2", "V1-V3": "V1--V3",
              "V3-V4": "V3--V4", "V4": "V4", "V4-V5": "V4--V5", "V4-V5_944R": "V4--V5 944R",
              "V7-V9": "V7--V9"}


# ============================================================================ input resolution
class Resolver:
    """Maps repo-relative input paths to files for one state and records what was read."""

    def __init__(self, state, repo, w, tag=TAG, tracked=None):
        self.state, self.repo, self.w, self.tag = state, repo, w, tag
        self.start = os.path.join(w, "start")
        self.pre_repair = os.path.join(w, "pre_repair")
        self._tracked = tracked
        self.inputs = {}
        self.warnings = []
        self._sha = {}

    @property
    def tracked(self):
        if self._tracked is None:
            out = subprocess.run(["git", "-C", self.repo, "ls-tree", "-r", "--name-only", self.tag],
                                 check=True, capture_output=True, text=True).stdout
            self._tracked = set(out.splitlines())
        return self._tracked

    def _record(self, rel, path, source, warn=None, sha=True):
        rec = {"path": path, "source": source}
        if sha:
            rec["sha256"] = self.sha256(path)
        if warn:
            rec["warning"] = warn
            self.warnings.append("%s: %s" % (rel, warn))
        self.inputs[rel] = rec

    def sha256(self, path):
        if path not in self._sha:
            h = hashlib.sha256()
            with open(path, "rb") as fh:
                for b in iter(lambda: fh.read(64 << 20), b""):
                    h.update(b)
            self._sha[path] = h.hexdigest()
        return self._sha[path]

    def path(self, rel, sha=True):
        if self.state == "post":
            p, src, warn = os.path.join(self.repo, rel), "live", None
        else:
            pr = os.path.join(self.pre_repair, rel)
            if os.path.exists(pr):
                p, src, warn = pr, "pre_repair", None
            elif rel in self.tracked:
                p, src, warn = os.path.join(self.start, rel), "start", None
            else:
                p, src = os.path.join(self.repo, rel), "live"
                warn = ("git-ignored and not in $W/pre_repair: read from the live tree, "
                        "valid as a pre-repair input only while the re-score has not rewritten it")
        if not os.path.exists(p):
            raise SystemExit("[derived] missing input %s (%s)" % (rel, p))
        self._record(rel, p, src, warn, sha)
        return p

    def glob(self, pattern, sha=True):
        """Repo-relative glob, resolved as a whole from one source."""
        if self.state == "post":
            bases = [("live", self.repo)]
        else:
            bases = [("pre_repair", self.pre_repair)]
            if any(fnmatch.fnmatch(t, pattern) for t in self.tracked):
                bases.append(("start", self.start))
            bases.append(("live", self.repo))
        for src, base in bases:
            hits = sorted(glob.glob(os.path.join(base, pattern)))
            if hits:
                warn = None
                if self.state == "pre" and src == "live":
                    warn = ("git-ignored and not in $W/pre_repair: read from the live tree, valid "
                            "as a pre-repair input only while the re-score has not rewritten it")
                for p in hits:
                    self._record(os.path.relpath(p, base), p, src, warn, sha)
                return hits
        raise SystemExit("[derived] no input matches %s" % pattern)

    def external_dir(self, label, path, digest):
        rec = {"path": path, "source": "external", "digest_sha256": digest,
               "digest_of": "sha256 over sorted '<file name> <file sha256>' lines"}
        self.inputs[label] = rec


# ============================================================================ formatting
def _fmt(kind, v):
    if v is None:
        return None
    if kind == "int":
        return "{:,}".format(int(v))
    if kind == "f0":
        return "%.0f" % v
    if kind == "f1":
        return "%.1f" % v
    if kind == "f2":
        return "%.2f" % v
    if kind == "f3":
        return "%.3f" % v
    if kind == "f4":
        return "%.4f" % v
    if kind == "sf2":
        return "%+.2f" % v
    if kind == "sf3":
        return "%+.3f" % v
    if kind == "pct0":
        return "%.0f" % (100.0 * v)
    if kind == "pct1":
        return "%.1f" % (100.0 * v)
    if kind == "e1":
        return "%.1e" % v
    if kind == "e0":
        return "%.0e" % v
    if kind == "bool":
        return "yes" if v else "no"
    if kind == "str":
        return str(v)
    raise ValueError(kind)


def number(key, value, kind, printed, tex_lines, snippet, recipe, inputs, deck=None, extra=None):
    """One reported number. `printed` is the manuscript (or deck) value in the formatted
    form; `snippet` is the literal text expected in prFBA.tex (None for deck-only)."""
    f = _fmt(kind, value)
    d = dict(key=key, tex_lines=tex_lines, deck_constant=deck, printed=printed, value=value,
             formatted=f, match=(f == printed) if printed is not None else None, recipe=recipe,
             inputs=sorted(set(inputs)),
             tex_snippet=snippet)
    if extra:
        d["extra"] = extra
    return d


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def _mean(xs):
    xs = list(xs)
    return sum(xs) / len(xs) if xs else None


def _median(xs):
    s = sorted(xs)
    n = len(s)
    if not n:
        return None
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


# ============================================================================ pure recipes
def selection_size_stats(records):
    """records: iterable of (md5, fold, n_reps); answered = n_reps > 0, one per (md5, fold)."""
    d = {}
    for md5, fold, n in records:
        n = int(n or 0)
        if n > 0:
            d.setdefault((md5, str(fold)), n)
    n = len(d)
    c = collections.Counter(d.values())
    f0 = [v for (m, f), v in d.items() if f == "0"]
    return dict(answered=n, single=c[1], single_share=c[1] / n if n else None,
                mean=sum(d.values()) / n if n else None,
                fold0_answered=len(f0), fold0_single=sum(1 for v in f0 if v == 1),
                fold0_two=sum(1 for v in f0 if v == 2), fold0_three_plus=sum(1 for v in f0 if v >= 3),
                by_key=d)


def tie_set_coherence(rep_sets, lineage):
    """rep_sets: list of rep tuples with >= 2 genomes (one per (md5, fold) selection).
    lineage: genome -> {'phylum':..,'genus':..}. A set spans phyla if its members carry >= 2
    distinct non-empty phyla; it shares a genus if every member has the same non-empty genus."""
    def lin(g, r):
        return ((lineage.get(g) or lineage.get(str(g).removeprefix("g_")) or {}).get(r)) or None

    def spans(s):
        return len({p for p in (lin(g, "phylum") for g in s) if p}) >= 2

    def one_genus(s):
        ge = {lin(g, "genus") for g in s}
        return len(ge) == 1 and None not in ge

    recs = [frozenset(s) for s in rep_sets]
    uniq = sorted(set(recs), key=lambda s: sorted(s))
    return dict(n_records=len(recs), n_unique=len(uniq),
                span_records=sum(map(spans, recs)) / len(recs) if recs else None,
                span_unique=sum(map(spans, uniq)) / len(uniq) if uniq else None,
                genus_unique=sum(map(one_genus, uniq)) / len(uniq) if uniq else None)


def pair_mag_amp(mag, amp, region):
    """compare_arms.paired() for metric f1. mag: gid -> (outcome, has_precision, f1);
    amp: (gid, region) -> (outcome, f1, ...). Returns conditional and population stats plus
    each arm's coverage over the population pairs."""
    cond, pop, mo, ao = [], [], 0, 0
    for g, (mout, mprec, mf) in mag.items():
        a = amp.get((g, region))
        if a is None or mout == "truth_empty" or a[0] == "truth_empty":
            continue
        mok = mout == "ok" and mprec
        aok = a[0] == "ok"
        if mok and aok and mf is not None and a[1] is not None:
            cond.append((mf, a[1]))
        mz = mf if mok else 0.0
        az = a[1] if aok else 0.0
        if mz is not None and az is not None:
            pop.append((mz, az))
            mo += bool(mok)
            ao += bool(aok)
    cm, ca = _mean(a for a, _ in cond), _mean(b for _, b in cond)
    return dict(n=len(cond), mag=cm, amp=ca,
                delta=(cm - ca) if cond else None,
                median_delta=_median([a - b for a, b in cond]),
                n_pop=len(pop), pop_mag=_mean(a for a, _ in pop), pop_amp=_mean(b for _, b in pop),
                cov_mag=mo / len(pop) if pop else None, cov_amp=ao / len(pop) if pop else None)


def _taxon(g):
    try:
        return int(str(g).split(".")[0])
    except ValueError:
        return None


def exclude_hits(hits, gid, tau=95.0):
    """score_mag.exclude_hits: drop the source, same-taxon references and hits at ANI >= tau."""
    t = _taxon(gid)
    return [h for h in hits if h["ref"] != gid and (t is None or _taxon(h["ref"]) != t)
            and h["ani"] < tau]


AF_VARIANTS = (  # (deck label, ANI floor, AF gate); None = any surviving hit
    ("as shipped — ANI ≥ 83, AF ≥ 15%", 83.0, 15.0),
    ("ANI floor 83 → 0", 0.0, 15.0),
    ("AF gate 15% → 5%", 83.0, 5.0),
    ("AF gate 15% → 2%", 83.0, 2.0),
    ("both relaxed — ANI ≥ 80, AF ≥ 2%", 80.0, 2.0),
    ("any surviving hit at all (ceiling)", None, None),
)


def af_sweep(scoreable, hits_by_mag, gid_by_mag):
    """Share of scoreable organisms (c100_x0, exclude) with a qualifying surviving hit."""
    out = []
    for label, floor, gate in AF_VARIANTS:
        k = 0
        for m in scoreable:
            hh = exclude_hits(hits_by_mag.get(m, []), gid_by_mag[m])
            if floor is None:
                k += bool(hh)
            else:
                k += any(h["af_query"] >= gate and h["ani"] >= floor for h in hh)
        out.append(dict(label=label, floor=floor, gate=gate, k=k, n=len(scoreable),
                        share=k / len(scoreable) if scoreable else None))
    return out


def marker_reach(truth_lin, hits_by_gid):
    """truth_lin: panel genome -> truth lineage (benchmark_truth). hits_by_gid: genome ->
    c100_x0 skani hits (only genomes with a synthetic MAG). Closest panel relative = the
    lowest rank whose (non-empty) taxon another panel genome shares."""
    ranks = ("Genus", "Family", "Order", "Class", "Phylum", "Kingdom")
    cnt = {r: collections.Counter(v.get(r) for v in truth_lin.values()) for r in ranks}

    def closest(g):
        for r in ranks:
            v = truth_lin[g].get(r)
            if v and cnt[r][v] >= 2:
                return r
        return None

    cl = {g: closest(g) for g in truth_lin}
    cong = [g for g in truth_lin if cl[g] == "Genus"]
    cong_mag = [g for g in cong if g in hits_by_gid]
    best_af, nohit = [], 0
    for g in cong_mag:
        hh = exclude_hits(hits_by_gid[g], g)
        if not hh:
            nohit += 1
            continue
        best_af.append(max(h["af_query"] for h in hh))   # best = highest aligned fraction

    def nohit_share(sel):
        G = [g for g in truth_lin if sel(cl[g]) and g in hits_by_gid]
        return (sum(1 for g in G if not exclude_hits(hits_by_gid[g], g)) / len(G)) if G else None, len(G)

    fam, n_fam = nohit_share(lambda r: r == "Family")
    # "nothing at order level or below": no other panel organism shares the organism's
    # (non-empty) Order -- the per-rank reading of RESULTS_v2.md, which reproduces 72%.
    # The closest-relative bucketing (closest panel relative above Order) is kept as a variant.
    G = [g for g in truth_lin if g in hits_by_gid and not (truth_lin[g].get("Order")
                                                          and cnt["Order"][truth_lin[g]["Order"]] >= 2)]
    beyond = (sum(1 for g in G if not exclude_hits(hits_by_gid[g], g)) / len(G)) if G else None
    n_beyond = len(G)
    beyond_cl, n_beyond_cl = nohit_share(lambda r: r not in ("Genus", "Family", "Order"))
    cls, n_cls = nohit_share(lambda r: r == "Class")
    return dict(n_congener=len(cong), n_congener_with_mag=len(cong_mag),
                median_best_af=_median(best_af),
                below_gate=sum(1 for a in best_af if a < 15.0) / len(best_af) if best_af else None,
                nohit=nohit / len(cong_mag) if cong_mag else None,
                nohit_family_only=fam, n_family_only=n_fam,
                nohit_beyond_order=beyond, n_beyond_order=n_beyond,
                nohit_closest_above_order=beyond_cl, n_closest_above_order=n_beyond_cl,
                nohit_class_only=cls, n_class_only=n_cls,
                closest_rank_counts=dict(collections.Counter(str(v) for v in cl.values())))


def wald_diff_ci(p1, n1, p2, n2, z=1.959964):
    se = math.sqrt(p1 * (1 - p1) / n1 + p2 * (1 - p2) / n2)
    d = p1 - p2
    return d, d - z * se, d + z * se


def modal_fraction(counter):
    """Mode of the fractional (< 1) weights; ties broken by the larger count then value."""
    if not counter:
        return None, 0
    v, c = max(counter.items(), key=lambda kv: (kv[1], kv[0]))
    return v, c


# ============================================================================ loaders
def load_amp_cells(path, regions=None):
    """gene_capture_cells.csv, exclude_0.987 -> (gid, region) -> (outcome, f1, bnsi, n_removed_self)."""
    out = {}
    with open(path) as fh:
        for r in csv.DictReader(fh):
            if r["self_mode"] != "exclude_0.987" or (regions and r["region"] not in regions):
                continue
            out[(r["src_genome_id"], r["region"])] = (r["outcome"], _num(r["f1"]),
                                                     _num(r["best_nonself_identity"]),
                                                     int(r["n_removed_self"] or 0))
    return out


def load_mag_by_comp(paths, regime="exclude"):
    """mag cells -> completeness -> gid -> (outcome, has_precision, f1); contamination 0."""
    out = collections.defaultdict(dict)
    for p in paths:
        with open(p) as fh:
            for r in csv.DictReader(fh):
                if r["regime"] != regime or float(r["contamination"]) != 0.0:
                    continue
                out[float(r["completeness"])][r["genome_id"]] = (r["outcome"], r["precision"] != "",
                                                                 _num(r["f1"]))
    return out


def load_truth(path):
    out = {}
    for v in json.load(open(path)).values():
        out.setdefault(v["src_genome_id"], v)
    return out


# ============================================================================ recipe groups
def group_sharedref(res):
    rel_cells = "region_validation/shared_reference/data/cells.csv.gz"
    rel_pe = "region_validation/shared_reference/data/panel_eligible.json"
    rel_pred = ["region_validation/shared_reference/data/prfba/f%02d_d100/predictions.json" % f
                for f in range(5)]
    p_cells = res.path(rel_cells)
    recs = {"prfba_exact_tie": [], "prfba_legacy": []}
    with gzip.open(p_cells, "rt") as fh:
        for r in csv.DictReader(fh):
            if r["density"] == "1.0" and r["space"] == "ec" and r["arm"] in recs:
                recs[r["arm"]].append((r["md5"], r["fold"], r["n_reps"]))
    et = selection_size_stats(recs["prfba_exact_tie"])
    lg = selection_size_stats(recs["prfba_legacy"])
    # legacy single share over the exact-tie answered (md5, fold) set ("the same amplicons")
    lg_same = sum(1 for k in et["by_key"] if lg["by_key"].get(k) == 1) / et["answered"]
    I = [rel_cells]
    rec = ("cells.csv.gz, arm prfba_exact_tie, density 1.0, space ec, n_reps > 0, "
           "de-duplicated on (md5, fold)")
    out = [
        number("sr_fold0_single", et["fold0_single"], "int", "9,846", "219",
               "9{,}846 of 11{,}751 answered amplicons in fold~0", rec + "; fold 0", I),
        number("sr_fold0_answered", et["fold0_answered"], "int", "11,751", "219",
               "9{,}846 of 11{,}751 answered amplicons in fold~0", rec + "; fold 0", I),
        number("sr_fold0_single_share", et["fold0_single"] / et["fold0_answered"], "pct1", "83.8",
               "219", "in fold~0 (83.8\\%)", rec + "; fold 0", I),
        number("sr_fold0_two", et["fold0_two"], "int", "1,354", "219", "two for 1{,}354 (11.5\\%)",
               rec + "; fold 0, n_reps == 2", I),
        number("sr_fold0_two_share", et["fold0_two"] / et["fold0_answered"], "pct1", "11.5", "219",
               "two for 1{,}354 (11.5\\%)", rec + "; fold 0", I),
        number("sr_fold0_three_plus", et["fold0_three_plus"], "int", "551", "219",
               "three or more for 551 (4.7\\%)", rec + "; fold 0, n_reps >= 3", I),
        number("sr_fold0_three_plus_share", et["fold0_three_plus"] / et["fold0_answered"], "pct1",
               "4.7", "219", "three or more for 551 (4.7\\%)", rec + "; fold 0", I),
        number("sr_single_share", et["single_share"], "pct1", "86.8", "219, 303, 352, 1170, 1229",
               "86.8\\% of 50{,}842 answered amplicons", rec + "; pooled over five folds", I),
        number("sr_answered", et["answered"], "int", "50,842", "219, 352, 1170",
               "86.8\\% of 50{,}842 answered amplicons", rec, I),
        number("sr_single_n", et["single"], "int", "44,132", "1170", "(44{,}132 of 50{,}842",
               rec + "; n_reps == 1", I),
        number("sr_mean_size", et["mean"], "f3", "1.199", "219, 1170", "mean selection size 1.199",
               rec + "; mean n_reps", I),
        number("sr_legacy_single_share", lg_same, "pct1", "61.7", "352, 1170, 1229",
               "61.7\\% of the same amplicons under the legacy reducer",
               "cells.csv.gz arm prfba_legacy n_reps == 1 over the exact-tie answered (md5, fold) set",
               I, extra=dict(legacy_answered=lg["answered"], legacy_single_share_own=lg["single_share"])),
    ]
    # tie-set coherence
    p_pe = res.path(rel_pe)
    lineage = json.load(open(p_pe))["genomes"]
    import ijson
    multi = []
    for f, rel in enumerate(rel_pred):
        p = res.path(rel)
        with open(p, "rb") as fh:
            for md5, r in ijson.kvitems(fh, "pred.exact_tie"):
                if et["by_key"].get((md5, str(f)), 0) >= 2:
                    multi.append(tuple(r["reps"]))
    n_multi_expected = sum(1 for v in et["by_key"].values() if v >= 2)
    if len(multi) != n_multi_expected:
        raise SystemExit("[derived] tie sets: %d multi-genome selections in predictions vs %d in cells"
                         % (len(multi), n_multi_expected))
    tc = tie_set_coherence(multi, lineage)
    I2 = I + [rel_pe] + rel_pred
    rec2 = ("exact_tie reps of every answered (md5, fold) selection with n_reps >= 2 "
            "(prfba/f0X_d100/predictions.json), lineages from panel_eligible.json")
    out += [
        number("sr_tieset_span_phyla_unique", tc["span_unique"], "pct1", "6.4", "232, 1171, 1236",
               "6.4\\% unite genomes from different phyla", rec2 + "; unique tie sets", I2,
               extra=dict(n_unique=tc["n_unique"])),
        number("deck_XPHYLUM_BENCH", tc["span_unique"], "pct1", "6.4", None, None,
               "= sr_tieset_span_phyla_unique", I2, deck="XPHYLUM_BENCH"),
        number("sr_tieset_share_genus_unique", tc["genus_unique"], "pct1", "35.8", "232",
               "only 35.8\\% share a genus", rec2 + "; unique tie sets, one non-empty genus", I2),
        number("sr_tieset_span_phyla_records", tc["span_records"], "pct1", "4.8", "1236",
               "4.8\\% of multi-genome selections (6.4\\% of unique tie sets)",
               rec2 + "; per selection", I2, extra=dict(n_records=tc["n_records"])),
    ]
    return out


def _scan_synthetic_dir(d):
    """Fractional feature weights of multi-genome syntheses in one build directory."""
    try:
        import orjson
        loads = orjson.loads
    except ImportError:
        loads = json.loads
    frac = collections.Counter()
    n_gen = n_multi = n_active = 0
    lines = []
    for f in sorted(glob.glob(os.path.join(d, "*.json"))):
        b = open(f, "rb").read()
        lines.append("%s %s" % (os.path.basename(f), hashlib.sha256(b).hexdigest()))
        g = loads(b)
        if not isinstance(g, dict) or "features" not in g:
            continue
        n_gen += 1
        if "|" in str(g.get("source_id", "")):
            n_multi += 1
        fr = [round(float(ft["probability"]), 6) for ft in g["features"]
              if ft.get("probability") is not None and float(ft["probability"]) < 1.0]
        if fr:
            n_active += 1
        frac.update(fr)
    digest = hashlib.sha256("\n".join(lines).encode()).hexdigest()
    return frac, n_gen, n_multi, n_active, digest


def group_recon(res):
    if res.state == "pre":
        d, label = EMILYKIN_SHIPPED_D100, "EmilyKin shipped d100 k/n build (mixed annotation)"
    else:
        d, label = os.path.join(res.repo, KN_EMILYKIN_BVBRC), "EmilyKin BV-BRC-only d100 k/n arm"
    frac, n_gen, n_multi, n_active, digest = _scan_synthetic_dir(d)
    res.external_dir("synthetic_genomes:" + label, d, digest)
    mode, cnt = modal_fraction(frac)
    tot = sum(frac.values())
    return [number("recon_modal_fractional_weight", mode, "f1", "0.5", "225",
                   "the modal fractional weight is 0.5",
                   "mode of the k/n feature probabilities < 1 over the deployment's synthetic genomes (%s)"
                   % label, ["synthetic_genomes:" + label],
                   extra=dict(n_genomes=n_gen, n_multi_source=n_multi, n_with_fractional=n_active,
                              n_fractional_features=tot, mode_count=cnt,
                              mode_share=cnt / tot if tot else None,
                              top5=[[k, v] for k, v in frac.most_common(5)]))]


def group_as10(res):
    rel = "region_validation_d100/picrust2_comparison/data/comparison_cells.csv.gz"
    p = res.path(rel)
    pr, pc = {}, {}
    with gzip.open(p, "rt") as fh:
        for r in csv.DictReader(fh):
            if r["space"] != "ec":
                continue
            if r["method"] == "prfba_exclude":
                pr[r["amplicon"]] = (r["outcome"], _num(r["f1"]))
            elif r["method"] == "picrust2":
                pc[r["amplicon"]] = (r["outcome"], _num(r["f1"]))
    both = [a for a in pr if pr[a][0] == "ok" and pc.get(a, ("",))[0] == "ok"]
    m1, m2 = _mean(pr[a][1] for a in both), _mean(pc[a][1] for a in both)
    rec = "comparison_cells.csv.gz space ec, amplicons with outcome ok for both prfba_exclude and picrust2"
    return [
        number("as10_n_matched", len(both), "int", "59,607", "521", "on the 59{,}607 amplicons both methods answer",
               rec, [rel]),
        number("as10_prfba_f1", m1, "f3", "0.802", "74, 521, 1187, 1216",
               "conditional EC F1 is 0.802 versus 0.727", rec + "; mean F1 prFBA", [rel]),
        number("as10_picrust2_f1", m2, "f3", "0.727", "74, 521, 1187, 1216",
               "conditional EC F1 is 0.802 versus 0.727", rec + "; mean F1 PICRUSt2", [rel]),
        number("as10_delta", m1 - m2, "sf3", "+0.075", "1187", "$+0.075$ on matched amplicons (0.802 against 0.727)",
               rec + "; difference of means", [rel]),
    ]


# tab:magvs16s as printed: n, MAG, 16S, delta, median delta, cov MAG, cov 16S, pop MAG, pop 16S
TAB_MAGVS16S = {
    "FullLength16S": ("1,916", "0.813", "0.710", "+0.103", "+0.056", "0.304", "0.921", "0.248", "0.605"),
    "V1-V2": ("1,992", "0.814", "0.698", "+0.116", "+0.059", "0.298", "0.920", "0.242", "0.593"),
    "V1-V3": ("2,020", "0.813", "0.704", "+0.110", "+0.052", "0.298", "0.927", "0.243", "0.602"),
    "V3-V4": ("1,885", "0.815", "0.629", "+0.186", "+0.148", "0.298", "0.875", "0.242", "0.516"),
    "V4": ("2,299", "0.810", "0.598", "+0.212", "+0.171", "0.313", "0.863", "0.253", "0.483"),
    "V4-V5": ("2,334", "0.809", "0.616", "+0.192", "+0.157", "0.314", "0.868", "0.254", "0.504"),
    "V4-V5_944R": ("2,255", "0.811", "0.623", "+0.188", "+0.157", "0.318", "0.873", "0.257", "0.510"),
    "V7-V9": ("2,301", "0.808", "0.641", "+0.168", "+0.126", "0.326", "0.880", "0.264", "0.523"),
}


def _tab_row_snippet(region, pr):
    n = pr[0].replace(",", "{,}")
    return "%s & %s & %s & %s & $%s$ ($%s$) & %s & %s & %s & \\textbf{%s}" % (
        TEX_REGION[region], n, pr[1], pr[2], pr[3], pr[4], pr[5], pr[6], pr[7], pr[8])


def _amp_d100(res, cache):
    rel = "region_validation_d100/data/gene_capture_cells.csv"
    p = res.path(rel)
    if "amp_d100" not in cache:
        cache["amp_d100"] = load_amp_cells(p)
    return cache["amp_d100"]


def group_magvs16s(res, cache):
    rel_amp = "region_validation_d100/data/gene_capture_cells.csv"
    rel_mag = "mag_validation/out/mag_cells_all.csv"
    rel_tau = "mag_validation/out/mag_cells_tau99.csv"
    rel_tau_cmp = "region_validation_d100/data/gene_capture_tau_compare.json"
    amp = _amp_d100(res, cache)
    lowc = res.glob("mag_validation/lowc_out_batches/b*/mag_cells.csv")
    mag = load_mag_by_comp([res.path(rel_mag)] + lowc)
    cache["mag_c100"] = mag[1.0]
    I = [rel_amp, rel_mag, "mag_validation/lowc_out_batches/b*/mag_cells.csv"]
    rec = ("compare_arms.py pairing: MAG cells exclude, contamination 0, ok = outcome ok with "
           "precision; amplicon gene_capture_cells.csv exclude_0.987, ok = outcome ok")
    out = []
    pairs = {}
    for comp in sorted(mag, reverse=True):
        for reg in REGIONS:
            pairs[(comp, reg)] = pair_mag_amp(mag[comp], amp, reg)
    # tab:magvs16s at 100% completeness
    for reg in REGIONS:
        s = pairs[(1.0, reg)]
        pr = TAB_MAGVS16S[reg]
        snip = _tab_row_snippet(reg, pr)
        cols = (("n", s["n"], "int"), ("mag", s["mag"], "f3"), ("amp", s["amp"], "f3"),
                ("delta", s["delta"], "sf3"), ("median_delta", s["median_delta"], "sf3"),
                ("cov_mag", s["cov_mag"], "f3"), ("cov_amp", s["cov_amp"], "f3"),
                ("pop_mag", s["pop_mag"], "f3"), ("pop_amp", s["pop_amp"], "f3"))
        for i, (col, v, kind) in enumerate(cols):
            out.append(number("magvs16s_%s_%s" % (reg, col), v, kind, pr[i], "842-865",
                              snip,
                              rec + "; completeness 1.0; column %s" % col, I))
    ds = {k: v["delta"] for k, v in pairs.items()}
    at1 = [v for (c, r), v in ds.items() if c == 1.0]
    at01 = [v for (c, r), v in ds.items() if c == 0.01]
    allv = list(ds.values())
    levels = sorted({c for c, _ in ds})
    out += [
        number("magvs16s_delta_min_c100", min(at1), "sf2", "+0.10", "63", "$+0.10$ to $+0.21$ F1 at full completeness",
               rec + "; min over 8 regions of the conditional paired delta at 100%", I),
        number("magvs16s_delta_max_c100", max(at1), "sf2", "+0.21", "63", "$+0.10$ to $+0.21$ F1 at full completeness",
               rec + "; max over 8 regions at 100%", I),
        number("magvs16s_delta_min_c1", min(at01), "sf2", "+0.08", "63", "narrowing to $+0.08$ to $+0.19$ at 1\\% completeness",
               rec + "; min over 8 regions at 1%", I),
        number("magvs16s_delta_max_c1", max(at01), "sf2", "+0.19", "63", "narrowing to $+0.08$ to $+0.19$ at 1\\% completeness",
               rec + "; max over 8 regions at 1%", I),
        number("magvs16s_positive_everywhere", all(v is not None and v > 0 for v in allv), "bool", "yes", "63",
               "positive in every region at every level tested",
               rec + "; every region x every completeness level (%d levels)" % len(levels), I,
               extra=dict(n_cells=len(allv), levels=levels, min=min(allv), max=max(allv),
                          argmin=min(ds, key=ds.get), argmax=max(ds, key=ds.get))),
        number("magvs16s_delta_min_all_levels", min(allv), "f2", "0.08", "882", "by 0.08--0.21 F1, at any completeness at or above 1\\%",
               rec + "; min over all levels >= 1% and regions", I),
        number("magvs16s_delta_max_all_levels", max(allv), "f2", "0.21", "882", "by 0.08--0.21 F1, at any completeness at or above 1\\%",
               rec + "; max over all levels >= 1% and regions", I),
    ]
    # scope claims (L882-885) at 100%
    m100 = mag[1.0]
    scoreable = {g: v for g, v in m100.items() if v[0] != "truth_empty"}
    resolved = sum(1 for v in scoreable.values() if v[0] == "ok" and v[1])
    rr = resolved / len(scoreable)
    notans = {}
    for reg in ("FullLength16S", "V4"):
        vals = [amp[(g, reg)][1] for g, v in scoreable.items()
                if not (v[0] == "ok" and v[1]) and (g, reg) in amp and amp[(g, reg)][0] == "ok"]
        notans[reg] = (_mean(vals), len(vals))
    fl = pairs[(1.0, "FullLength16S")]
    all_fl_ok = [v[1] for (g, r), v in amp.items() if r == "FullLength16S" and v[0] == "ok"]
    easier = fl["amp"] - _mean(all_fl_ok)
    out += [
        number("mag_resolve_rate_pct", rr, "pct0", "31", "882, 884", "covers the 31\\% of organisms the MAG arm answers",
               "mag_cells_all.csv c100_x0 exclude: ok with precision / scoreable (outcome != truth_empty)", [rel_mag],
               extra=dict(resolved=resolved, scoreable=len(scoreable))),
        number("mag_abstain_rate_pct", 1 - rr, "pct0", "69", "870, 884", "on the remaining 69\\% it abstains",
               "1 - mag_resolve_rate", [rel_mag]),
        number("mag_resolve_rate", rr, "f3", "0.308", "878, 890", "the resolve rate itself, after selection, is 0.308",
               "as mag_resolve_rate_pct", [rel_mag]),
        number("magvs16s_amp_f1_where_mag_abstains_FL", notans["FullLength16S"][0], "f2", "0.64", "884",
               "answers at F1 0.64 (full-length) and 0.54 (V4)",
               "mean amplicon ok-cell F1 (FL) over scoreable organisms the MAG arm (100%) does not answer", I,
               extra=dict(n=notans["FullLength16S"][1])),
        number("magvs16s_amp_f1_where_mag_abstains_V4", notans["V4"][0], "f2", "0.54", "884",
               "answers at F1 0.64 (full-length) and 0.54 (V4)", "as above, V4", I,
               extra=dict(n=notans["V4"][1])),
        number("magvs16s_pop_amp_FL", fl["pop_amp"], "f3", "0.605", "884", "the amplicon arm wins, 0.605 against 0.248",
               rec + "; population FL", I),
        number("magvs16s_pop_mag_FL", fl["pop_mag"], "f3", "0.248", "884", "the amplicon arm wins, 0.605 against 0.248",
               rec + "; population FL", I),
        number("magvs16s_mag_answered_easier", easier, "f3", "0.053", "885", "are also 0.053 easier for the amplicon arm",
               "FL amplicon ok-cell mean F1 on MAG-answered organisms (paired set) minus over all FL ok cells", I,
               extra=dict(paired=fl["amp"], all_fl_ok=_mean(all_fl_ok), n_all=len(all_fl_ok))),
    ]
    # tau = 99 (L892-895)
    magt = load_mag_by_comp([res.path(rel_tau)])[1.0]
    t = pair_mag_amp(magt, amp, "FullLength16S")
    tc = json.load(open(res.path(rel_tau_cmp)))
    f0995 = tc["regions"]["FullLength16S"]["exclude_0.995"]["f1_population"]["value"]
    It = [rel_amp, rel_tau, rel_tau_cmp]
    out += [
        number("magvs16s_tau99_delta_FL", t["delta"], "sf3", "+0.127", "892", "widens from $+0.103$ to $+0.127$",
               "mag_cells_tau99.csv (exclude, c1, x0) paired with d100 FL cells; conditional delta", It,
               extra=dict(n=t["n"], mag=t["mag"], amp=t["amp"])),
        number("magvs16s_tau99_pop_mag_FL", t["pop_mag"], "f3", "0.301", "892, 894", "population F1 rises from 0.248 to 0.301",
               "as above; population MAG", It),
        number("magvs16s_tau99_pop_amp_FL", t["pop_amp"], "f3", "0.605", "892", "against the amplicon's unchanged 0.605",
               "as above; population 16S", It),
        number("magvs16s_tau0995_FL_pop", f0995, "f3", "0.655", "894", "(0.301 against 0.655)",
               "gene_capture_tau_compare.json regions.FullLength16S.exclude_0.995.f1_population.value", It),
        number("magvs16s_gap_matched", f0995 - t["pop_mag"], "f3", "0.354", "894", "(0.354 against 0.357",
               "0.995 FL population F1 minus tau-99 MAG population F1 (rounded inputs in the text)", It),
        number("magvs16s_gap_conservative", fl["pop_amp"] - fl["pop_mag"], "f3", "0.357", "894",
               "(0.354 against 0.357, from 0.248 against 0.605)", "FL population 16S minus MAG at tau 95", I),
    ]
    return out


def _load_c100_hits(res):
    meta = {}
    for p in res.glob("mag_validation/manifests/b*.json"):
        for r in json.load(open(p))["rows"]:
            if r["completeness_target"] == 1.0 and r["contamination_target"] == 0:
                meta[r["mag"]] = r
    hits = collections.defaultdict(list)
    for p in res.glob("mag_validation/hits/b*.tsv"):
        with open(p) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                q = os.path.basename(row["Query_file"])[:-4]
                if q in meta:
                    hits[q].append(dict(ref=os.path.basename(row["Ref_file"])[:-4],
                                        ani=float(row["ANI"]),
                                        af_query=float(row["Align_fraction_query"])))
    return meta, hits


def group_mag_hits(res, cache):
    meta, hits = _load_c100_hits(res)
    gid_by_mag = {m: r["genome_id"] for m, r in meta.items()}
    rel_mag = "mag_validation/out/mag_cells_all.csv"
    scoreable = []
    with open(res.path(rel_mag)) as fh:
        for r in csv.DictReader(fh):
            if (r["regime"] == "exclude" and float(r["completeness"]) == 1.0
                    and float(r["contamination"]) == 0 and r["outcome"] != "truth_empty"):
                scoreable.append(r["mag"])
    missing = [m for m in scoreable if m not in meta]
    if missing:
        raise SystemExit("[derived] %d scoreable MAGs absent from the manifests" % len(missing))
    sw = af_sweep(scoreable, hits, gid_by_mag)
    I = ["mag_validation/manifests/b*.json", "mag_validation/hits/b*.tsv", rel_mag]
    rec = ("share of scoreable c100_x0 organisms (mag_cells_all.csv outcome != truth_empty) with a "
           "surviving hit (score_mag.exclude_hits, tau 95) at ANI >= floor and AF_query >= gate")
    deck_printed = ("0.311", "0.316", "0.379", "0.452", "0.567", "0.655")
    out = []
    for i, (s, pv) in enumerate(zip(sw, deck_printed)):
        out.append(number("deck_AF_SWEEP_%d" % i, s["share"], "f3", pv, None, None,
                          rec + "; " + s["label"], I, deck="AF_SWEEP[%d]" % i,
                          extra=dict(k=s["k"], n=s["n"], floor=s["floor"], gate=s["gate"])))
    out += [
        number("mag_af_shipped", sw[0]["share"], "f3", "0.311", "878", "from 0.311 to 0.316", rec + "; as shipped", I),
        number("mag_af_floor0", sw[1]["share"], "f3", "0.316", "878", "from 0.311 to 0.316", rec + "; ANI floor 0", I),
        number("mag_af_gate2", sw[3]["share"], "f3", "0.452", "878", "moves it to 0.452", rec + "; AF gate 2%", I),
    ]
    # marker reach
    rel_truth = "region_validation/data/benchmark_truth.json"
    truth = {g: v["truth_lineage"] for g, v in load_truth(res.path(rel_truth)).items()}
    hits_by_gid = {gid_by_mag[m]: hits.get(m, []) for m in meta}
    mr = marker_reach(truth, hits_by_gid)
    Ir = I[:2] + [rel_truth]
    rec_r = ("panel organisms (benchmark_truth.json) sharing a non-empty Genus with another panel "
             "organism; best surviving hit = highest AF_query after score_mag exclusion")
    out += [
        number("mag_reach_n_congener", mr["n_congener"], "int", "6,011", "871", "among the 6{,}011 panel organisms",
               rec_r, Ir, extra=dict(n_with_mag=mr["n_congener_with_mag"])),
        number("deck_MAG_REACH_0", mr["n_congener"], "str", "6011", None, None, "= mag_reach_n_congener", Ir,
               deck="MAG_REACH[0]"),
        number("mag_reach_median_af", mr["median_best_af"], "f1", "20.3", "872",
               "aligns over a median of just 20.3\\% of the query", rec_r + "; median over those with a surviving hit", Ir),
        number("deck_MAG_REACH_1", mr["median_best_af"], "f1", "20.3", None, None, "= mag_reach_median_af", Ir,
               deck="MAG_REACH[1]"),
        number("mag_reach_below_gate", mr["below_gate"], "pct1", "45.1", "872", "45.1\\% fall below the 15\\% AF gate",
               rec_r + "; share with best AF_query < 15 among those with a surviving hit", Ir),
        number("deck_MAG_REACH_2", mr["below_gate"], "pct1", "45.1", None, None, "= mag_reach_below_gate", Ir,
               deck="MAG_REACH[2]"),
        number("mag_reach_nohit_congener", mr["nohit"], "pct1", "15.9", "875",
               "leaves 15.9\\% of congener-bearing organisms with no surviving hit",
               rec_r + "; share with no surviving hit", Ir),
        number("mag_reach_nohit_beyond_order", mr["nohit_beyond_order"], "pct0", "72", "875",
               "rising to 72\\% for organisms with no relative at order level or below",
               "share with no surviving hit among organisms with no other panel organism in their "
               "(non-empty) Order; family-only relatives give 51% (RESULTS_v2.md)", Ir,
               extra=dict(n=mr["n_beyond_order"], family_only=mr["nohit_family_only"],
                          closest_relative_above_order=mr["nohit_closest_above_order"],
                          n_closest_relative_above_order=mr["n_closest_above_order"],
                          n_family_only=mr["n_family_only"], class_only=mr["nohit_class_only"],
                          n_class_only=mr["n_class_only"], closest_rank_counts=mr["closest_rank_counts"])),
    ]
    return out


def group_identity_bins(res, cache):
    rel100 = "region_validation_d100/data/gene_capture_cells.csv"
    rel20 = "region_validation/data/gene_capture_cells.csv"
    rel_mag = "mag_validation/out/mag_cells_all.csv"
    amp100 = _amp_d100(res, cache)
    amp20 = load_amp_cells(res.path(rel20), regions=("FullLength16S",))
    p_mag = res.path(rel_mag)
    mag = cache.get("mag_c100") or load_mag_by_comp([p_mag])[1.0]
    pop = {g: (v[0] == "ok" and v[1]) for g, v in mag.items() if v[0] in ("ok", "abstained")}
    bins = (("0.97-0.987", 0.97, 0.987), ("0.945-0.97", 0.945, 0.97))

    def rates(amp):
        r = {}
        for lab, lo, hi in bins:
            G = [g for g in pop if (g, "FullLength16S") in amp and amp[(g, "FullLength16S")][2] is not None
                 and lo <= amp[(g, "FullLength16S")][2] < hi]
            A = [g for g in G if amp[(g, "FullLength16S")][0] != "truth_empty"]
            r[lab] = dict(n=len(G), mag=sum(pop[g] for g in G) / len(G),
                          amp=sum(amp[(g, "FullLength16S")][0] == "ok" for g in A) / len(A), n_amp=len(A))
        return r
    r100, r20 = rates(amp100), rates(amp20)
    rec = ("MAG population c100_x0 exclude (ok/abstained), binned by the amplicon FL exclude_0.987 "
           "best_nonself_identity; MAG answers = ok with precision; amplicon answers = outcome ok "
           "(pred_empty counts as no answer)")
    I100, I20 = [rel100, rel_mag], [rel20, rel_mag]
    lo_amp = min(v["amp"] for v in r100.values())
    hi_amp = max(v["amp"] for v in r100.values())
    return [
        number("mag_band_amp_min", lo_amp, "pct0", "98", "874", "the amplicon answers 98\\% of the time",
               rec + "; depth-100 bins, min over 0.945-0.97 / 0.97-0.987", I100, extra=r100),
        number("mag_band_amp_max", hi_amp, "pct0", "98", "874", "the amplicon answers 98\\% of the time",
               rec + "; depth-100 bins, max", I100),
        number("mag_band_mag_low", r100["0.945-0.97"]["mag"], "pct0", "8", "874", "whole-genome ANI answers 8--20\\%",
               rec + "; depth-100 bins, 0.945-0.97", I100),
        number("mag_band_mag_high", r100["0.97-0.987"]["mag"], "pct0", "20", "874", "whole-genome ANI answers 8--20\\%",
               rec + "; depth-100 bins, 0.97-0.987", I100),
        number("deck_BAND_MAG_0", r20["0.945-0.97"]["mag"], "f3", "0.082", None, None,
               rec + "; DEPTH-20 bins (region_validation, as RESULTS_v2.md)", I20, deck="BAND_MAG[0]",
               extra=dict(depth100=r100["0.945-0.97"]["mag"], depth20=r20)),
        number("deck_BAND_MAG_1", r20["0.97-0.987"]["mag"], "f3", "0.202", None, None,
               rec + "; DEPTH-20 bins", I20, deck="BAND_MAG[1]",
               extra=dict(depth100=r100["0.97-0.987"]["mag"])),
        # The deck constant was computed on the depth-20 bins (RESULTS_v2.md); the deck and the paper
        # are depth-100 throughout since c91cf19, so the depth-100 bins are carried as first-class
        # numbers for the deck update (the printed deck value is the d20 one; no d100 print exists).
        number("deck_BAND_MAG_d100_0", r100["0.945-0.97"]["mag"], "f3", None, None, None,
               rec + "; DEPTH-100 bins (region_validation_d100); no printed value: the deck prints the "
               "depth-20 deck_BAND_MAG_0", I100, deck="BAND_MAG[0] (depth-100 variant)"),
        number("deck_BAND_MAG_d100_1", r100["0.97-0.987"]["mag"], "f3", None, None, None,
               rec + "; DEPTH-100 bins; no printed value", I100, deck="BAND_MAG[1] (depth-100 variant)"),
    ]


def group_panels(res, cache):
    from scipy.stats import ttest_ind
    trees = {"original": "region_validation_d100", "kmeans": "region_validation_kmeans_d100",
             "random": "region_validation_random_d100"}
    cells = {}
    for p, t in trees.items():
        rel = t + "/data/gene_capture_cells.csv"
        cells[p] = _amp_d100(res, cache) if p == "original" else load_amp_cells(res.path(rel))
    It = [t + "/data/gene_capture_cells.csv" for t in trees.values()]
    out = []
    # Welch, original vs k-means, exclude_0.987 ok-cell F1 per region
    wel = {}
    for reg in REGIONS:
        a = [v[1] for (g, r), v in cells["original"].items() if r == reg and v[0] == "ok"]
        b = [v[1] for (g, r), v in cells["kmeans"].items() if r == reg and v[0] == "ok"]
        t = ttest_ind(a, b, equal_var=False)
        wel[reg] = dict(delta=_mean(a) - _mean(b), p=float(t.pvalue), n_original=len(a), n_kmeans=len(b))
    rec = ("scipy ttest_ind(equal_var=False) of exclude_0.987 ok-cell F1, region_validation_d100 vs "
           "region_validation_kmeans_d100 gene_capture_cells.csv")
    I2 = It[:2]
    below8 = sorted(r for r in REGIONS if wel[r]["p"] < 1e-8)
    out += [
        number("panels_welch_all_regions_p_lt_0.01", all(w["p"] < 0.01 for w in wel.values()), "bool", "yes", "1020",
               "significant in all eight regions (Welch $p<0.01$", rec, I2,
               extra=dict(max_p=max(w["p"] for w in wel.values()), per_region=wel)),
        number("panels_welch_regions_p_lt_1e-8", ",".join(below8), "str", "V4,V4-V5,V7-V9", "1020",
               "$p<10^{-8}$ for V4, V4--V5 and V7--V9", rec + "; regions with p < 1e-8", I2),
        number("panels_welch_FL_delta", wel["FullLength16S"]["delta"], "sf3", "+0.011", "1020",
               "full-length 16S included ($+0.011$, $p=2.3\\times10^{-4}$)", rec + "; FL mean difference", I2),
        number("panels_welch_FL_p", wel["FullLength16S"]["p"], "e1", "2.3e-04", "1020",
               "full-length 16S included ($+0.011$, $p=2.3\\times10^{-4}$)", rec + "; FL p", I2),
    ]
    # MAG arm: original n500 vs k-means, exclude, contamination 0, resolved cells
    import pandas as pd
    og_paths = res.glob("mag_validation/out_batches_n500/b*/mag_cells.csv")
    km_paths = res.glob("mag_validation_km/out_batches/b*/mag_cells.csv")

    def mcells(paths):
        d = pd.concat([pd.read_csv(f, dtype={"genome_id": str}) for f in paths], ignore_index=True)
        d = d[(d.contamination.astype(float) == 0) & d.outcome.isin(("ok", "abstained"))].copy()
        d["resolved"] = d.outcome.eq("ok") & d.precision.notna()
        return d
    og, km = mcells(og_paths), mcells(km_paths)
    Im = ["mag_validation/out_batches_n500/b*/mag_cells.csv", "mag_validation_km/out_batches/b*/mag_cells.csv"]
    mw = {}
    for comp in sorted(set(og.completeness) & set(km.completeness), reverse=True):
        a = og[(og.regime == "exclude") & (og.completeness == comp) & og.resolved].f1
        b = km[(km.regime == "exclude") & (km.completeness == comp) & km.resolved].f1
        t = ttest_ind(a, b, equal_var=False)
        mw[str(comp)] = dict(delta=float(a.mean() - b.mean()), p=float(t.pvalue), n=[int(len(a)), int(len(b))])
    out.append(number("panels_welch_mag_every_completeness_p_lt_1e-3", all(v["p"] < 1e-3 for v in mw.values()),
                      "bool", "yes", "1020", "and at every MAG completeness ($p<10^{-3}$)",
                      "Welch of resolved-MAG F1 (exclude, contamination 0), out_batches_n500 vs mag_validation_km",
                      Im, extra=dict(max_p=max(v["p"] for v in mw.values()), per_completeness=mw)))
    # singleton abstention and multiplicity (L1083-1085)
    rel_ot = "region_validation/data/benchmark_truth.json"
    rel_kt = "region_validation_kmeans/data/benchmark_truth.json"
    ab = {}
    for label, rel_t, d in (("original", rel_ot, og), ("kmeans", rel_kt, km)):
        tr = load_truth(res.path(rel_t))
        sp = {g: (v.get("truth_lineage") or {}).get("Species") for g, v in tr.items()}
        cnt = collections.Counter(sp.values())
        x = d[(d.regime == "exclude") & (d.completeness.astype(float) == 1.0)].copy()
        x["sp"] = x.genome_id.map(sp)
        x["mult"] = x.sp.map(lambda s: cnt[s] if isinstance(s, str) else None)
        x["abst"] = ~x.resolved
        spl = x[x.sp.notna()]
        ab[label] = dict(single=float(x[x.mult == 1].abst.mean()), n_single=int((x.mult == 1).sum()),
                         species_level=float(spl.abst.mean()), n_species_level=int(len(spl)),
                         ge20=float(spl[spl.mult >= 20].abst.mean()) if (spl.mult >= 20).any() else None,
                         n_ge20=int((spl.mult >= 20).sum()),
                         without_ge20=float(spl[spl.mult < 20].abst.mean()))
    k, o = ab["kmeans"], ab["original"]
    d, lo, hi = wald_diff_ci(k["single"], k["n_single"], o["single"], o["n_single"])
    recs = ("MAG cells exclude, completeness 1.0, contamination 0, outcome ok/abstained; abstain = not "
            "(ok with precision); species multiplicity from each panel's benchmark_truth.json")
    Is = Im + [rel_ot, rel_kt]
    out += [
        number("panels_mag_single_abst_kmeans", k["single"], "f4", "0.6895", "1083", "0.6895 on the $k$-means panel ($n=3{,}942$)",
               recs + "; species with one panel member", Is),
        number("panels_mag_single_n_kmeans", k["n_single"], "int", "3,942", "1083", "($n=3{,}942$)", recs, Is),
        number("panels_mag_single_abst_original", o["single"], "f4", "0.6883", "1083", "against 0.6883 on the original ($n=9{,}580$)",
               recs, Is),
        number("panels_mag_single_n_original", o["n_single"], "int", "9,580", "1083", "($n=9{,}580$)", recs, Is),
        number("panels_mag_single_diff", d, "sf3", "+0.001", "1083", "a difference of $+0.001$", recs + "; k-means minus original", Is),
        number("panels_mag_single_ci_lo", lo, "sf3", "-0.016", "1083", "(95\\% CI $-0.016$ to $+0.018$)",
               recs + "; Wald 95% CI of the difference of proportions", Is),
        number("panels_mag_single_ci_hi", hi, "sf3", "+0.018", "1083", "(95\\% CI $-0.016$ to $+0.018$)", recs, Is),
        number("panels_mag_ge20_abst_kmeans", k["ge20"], "f3", "0.559", "1084", "falls to 0.559 for organisms whose species has 20 or more members",
               recs + "; k-means species with >= 20 panel members", Is, extra=dict(n=k["n_ge20"])),
        number("panels_mag_species_level_abst_kmeans", k["species_level"], "f3", "0.658", "1085",
               "(0.658 against 0.688 over organisms with a species-level truth lineage)", recs, Is,
               extra=dict(n=k["n_species_level"])),
        number("panels_mag_species_level_abst_original", o["species_level"], "f3", "0.688", "1085",
               "(0.658 against 0.688 over organisms with a species-level truth lineage)", recs, Is,
               extra=dict(n=o["n_species_level"])),
        number("panels_mag_without_ge20_abst_kmeans", k["without_ge20"], "f3", "0.687", "1085",
               "without it the $k$-means panel abstains at 0.687", recs, Is),
    ]
    # emptied-record rates (L742): emptied = n_removed_self >= 100 (the record depth), the definition
    # of three_panel/emptied_strata.py, over outcome != truth_empty; FL and pooled over eight regions
    em = {}
    for p, c in cells.items():
        rows = [v for v in c.values() if v[0] != "truth_empty"]
        fl = [v for (g, r), v in c.items() if r == "FullLength16S" and v[0] != "truth_empty"]
        em[p] = dict(pooled=sum(v[3] >= 100 for v in rows) / len(rows), n_pooled=len(rows),
                     fl=sum(v[3] >= 100 for v in fl) / len(fl), n_fl=len(fl),
                     fl_emptied=sum(v[3] >= 100 for v in fl))
    rece = ("gene_capture_cells.csv exclude_0.987, outcome != truth_empty, emptied = n_removed_self >= 100 "
            "(three_panel/emptied_strata.py, depth 100)")
    for p, pr_fl, pr_pool in (("original", "4.7", "8.3"), ("kmeans", "22.9", "29.4"), ("random", None, None)):
        out.append(number("panels_emptied_pooled_%s" % p, em[p]["pooled"], "pct1", pr_pool,
                          "742" if pr_pool else None,
                          "(8.3\\% to 29.4\\% pooled over all eight regions)" if pr_pool else None,
                          rece + "; 8 regions pooled", It, extra=dict(n=em[p]["n_pooled"])))
        out.append(number("panels_emptied_FL_%s" % p, em[p]["fl"], "pct1", pr_fl,
                          "742" if pr_fl else None,
                          "emptied-record rate rises from 4.7\\% to 22.9\\%" if pr_fl else None,
                          rece + "; FullLength16S (equals emptied_strata_d100.json emptied_pct)", It,
                          extra=dict(n=em[p]["n_fl"], emptied=em[p]["fl_emptied"])))
    # surviving-uncovered rates (L1057-1060)
    sv = {}
    for p, c in cells.items():
        rows = [v for v in c.values() if v[0] != "truth_empty"]
        surv = [v for v in rows if v[3] < 100]
        fl = [v for (g, r), v in c.items() if r == "FullLength16S" and v[0] != "truth_empty" and v[3] < 100]
        strata = {}
        for lab, a, b in (("0-19", 0, 19), ("20-49", 20, 49), ("50-99", 50, 99)):
            s = [v for v in fl if a <= v[3] <= b]
            strata[lab] = sum(v[0] != "ok" for v in s) / len(s) if s else None
        sv[p] = dict(pooled=sum(v[0] != "ok" for v in surv) / len(surv), n_pooled=len(surv),
                     fl=sum(v[0] != "ok" for v in fl) / len(fl), n_fl=len(fl),
                     fl_nonok=sum(v[0] != "ok" for v in fl), strata=strata)
    strata_all = [x for p in sv.values() for x in p["strata"].values() if x is not None]
    recu = ("gene_capture_cells.csv exclude_0.987, outcome != truth_empty, surviving = n_removed_self < 100, "
            "uncovered = outcome != ok")
    pctp = {"original": ("3.4", "3.1"), "kmeans": ("2.8", "2.4"), "random": ("2.3", "1.8")}
    for p in ("original", "kmeans", "random"):
        out.append(number("panels_surv_uncovered_pooled_%s" % p, sv[p]["pooled"], "pct1", pctp[p][1], "1057",
                          "3.1\\%, 2.4\\% and 1.8\\% pooled over all eight regions", recu + "; 8 regions pooled",
                          It, extra=dict(n=sv[p]["n_pooled"])))
        out.append(number("panels_surv_uncovered_FL_%s" % p, sv[p]["fl"], "pct1", pctp[p][0], "1057",
                          "(3.4\\%, 2.8\\% and 2.3\\% of surviving full-length records", recu + "; FullLength16S",
                          It, extra=dict(n=sv[p]["n_fl"], nonok=sv[p]["fl_nonok"], strata=sv[p]["strata"])))
    out += [
        number("panels_surv_uncovered_FL_strata_min", min(strata_all), "pct1", "1.5", "1057",
               "1.5--3.7\\% in every full-length stratum", recu + "; FL strata 0-19/20-49/50-99, all panels", It),
        number("panels_surv_uncovered_FL_strata_max", max(strata_all), "pct1", "3.7", "1057",
               "1.5--3.7\\% in every full-length stratum", recu, It),
        number("panels_emptied_fold_min", 1.0 / max(sv[p]["fl"] for p in sv), "f0", "29", "1057",
               "some 29- to 43-fold below", "1 / largest FL surviving-uncovered rate (the emptied-record rate is 1.000)", It),
        number("panels_emptied_fold_max", 1.0 / min(sv[p]["fl"] for p in sv), "f0", "43", "1057",
               "some 29- to 43-fold below", "1 / smallest FL surviving-uncovered rate", It),
        number("panels_surv_uncovered_FL_spread_pts", 100 * (max(sv[p]["fl"] for p in sv) - min(sv[p]["fl"] for p in sv)),
               "f1", "1.1", "1058", "(1.1 points, and 1.3 pooled over all eight regions)",
               "max - min of the FL surviving-uncovered rate over the three panels, in points", It),
        number("panels_surv_uncovered_pooled_spread_pts",
               100 * (max(sv[p]["pooled"] for p in sv) - min(sv[p]["pooled"] for p in sv)),
               "f1", "1.3", "1058", "(1.1 points, and 1.3 pooled over all eight regions)",
               "as above, pooled over 8 regions", It),
    ]
    return out


PHENO_MAIN_ARMS = ("truth", "prfba_exact_tie", "picrust2", "base_nn1", "base_family_donor", "base_constant")
ATP_T = 1e-6


def _jsonl(path):
    with open(path) as fh:
        return [json.loads(l) for l in fh if l.strip()]


def _atp_informative_pairs(by, genomes, other, pipeline="controlled"):
    """phenotype_fold's paired set: truth, prfba_exact_tie and `other` all status ok and the truth ATP
    vector has at least one positive medium (the informative-truth-only denominator)."""
    out = []
    for g in genomes:
        t, p, q = (by.get((g, a, pipeline)) for a in ("truth", "prfba_exact_tie", other))
        if not (t and p and q) or any(x.get("status") != "ok" for x in (t, p, q)):
            continue
        tv = (t.get("atp") or {}).get("vector") or {}
        if not any(v is not None and v > ATP_T for v in tv.values()):
            continue
        out.append((g, p, q))
    return out


def group_phenotype(res):
    """Phenotype (A-C13 configuration) numbers no summary key stores: PH-01..PH-05, abs-phenotype-readout,
    disc.refusal_credit, disc.atp_vs_nn1, lim.phenotype_scope, and the amended-run counts of
    abs-phenotype-controls / AM-02 (C12_SCOPE.md). Ported from the re-score's
    tools/r4_pheno_manual_numbers.py, which reproduced every printed value on the pre-repair pass."""
    from scipy.stats import spearmanr, wilcoxon
    sr = "region_validation/shared_reference/data/"
    rel_atp, rel_atps = sr + "phenotype_atp_c13/per_genome.jsonl", sr + "phenotype_atp_c13/phenotype_summary.json"
    rel_gf, rel_gfs = sr + "phenotype_gf_c13/per_genome.jsonl", sr + "phenotype_gf_c13/phenotype_summary.json"
    rel_cells = sr + "cells.csv.gz"
    rows = _jsonl(res.path(rel_atp))
    S = json.load(open(res.path(rel_atps)))
    by = {(r["src_genome_id"], r["arm"], r["pipeline"]): r for r in rows}
    genomes = sorted({r["src_genome_id"] for r in rows})
    media = sorted({m for r in rows for m in ((r.get("atp") or {}).get("vector") or {})})
    sig = sorted({m for r in rows if r["arm"] in PHENO_MAIN_ARMS
                  for m, v in ((r.get("atp") or {}).get("vector") or {}).items() if v is not None and v > ATP_T})
    pairs = _atp_informative_pairs(by, genomes, "picrust2")
    f1p = [(p.get("ec_f1"), q.get("ec_f1"), (p.get("vs_truth") or {}).get("atp_hamming"),
            (q.get("vs_truth") or {}).get("atp_hamming")) for g, p, q in pairs]
    f1p = [x for x in f1p if x[0] is not None and x[1] is not None]
    rho, sp_p = spearmanr([a - b for a, b, _, _ in f1p], [c - d for _, _, c, d in f1p])
    c = S["contrasts"]["prfba_exact_tie__vs__picrust2"]["controlled"]
    mc = c.get("mcnemar_per_medium") or {}
    sigm = {m: v for m, v in mc.items() if isinstance(v, dict) and v.get("q_bh") is not None and v["q_bh"] <= 0.05}
    # McNemar discordant counts: b = prFBA matches truth on the medium and PICRUSt2 does not, c = the reverse
    # (the pre-registered pass's two media, printed "favor PICRUSt2", both have c > b)
    fav_p2 = sorted(m for m, v in sigm.items() if v.get("c", 0) > v.get("b", 0))
    fav_pr = sorted(m for m, v in sigm.items() if v.get("b", 0) > v.get("c", 0))
    I = [rel_atp, rel_atps]
    rec = ("phenotype_atp_c13/per_genome.jsonl, controlled pipeline; paired set = truth, prfba_exact_tie and "
           "picrust2 all status ok, truth ATP vector with >= 1 medium > 1e-6")
    out = [
        number("pheno_n_paired", len(pairs), "int", "755", "79, 632, 1193, 1222", "755",
               rec + "; equals the summary's contrasts.prfba_exact_tie__vs__picrust2.controlled.n_pairs", I,
               extra=dict(n_pairs_summary=c.get("n_pairs"))),
        number("pheno_n_paired_v4", sum(1 for g, p, q in pairs if p.get("region") == "V4"), "int", "747", "79",
               "747", rec + "; pairs whose amplicon is V4", I),
        number("pheno_ecf1_delta_vs_picrust2", sum(a - b for a, b, _, _ in f1p) / len(f1p), "sf3", "-0.028",
               "79, 1193", "only $-0.028$", rec + "; mean per-genome ec_f1, prfba_exact_tie minus picrust2", I,
               extra=dict(n=len(f1p))),
        number("pheno_media_in_vector", len(media), "int", "54", "79, 638, 1222", "22 of them", rec +
               "; media keys of the ATP vector (53 media + the empty-medium control)", I),
        number("pheno_media_with_signal", len(sig), "int", "22", "79, 638, 1222", "only 22 of them ever carry signal",
               "media with ATP > 1e-6 in any row of truth / prfba_exact_tie / picrust2 / base_nn1 / "
               "base_family_donor / base_constant", I, extra=dict(media=sig)),
        number("pheno_spearman_rho", rho, "sf2", "-0.14", "634, 1190", "$\\rho = -0.14$",
               rec + "; Spearman of (EC-F1 delta, ATP-Hamming delta), prFBA minus PICRUSt2", I,
               extra=dict(summary=c.get("ec_f1_delta_vs_phenotype_spearman"))),
        number("pheno_spearman_p", float(sp_p), "e0", "2e-04", "1190", "$p = 2\\times10^{-4}$",
               rec + "; p of that Spearman", I, extra=dict(p=float(sp_p))),
        number("pheno_mcnemar_n_q005", len(sigm), "int", "2", "79, 631, 639, 1185", "on two of the 22 testable media",
               "phenotype_summary.json contrasts.prfba_exact_tie__vs__picrust2.controlled.mcnemar_per_medium "
               "with q_bh <= 0.05", I, extra=dict(media={m: v for m, v in sigm.items()},
                                                    n_media_tested=len(mc))),
        number("pheno_mcnemar_n_favour_picrust2", len(fav_p2), "int", "2", "79, 631, 639, 1185",
               "on two of the 22 testable media", "as above; of those, media with c > b (PICRUSt2 matches truth "
               "where prFBA does not more often than the reverse)", I),
        number("pheno_mcnemar_n_favour_prfba", len(fav_pr), "int", "0", "631", None, "as above; b > c", I),
        number("pheno_mcnemar_media_favour_picrust2", ",".join(fav_p2), "str", "Ac.NO,Etho.O2", "631, 1185",
               "(Ac.NO and Etho.O2, $q = 0.002$)", "as above; the media, sorted", I),
        number("pheno_mcnemar_q_max", max((v["q_bh"] for v in sigm.values()), default=None), "f3", "0.002",
               "631, 1185", "$q = 0.002$", "as above; the largest q of those media (the text prints one q)", I),
    ]
    # EC F1 prFBA vs nn1 on the nn1 contrast's 755 genomes (PH-04, disc.atp_vs_nn1)
    pn = _atp_informative_pairs(by, genomes, "base_nn1")
    want = {(p["md5"], str(p["fold"])) for g, p, q in pn}
    f1 = {}
    p_cells = res.path(rel_cells)
    with gzip.open(p_cells, "rt") as fh:
        for r in csv.DictReader(fh):
            if (r["density"] == "1.0" and r["space"] == "ec" and r["arm"] in ("prfba_exact_tie", "base_nn1")
                    and (r["md5"], r["fold"]) in want):
                f1[(r["md5"], r["fold"], r["arm"])] = _num(r["f1"])
    a = [f1.get((p["md5"], str(p["fold"]), "prfba_exact_tie")) for g, p, q in pn]
    b = [f1.get((p["md5"], str(p["fold"]), "base_nn1")) for g, p, q in pn]
    ab = [(x, y) for x, y in zip(a, b) if x is not None and y is not None]
    w = wilcoxon([x for x, _ in ab], [y for _, y in ab])
    In = I + [rel_cells]
    recn = ("cells.csv.gz density 1.0, space ec, f1 of prfba_exact_tie and base_nn1 at each paired genome's "
            "(md5, fold), over the prfba_exact_tie__vs__base_nn1 paired set (truth, prFBA, nn1 ok; "
            "informative truth)")
    out += [
        number("pheno_nn1_n_paired", len(pn), "int", "755", "635-637, 1186", "($n = 755$",
               recn + "; n", In, extra=dict(n_with_both_f1=len(ab))),
        number("pheno_nn1_ecf1_prfba", _mean(x for x, _ in ab), "f3", "0.782", "1186", "EC F1 0.782 against 0.781",
               recn + "; mean prFBA", In),
        number("pheno_nn1_ecf1_nn1", _mean(y for _, y in ab), "f3", "0.781", "1186", "EC F1 0.782 against 0.781",
               recn + "; mean nn1", In),
        number("pheno_nn1_ecf1_p", float(w.pvalue), "f2", "0.64", "1186", "$p = 0.64$",
               recn + "; scipy wilcoxon (paired, zero_method wilcox)", In),
    ]
    # amended-run (A-C13 configuration) counts: C12_SCOPE.md, abs-phenotype-controls, AM-02
    S_gf = json.load(open(res.path(rel_gfs)))
    gf = _jsonl(res.path(rel_gf))
    gf_models = [r for r in gf if (r.get("gapfilled") or {}).get("growth") is not None]
    Ig = [rel_atp, rel_gf, rel_gfs]
    out += [
        number("pheno_atp_rows", len(rows), "int", "43,140", "80, 671, 679", "43{,}140",
               "rows of phenotype_atp_c13/per_genome.jsonl (genome x arm x pipeline)", Ig),
        number("pheno_atp_rows_empty_medium_atp", sum(1 for r in rows if (r.get("atp") or {}).get("empty_medium_atp")
                                                      not in (None, 0, 0.0)), "int", "0", "679",
               None, "ATP-pass rows whose empty-medium control makes ATP", Ig),
        number("pheno_gf_rows", len(gf), "int", "8,614", "80, 668, 679", "8{,}614",
               "rows of phenotype_gf_c13/per_genome.jsonl", Ig),
        number("pheno_gf_gated", sum((S_gf["meta"].get("n_gated_out") or {}).values()), "int", "0",
               "80, 668, 679", "0 of 8{,}614",
               "phenotype_gf_c13/phenotype_summary.json meta.n_gated_out summed over arms (C12 leg c)", Ig,
               extra=dict(n_gated_out=S_gf["meta"].get("n_gated_out"),
                          rows_sanity_gate_false=sum(1 for r in gf if (r.get("gapfilled") or {}).get("sanity_gate_ok")
                                                     is False),
                          rows_atp_gate_false=sum(1 for r in gf if (r.get("atp") or {}).get("gate_ok") is False))),
        number("pheno_gf_models", len(gf_models), "int", "2,362", None, None,
               "gapfilled models (rows with gapfilled.growth) in phenotype_gf_c13/per_genome.jsonl "
               "(C12_SCOPE.md: 'every one of the 2,362 gapfilled models')", Ig),
        number("pheno_gf_models_growing", sum(1 for r in gf_models if r["gapfilled"]["growth"] > ATP_T), "int",
               "2,362", None, None, "of those, growth > 1e-6 (C12_SCOPE.md)", Ig),
        number("pheno_gf_models_reversal", sum(1 for r in gf if (r.get("gapfilled") or {}).get("reversed_applied")),
               "int", "114", "673", "114 apply a reversal",
               "rows whose gapfilled.reversed_applied is non-empty (models applying a gapfiller reversal)", Ig),
        number("pheno_gf_n_genomes", S_gf["meta"].get("n_genomes"), "int", "200", "645, 1222", None,
               "phenotype_gf_c13 meta.n_genomes", Ig),
    ]
    return out


def group_mag_depth500(res):
    """The depth-500 MAG re-run comparison (L943-945, L1009-1013; ad hoc in paper f8047c5), ported from the
    re-score's tools/mag_depth500.py: depth-100 cells (mag_validation/out_batches, restricted to the MAGs
    of the depth-500 run) joined to out_batches_n500 on (mag, regime)."""
    def load(paths):
        d = {}
        for p in paths:
            with open(p) as fh:
                for r in csv.DictReader(fh):
                    d[(r["mag"], r["regime"])] = r
        return d
    A = load(res.glob("mag_validation/out_batches/b*/mag_cells.csv"))
    B = load(res.glob("mag_validation/out_batches_n500/b*/mag_cells.csv"))
    ks = [k for k in sorted(B) if k[1] == "exclude" and k in A]
    sc = [k for k in ks if B[k]["outcome"] != "truth_empty"]
    chg = collections.Counter()
    for k in ks:
        o1, o2 = A[k]["outcome"], B[k]["outcome"]
        if o1 != o2:
            chg["%s->%s%s" % (o1, o2, " (gene-less reference)" if o2 == "ok" and B[k]["precision"] in ("", "None")
                              else "")] += 1
    ninc = sum(1 for k in sorted(B) if k[1] == "include" and k in A and A[k]["outcome"] != B[k]["outcome"])
    sw = sum(1 for k in ks if A[k]["outcome"] == "ok" and B[k]["outcome"] == "ok"
             and (A[k]["anchor"], A[k]["n_selected"]) != (B[k]["anchor"], B[k]["n_selected"]))
    I = ["mag_validation/out_batches/b*/mag_cells.csv", "mag_validation/out_batches_n500/b*/mag_cells.csv"]
    rec = ("out_batches (depth 100) vs out_batches_n500 (depth 500) mag_cells.csv joined on (mag, regime), "
           "regime exclude; scoreable = depth-500 outcome != truth_empty")
    n_res = sum(v for t, v in chg.items() if t.endswith("ok") and "gene-less" not in t)
    return [
        number("mag_d500_n_scoreable", len(sc), "int", "49,185", "943-945, 1011", "49{,}185 scoreable cells",
               rec, I, extra=dict(n_cells=len(ks), n_missing_in_d100=sum(1 for k in B if k not in A))),
        number("mag_d500_n_outcome_changed", sum(chg.values()), "int", "88", "944, 1012", "88",
               rec + "; cells whose outcome differs", I, extra=dict(changes=dict(chg), n_include_changed=ninc)),
        number("mag_d500_n_newly_resolved", n_res, "int", "87", "944, 1012", "87",
               rec + "; changed to ok with a precision (newly resolved)", I),
        number("mag_d500_n_geneless", sum(v for t, v in chg.items() if "gene-less" in t), "int", "1", "944, 1012", "1",
               rec + "; changed to ok on a gene-less reference (empty precision)", I),
        number("mag_d500_n_selection_switches", sw, "int", "0", "944, 1012", None,
               rec + "; both resolved and (anchor, n_selected) differ", I),
    ]


def _sr_cells(path, arms):
    """SR cells table -> arm -> {(amplicon, fold): (outcome, f1, n_reps, novelty_bin, md5)} at density 1.0,
    space ec (the scorer's primary cell)."""
    out = collections.defaultdict(dict)
    with gzip.open(path, "rt") as fh:
        for r in csv.DictReader(fh):
            if r["density"] == "1.0" and r["space"] == "ec" and r["arm"] in arms:
                out[r["arm"]][(r["amplicon"], r["fold"])] = (r["outcome"], r["f1"], r["n_reps"],
                                                              r["novelty_bin"], r["md5"])
    return out


def _sr_pop(cells, arm, mask=None):
    """Population F1 as the scorer computes it: unscoreable cells dropped, f1 if scored else 0; with `mask`
    (key -> outcome of the random-abstain arm) a cell the mask does not score counts as 0."""
    v = []
    for k, (o, f, *_ ) in cells[arm].items():
        if o.startswith("unscoreable"):
            continue
        x = float(f) if o == "scored" else 0.0
        if mask is not None and mask.get(k) != "scored":
            x = 0.0
        v.append(x)
    return sum(v) / len(v) if v else None


def group_sr_masks(res):
    """Controlled-benchmark numbers recomputed from the SR cells tables (no summary contrast stores them):
    the random-abstention mask applied to each policy's own forced run (SR-09, abs-controlled-random-mask,
    FR-03; recipe UNION_THROUGHOUT_CHANGES.md 3.5 / 3.17: base_random_abstain is the exact-tie forced run
    masked to prFBA's coverage) and P(n_genomes >= 2) by novelty bin (sharedref_ngenomes_vs_identity)."""
    out = []
    arms = ("prfba_exact_tie", "prfba_exact_tie_forced", "prfba_legacy", "prfba_legacy_forced",
            "base_random_abstain")
    for rank, rel, pr in (("genus", "region_validation/shared_reference/data/cells.csv.gz",
                           dict(legacy_masked="0.668", lfl="+0.019", et_mask="+0.019")),
                          ("family", "region_validation/shared_reference/data_family/cells.csv.gz",
                           dict(legacy_masked=None, lfl="+0.016", et_mask="+0.015"))):
        c = _sr_cells(res.path(rel), arms)
        mask = {k: v[0] for k, v in c["base_random_abstain"].items()}
        pop = {a: _sr_pop(c, a) for a in arms}
        lm = _sr_pop(c, "prfba_legacy_forced", mask)
        rec = ("%s cells.csv.gz density 1.0, space ec, unscoreable dropped, population F1 = f1 if scored else 0; "
               "mask = base_random_abstain's scored cells" % rank)
        tl = "612" if rank == "genus" else "78, 626"
        out += [
            number("sr_%s_legacy_forced_masked" % rank, lm, "f3", pr["legacy_masked"],
                   "612" if rank == "genus" else None,
                   "the same mask scores 0.668" if rank == "genus" else None,
                   rec + "; prfba_legacy_forced with the random-abstain mask", [rel]),
            number("sr_%s_legacy_like_for_like" % rank, pop["prfba_legacy"] - lm, "sf3", pr["lfl"], tl,
                   "a like-for-like credit of $+0.019$" if rank == "genus" else "$+0.016$ like-for-like",
                   rec + "; prfba_legacy minus its masked forced run", [rel],
                   extra=dict(legacy=pop["prfba_legacy"], masked=lm)),
            number("sr_%s_exact_tie_vs_mask" % rank, pop["prfba_exact_tie"] - pop["base_random_abstain"], "sf3",
                   pr["et_mask"], "78, 612" if rank == "genus" else "78, 626",
                   None, rec + "; prfba_exact_tie minus base_random_abstain (its own masked forced run)", [rel],
                   extra=dict(exact_tie=pop["prfba_exact_tie"], random_abstain=pop["base_random_abstain"])),
            number("sr_%s_exact_tie_forced_gain" % rank, pop["prfba_exact_tie_forced"] - pop["prfba_exact_tie"], "sf3",
                   "+0.061" if rank == "genus" else None, "615" if rank == "genus" else None, None,
                   rec + "; prfba_exact_tie_forced minus prfba_exact_tie", [rel]),
        ]
        if rank == "genus":
            # P(n_genomes >= 2) by novelty bin over the exact-tie answered (md5, fold) selections
            sel = {}
            for (amp, fold), (o, f, n, nb, md5) in c["prfba_exact_tie"].items():
                if int(n or 0) > 0:
                    sel.setdefault((md5, fold), (int(n), nb))
            by = collections.defaultdict(lambda: [0, 0])
            for n, nb in sel.values():
                by[nb][0] += 1
                by[nb][1] += n >= 2
            bins = sorted(by)
            share = [by[b][1] / by[b][0] for b in bins]
            recb = ("SR genus cells.csv.gz prfba_exact_tie density 1.0 ec, n_reps > 0, one per (md5, fold); "
                    "share with n_reps >= 2 by novelty_bin")
            out += [
                number("sr_p_multi_monotone", all(a < b for a, b in zip(share, share[1:])), "bool", "yes", "228",
                       "rises monotonically with rank-1 identity", recb + "; strictly increasing over the bins", [rel],
                       extra=dict(bins={b: dict(n=by[b][0], share=s) for b, s in zip(bins, share)})),
                number("sr_p_multi_low_bin", share[0], "f3", "0.100", None, None, recb + "; lowest bin " + bins[0],
                       [rel]),
                number("sr_p_multi_high_bin", share[-1], "f3", "0.180", None, None, recb + "; highest bin " + bins[-1],
                       [rel]),
            ]
    return out


def group_as_selection(res):
    """As-shipped (depth-100) numbers that need prfba_selection.json or a per-genome join (AS-17, AS-18,
    AS-19, AS-23, fig:picrust2, tab:picrust2); recipes reproduce the printed values on the pre-repair files."""
    import ijson
    rel_sel = "region_validation_d100/picrust2_comparison/data/prfba_selection.json"
    rel_cc = "region_validation_d100/picrust2_comparison/data/comparison_cells.csv.gz"
    rel_ceil = "region_validation/picrust2_comparison/data/annotation_ceiling.json"
    S = {}
    with open(res.path(rel_sel), "rb") as fh:
        for k, v in ijson.kvitems(fh, ""):
            if k.startswith("_"):
                continue
            e = v.get("exclude_0.987") or {}
            S[k] = (int(e.get("n_selected") or 0),
                    float(e["best_identity"]) if e.get("best_identity") is not None else None,
                    float(v["best_nonself_identity"]) if v.get("best_nonself_identity") is not None else None)
    C, P2 = {}, []
    with gzip.open(res.path(rel_cc), "rt") as fh:
        for r in csv.DictReader(fh):
            if r["space"] != "ec":
                continue
            if r["method"] == "prfba_exclude":
                C[r["amplicon"]] = (r["outcome"], r["nsti_bin"], r["ident_bin"], int(r["n_reps"] or 0))
            elif r["method"] == "picrust2" and r["src_in_picrust2_ref"] in ("1", "True", "true"):
                P2.append((r["amplicon"], r["outcome"], r["precision"]))
    I = [rel_sel, rel_cc]
    hi = [a for a, v in C.items() if v[0] == "ok" and v[1] == ">=0.5"]
    abst = [a for a, v in C.items() if v[0] not in ("ok", "truth_empty") and v[3] == 0]
    floor = [a for a in abst if C[a][2] == "<0.865"]
    further = [a for a in abst if C[a][2] != "<0.865" and S.get(a, (0, None))[1] is not None and 0 < S[a][1] < 0.865]
    ans_n = [S[a][0] for a, v in C.items() if v[0] == "ok" and a in S]
    rec = "prfba_selection.json exclude_0.987 joined to comparison_cells.csv.gz prfba_exclude ec"
    out = [
        number("as17_n_high_nsti_answered", len(hi), "int", "603", "529-531", "all 603 answers", rec +
               "; outcome ok, nsti_bin >= 0.5", I),
        number("as17_n_nonself_ge_0995", sum(1 for a in hi if S[a][2] is not None and S[a][2] >= 0.995), "int", "248",
               "531", "for 248 of the 603", rec + "; best_nonself_identity >= 0.995", I),
        number("as17_all_below_0987", all(S[a][1] is not None and S[a][1] < 0.987 for a in hi), "bool", "yes", "531",
               "all 603 answers rest on references below 0.987 identity", rec + "; every exclude best_identity < 0.987", I),
        number("as17_median_identity", statistics.median(S[a][1] for a in hi), "f3", "0.956", "531", "(median 0.956)",
               rec + "; median exclude best_identity", I),
        number("as18_n_further_below_floor", len(further), "int", "145", "575",
               None, rec + "; outright abstentions (no reps) outside the <0.865 bin whose best surviving hit "
               "(exclude best_identity, 0 = emptied) is above 0 and below 0.865", I),
        number("as18_n_floor_explains", len(floor) + len(further), "int", "692", "575", None,
               rec + "; 547-bin abstentions plus the further ones", I),
        number("as18_floor_share_of_abstentions", (len(floor) + len(further)) / len(abst), "pct1", "10.8", "575",
               None, rec + "; share of outright abstentions (no reps)", I, extra=dict(n_abstentions=len(abst))),
        number("as23_median_selection_size", statistics.median(ans_n), "f0", "1", "540",
               "a median size of one reference", rec + "; median exclude n_selected over answered amplicons", I),
    ]
    ceil = {r["genome_id"]: r for r in json.load(open(res.path(rel_ceil)))["per_genome"]}
    n = k = 0
    for amp, o, p in P2:
        g = amp.split("__")[0]
        if o != "ok" or p == "" or g not in ceil:
            continue
        n += 1
        k += float(p) > ceil[g]["shared"]["precision"]
    out.append(number("as19_share_more_precise_than_reference", k / n if n else None, "pct0", "43", "548, 565",
                      None, "comparison_cells.csv.gz picrust2 ec, source in PICRUSt2's reference, outcome ok: "
                      "share whose precision exceeds annotation_ceiling.json per_genome[source].shared.precision",
                      [rel_cc, rel_ceil], extra=dict(n=n, k=k)))
    return out


GROUPS = (("sharedref", group_sharedref), ("recon", group_recon), ("as10", group_as10),
          ("magvs16s", group_magvs16s), ("mag_hits", group_mag_hits),
          ("identity_bins", group_identity_bins), ("panels", group_panels),
          ("phenotype", group_phenotype), ("mag_depth500", group_mag_depth500),
          ("sr_masks", group_sr_masks), ("as_selection", group_as_selection))

# Why a number does not reproduce the printed value on --state pre (filled from the control).
KNOWN_PRE_MISMATCH = {
    "pheno_nn1_ecf1_nn1": (
        "the printed 0.781 (paper f7fdc38, 2026-10-18) has no surviving recipe: on the A-C13 pass's 755 paired "
        "genomes the cells-table nn1 mean is 0.7817 (prints 0.782); the pre-registered per_genome.jsonl "
        "(10-16), its .bak-1015 and per_genome ec_f1 for prFBA were all tried and give 0.782 too. The two arms "
        "tie either way ('the two tie on ECs' holds)"),
    "pheno_nn1_ecf1_p": (
        "no paired test tried reproduces 0.64 (Wilcoxon wilcox/pratt/zsplit 0.45/0.44/0.52, paired t 0.67, "
        "sign test 0.49; also on rxn space, f1_accessory, jaccard, mcc and both pre-registered per_genome "
        "files). The recipe keeps the Wilcoxon signed-rank test the phenotype contrasts use; replace the printed "
        "p with this recipe's value (pre 0.45) or drop it"),
}


def check_snippets(numbers, repo, tag):
    """printed values must still be in the manuscript at PAPER_PRE_COMMIT / the deck at the tag."""
    try:
        tex = subprocess.run(["git", "-C", PAPER, "show", "%s:prFBA.tex" % PAPER_PRE_COMMIT],
                             check=True, capture_output=True, text=True).stdout
    except (subprocess.CalledProcessError, FileNotFoundError):
        tex = None
    try:
        deck = subprocess.run(["git", "-C", repo, "show", "%s:presentations/make_deck.py" % tag],
                              check=True, capture_output=True, text=True).stdout
    except subprocess.CalledProcessError:
        deck = None
    norm = (lambda x: re.sub(r"\s+", " ", x)) if tex is not None else None
    ntex = norm(tex) if tex is not None else None
    for n in numbers:
        if n.get("tex_snippet"):
            n["printed_found_in_tex"] = (norm(n["tex_snippet"]) in ntex) if tex is not None else None
    if deck is not None:
        found = {
            "XPHYLUM_BENCH": "XPHYLUM_BENCH = 6.4" in deck,
            "MAG_REACH": "MAG_REACH = (6011, 20.3, 45.1)" in deck,
            "BAND_MAG": "BAND_MAG = (0.082, 0.202)" in deck,
            "AF_SWEEP": all(s in deck for s in ("0.311)", "0.316)", "0.379)", "0.452)", "0.567)", "0.655)")),
        }
        for n in numbers:
            if n.get("deck_constant"):
                n["printed_found_in_deck"] = found.get(n["deck_constant"].split("[")[0])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--state", choices=("pre", "post"), required=True)
    ap.add_argument("--repo", default=R_DEFAULT)
    ap.add_argument("--w", default=W_DEFAULT, help="durable re-score work dir (pre_repair/, start/)")
    ap.add_argument("--tag", default=TAG)
    ap.add_argument("--out", default=None)
    ap.add_argument("--only", default=None, help="comma-separated group names: " +
                    ",".join(g for g, _ in GROUPS))
    a = ap.parse_args(argv)
    if a.state == "pre" and not os.path.isdir(os.path.join(a.w, "start")):
        raise SystemExit("[derived] %s/start is missing: git worktree add %s/start %s" % (a.w, a.w, a.tag))
    res = Resolver(a.state, a.repo, a.w, a.tag)
    only = set(a.only.split(",")) if a.only else None
    cache, numbers, timing = {}, [], {}
    import time
    for name, fn in GROUPS:
        if only and name not in only:
            continue
        t0 = time.time()
        numbers += fn(res, cache) if fn.__code__.co_argcount == 2 else fn(res)
        timing[name] = round(time.time() - t0, 1)
        print("[derived] %-14s %6.1fs" % (name, timing[name]), file=sys.stderr, flush=True)
    check_snippets(numbers, a.repo, a.tag)
    for n in numbers:
        if a.state == "pre" and n["match"] is False:
            n["mismatch_reason"] = KNOWN_PRE_MISMATCH.get(n["key"]) or "unexplained"
    head = subprocess.run(["git", "-C", a.repo, "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    mism = [n["key"] for n in numbers if n["match"] is False]
    unprinted = [n["key"] for n in numbers if n["match"] is None]
    doc = {"_meta": dict(state=a.state, generated=datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
                         repo=a.repo, repo_head=head, tag=a.tag, w=a.w, paper_commit_printed=PAPER_PRE_COMMIT,
                         generator="manuscript_numbers/derived_numbers.py", seconds=timing,
                         inputs=res.inputs, warnings=sorted(set(res.warnings))),
           "summary": dict(n_numbers=len(numbers), n_match=len(numbers) - len(mism) - len(unprinted),
                           n_mismatch=len(mism), mismatches=mism,
                           n_no_printed_value=len(unprinted), no_printed_value=unprinted,
                           mismatches_explained={k: KNOWN_PRE_MISMATCH[k] for k in mism if k in KNOWN_PRE_MISMATCH},
                           printed_not_found_in_tex=[n["key"] for n in numbers if n.get("printed_found_in_tex") is False],
                           printed_not_found_in_deck=[n["key"] for n in numbers if n.get("printed_found_in_deck") is False]),
           "numbers": numbers}
    out = a.out or os.path.join(a.repo, "manuscript_numbers", "derived_numbers.%s.json" % a.state)
    tmp = out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump(doc, fh, indent=1, sort_keys=False, default=str)
        fh.write("\n")
    os.replace(tmp, out)
    print("%-52s %-12s %-12s %-5s %s" % ("key", "printed", "computed", "match", "tex lines"))
    for n in numbers:
        print("%-52s %-12s %-12s %-5s %s" % (n["key"][:52], n["printed"], n["formatted"],
                                             {True: "ok", False: "DIFF", None: "-"}[n["match"]],
                                             n["tex_lines"] or n["deck_constant"]))
    print("\n%d numbers, %d match, %d differ%s; %d have no printed value" % (
        len(numbers), len(numbers) - len(mism) - len(unprinted), len(mism),
        (": " + ", ".join(mism)) if mism else "", len(unprinted)))
    if doc["summary"]["printed_not_found_in_tex"]:
        print("printed snippet not found in prFBA.tex@%s: %s" % (PAPER_PRE_COMMIT,
                                                                 ", ".join(doc["summary"]["printed_not_found_in_tex"])))
    if res.warnings:
        print("%d inputs read from the live tree (git-ignored, not in pre_repair)" % len(set(res.warnings)))
    print("wrote", out)
    return doc


if __name__ == "__main__":
    main()
