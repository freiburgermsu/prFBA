#!/usr/bin/env python
"""Gene capture at the 0.995 self-exclusion threshold, beside 0.987 (2026-11-11).

score_taxacc.py already reports exclude_0.995; score_genecap.py reports only
exclude_0.987. This driver scores both thresholds in one pass by calling
score_genecap's own functions (exclude_self_record, select_representatives via
score_amplicon, score_one, summarize and the writers) unchanged, so the rule is
identical: drop hits with genome_id == source OR taxon_id == source OR
identity >= tau, re-rank, re-select.

What differs from score_genecap.run() is only how the inputs are held in memory,
because the machine is shared and run() loads three multi-GB JSONs whole:
  * the 3.3 GB selection.json is streamed (ijson) and slimmed to the four fields
    score_one reads (tier, n_selected, best_identity, selected[].genome_id);
  * the 4.0 GB hits JSON is streamed record by record, twice (candidates, scoring);
  * the 5.1 GB shared PGFam cache is streamed once and only the genomes the
    scoring can consult are kept, as sets of interned strings.
The provider never fetches: a genome it is asked for but does not hold is
recorded as a miss, fetched, and the scoring pass is repeated until no misses
remain (fixed point), mirroring close_exclude_rep_gap.py.

Genomes absent from the shared cache are fetched into a supplement file, never
cached on failure. --merge-supplement later appends them to the shared cache by a
streaming rewrite (no whole-cache load).

Validation: the include and exclude_0.987 cells are written to repro_0987/ and
must be byte-identical to data/gene_capture_{cells.csv,long.csv,by_region.json}.

Outputs (data/):
  gene_capture_cells_tau0995.csv, gene_capture_long_tau0995.csv,
  gene_capture_by_region_tau0995.json   -- include + exclude_0.995, same schema
  gene_capture_tau_compare.json         -- per region, 0.987 vs 0.995 side by side
  genecap_tau0995/                      -- supplement, repro_0987/, run meta

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import filecmp
import json
import math
import os
import resource
import sys
import time
from concurrent.futures import ThreadPoolExecutor

import ijson

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _common as C  # noqa: E402
from _common import PATHS  # noqa: E402
import score_genecap as S  # noqa: E402

sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import prefetch_gene_families as PF  # noqa: E402
import select_references as SR  # noqa: E402

TAUS = (0.987, 0.995)
MODE87, MODE95 = "exclude_0.987", "exclude_0.995"
DATA = os.path.dirname(PATHS.selection_out)
WORK = os.path.join(DATA, "genecap_tau0995")
SUPP = os.path.join(WORK, "pgfam_supplement.json")
REPRO = os.path.join(WORK, "repro_0987")
OUT_CELLS = os.path.join(DATA, "gene_capture_cells_tau0995.csv")
OUT_LONG = os.path.join(DATA, "gene_capture_long_tau0995.csv")
OUT_SUMMARY = os.path.join(DATA, "gene_capture_by_region_tau0995.json")
OUT_COMPARE = os.path.join(DATA, "gene_capture_tau_compare.json")
OUT_META = os.path.join(WORK, "run_meta.json")
SEL_KEEP = ("tier", "n_selected", "best_identity")
FETCH_WORKERS = 16


def log(m):
    rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576
    print("[genecap-0.995 %s | maxrss %.1f GB] %s" % (time.strftime("%H:%M:%S"), rss, m), flush=True)


# --------------------------------------------------------------------------- #
# Streaming loaders
# --------------------------------------------------------------------------- #
def stream_items(path):
    with open(path, "rb") as fh:
        for k, v in ijson.kvitems(fh, "", use_float=True):
            yield k, v


def slim_selection(rec):
    if not isinstance(rec, dict):
        return rec
    out = {k: rec[k] for k in SEL_KEEP if k in rec}
    if "selected" in rec:
        sel = rec["selected"]
        out["selected"] = ([{"genome_id": s.get("genome_id")} for s in sel]
                           if isinstance(sel, list) else sel)
    return out


def load_selection_slim():
    sel = {}
    for k, v in stream_items(PATHS.selection_out):
        if k != "_meta" and (LIMIT_KEYS is None or k in LIMIT_KEYS):
            sel[k] = slim_selection(v)
    return sel


LIMIT_KEYS = None   # --limit: restrict to the first N hit records (smoke test)


def stream_hits():
    n = 0
    for k, v in stream_items(PATHS.hits_json):
        if k != "_meta" and isinstance(v, dict):
            if LIMIT_KEYS is not None:
                if k not in LIMIT_KEYS:
                    continue
                n += 1
            yield k, v
            if LIMIT_KEYS is not None and n >= len(LIMIT_KEYS):
                return


def first_hit_keys(n):
    out = []
    for k, v in stream_items(PATHS.hits_json):
        if k != "_meta" and isinstance(v, dict):
            out.append(k)
            if len(out) >= n:
                break
    return set(out)


# --------------------------------------------------------------------------- #
# Candidate universe (what the exclude re-selections can consult)
# --------------------------------------------------------------------------- #
def candidates_of(rec):
    res = SR.select_representatives(rec, knobs=S._select_knobs(), gene_provider=None)
    return {str(d["genome_id"]) for d in res["decisions"]
            if d["disposition"] == "selected" or d["stage"] in PF.STAGE7}


def needed_genomes(selection):
    need = set()
    for key, inc in selection.items():
        src_gid, _ = S.split_key(key)
        need.add(src_gid)
        need.update(S.selected_rep_gids(inc or {}))
    n_hits = 0
    for key, hr in stream_hits():
        if key not in selection or "top20" not in hr:
            continue
        n_hits += 1
        src_gid, _ = S.split_key(key)
        src_tid = S._src_taxon_id(src_gid)
        for tau in TAUS:
            rec, _ = S.exclude_self_record(hr, src_gid, src_tid, tau)
            need |= candidates_of(rec)
    need.discard("None")
    return need, n_hits


# --------------------------------------------------------------------------- #
# PGFam provider over a streamed subset of the shared cache + supplement
# --------------------------------------------------------------------------- #
class Provider:
    def __init__(self):
        self.sets = {}
        self.intern = {}
        self.misses = set()

    def _iset(self, fams):
        it = self.intern
        return {it.setdefault(f, f) for f in fams}

    def add(self, gid, fams):
        self.sets[str(gid)] = self._iset(fams)

    def __call__(self, gid):
        gid = str(gid)
        s = self.sets.get(gid)
        if s is None:
            self.misses.add(gid)
            return set()
        return s


def load_supplement():
    return json.load(open(SUPP)) if os.path.exists(SUPP) else {}


def save_supplement(supp):
    tmp = SUPP + ".tmp"
    json.dump(supp, open(tmp, "w"))
    os.replace(tmp, SUPP)


def fill_provider(prov, want, cache_keys_out):
    """Stream the shared cache; keep the genomes in `want`; record every key."""
    n = 0
    for gid, fams in stream_items(PATHS.pgfam_cache):
        cache_keys_out.add(gid)
        if gid in want:
            prov.add(gid, fams)
            n += 1
    return n


def fetch_missing(gids):
    got, err = {}, 0

    def one(g):
        try:
            return g, PF.fetch_pgfams(g)
        except Exception:
            return g, None
    with ThreadPoolExecutor(FETCH_WORKERS) as ex:
        for i, (g, fams) in enumerate(ex.map(one, sorted(gids)), 1):
            if fams is None:
                err += 1                  # leave absent, never cache a failure
            else:
                got[g] = fams
            if i % 200 == 0:
                log("  fetched %d/%d (err=%d)" % (i, len(gids), err))
    return got, err


# --------------------------------------------------------------------------- #
# Scoring pass
# --------------------------------------------------------------------------- #
def _same(a, b):
    return json.dumps(a, sort_keys=True, default=str) == json.dumps(b, sort_keys=True, default=str)


def scoring_pass(selection, prov, cache_keys):
    cells, best_nonself, seen = {}, {}, set()

    def score_key(key, hr):
        inc = selection[key]
        src_gid, _ = S.split_key(key)
        src_tid = S._src_taxon_id(src_gid)
        best_nonself[key] = S.best_nonself_identity(hr, src_gid, src_tid)
        inc_a, e87 = S.score_amplicon(key, inc, hr, prov, cache_keys, tau=0.987)
        inc_b, e95 = S.score_amplicon(key, inc, hr, prov, cache_keys, tau=0.995)
        if not _same(inc_a, inc_b):
            raise SystemExit("include cell differs between the two calls for %s" % key)
        e95["self_mode"] = MODE95
        cells[key] = (inc_a, e87, e95)

    for key, hr in stream_hits():
        if key in selection:
            seen.add(key)
            score_key(key, hr)
    for key in selection:
        if key not in seen:
            score_key(key, None)          # no hit record -> exclude modes no_hits
    return cells, best_nonself


def write_set(cells_by_key, best_nonself, idx, mode, cells_path, long_path, summary_path):
    ordered = []
    for key in sorted(cells_by_key):
        trip = cells_by_key[key]
        ordered.extend([trip[0], trip[idx]])
    S.TAU = 0.987 if mode == MODE87 else 0.995
    S.MODE_EXCLUDE = mode
    S.SELF_MODES = [S.MODE_INCLUDE, mode]
    S.write_cells_csv(cells_path, ordered, best_nonself)
    S.write_long_csv(long_path, ordered, best_nonself)
    summary = S.summarize(ordered, best_nonself)
    if mode == MODE95:
        summary["_meta"]["headline_note"] = summary["_meta"]["headline_note"].replace(
            "exclude_0.987", MODE95).replace("~24% of ", "")
        summary["_meta"]["produced_by"] = (
            "region_validation_d100/scripts/genecap_tau0995.py (score_genecap functions, "
            "tau 0.995); the 0.987 cells it reproduces are byte-identical to "
            "gene_capture_cells.csv / gene_capture_by_region.json")
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=1)
    return summary


def compare_block(s87, s95):
    out = {}
    for region in s87["regions"]:
        r87, r95 = s87["regions"][region], s95["regions"][region]
        row = {}
        for label, blk in (("include", r87[S.MODE_INCLUDE]), (MODE87, r87[MODE87]),
                           (MODE95, r95[MODE95])):
            row[label] = dict(
                outcome_counts=blk["outcome_counts"],
                n_scorable_ok=blk["n_scorable_ok"], n_scorable_pop=blk["n_scorable_pop"],
                recall_resolved=blk["macro"]["recall"], recall_population=blk["macro_pop"]["recall"],
                f1_resolved=blk["macro"]["f1"], f1_population=blk["macro_pop"]["f1"],
                precision_resolved=blk["macro"]["precision"])
        v = lambda m, q: row[m][q]["value"]
        row["delta_0995_minus_0987"] = {
            q: (None if v(MODE95, q) is None or v(MODE87, q) is None
                else round(v(MODE95, q) - v(MODE87, q), 6))
            for q in ("recall_resolved", "recall_population", "f1_resolved", "f1_population")}
        n = lambda m: row[m]["outcome_counts"]
        pop = lambda m: row[m]["n_scorable_pop"]
        row["abstention_rate"] = {m: (n(m)["abstained"] / pop(m) if pop(m) else None)
                                  for m in (MODE87, MODE95)}
        out[region] = row
    return out


# --------------------------------------------------------------------------- #
# Supplement merge (streaming rewrite of the shared cache; run when no other writer)
# --------------------------------------------------------------------------- #
def merge_supplement():
    supp = load_supplement()
    if not supp:
        log("supplement empty; nothing to merge")
        return
    path = PATHS.pgfam_cache
    tmp = path + ".tmp"
    if os.path.exists(tmp):
        raise SystemExit("%s exists: another writer may be active; not merging" % tmp)
    n_old = n_new = 0
    seen = set()
    with open(tmp, "w") as out:
        out.write("{")
        first = True
        for gid, fams in stream_items(path):
            out.write(("" if first else ", ") + json.dumps(gid) + ": " + json.dumps(fams))
            first = False
            seen.add(gid)
            n_old += 1
        for gid in supp:
            if gid in seen:               # never clobber an existing value
                continue
            out.write(("" if first else ", ") + json.dumps(gid) + ": " + json.dumps(supp[gid]))
            first = False
            n_new += 1
        out.write("}")
    # verify the new file parses to old + new keys before replacing
    n_check = sum(1 for _ in stream_items(tmp))
    if n_check != n_old + n_new:
        os.remove(tmp)
        raise SystemExit("merge check failed: %d != %d + %d" % (n_check, n_old, n_new))
    os.replace(tmp, path)
    log("merged %d supplement genomes into the shared cache (%d -> %d)" % (n_new, n_old, n_old + n_new))


# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merge-supplement", action="store_true",
                    help="append the fetched supplement to the shared cache and exit")
    ap.add_argument("--limit", type=int, default=0,
                    help="smoke test on the first N hit records; writes under "
                         "genecap_tau0995/smoke/ and checks the 0.987 rows against "
                         "gene_capture_cells.csv")
    args = ap.parse_args()
    global LIMIT_KEYS, WORK, SUPP, REPRO, OUT_CELLS, OUT_LONG, OUT_SUMMARY, OUT_COMPARE, OUT_META
    if args.limit:
        LIMIT_KEYS = first_hit_keys(args.limit)
        WORK = os.path.join(WORK, "smoke")
        SUPP = os.path.join(WORK, "pgfam_supplement.json")
        REPRO = os.path.join(WORK, "repro_0987")
        OUT_CELLS = os.path.join(WORK, "gene_capture_cells_tau0995.csv")
        OUT_LONG = os.path.join(WORK, "gene_capture_long_tau0995.csv")
        OUT_SUMMARY = os.path.join(WORK, "gene_capture_by_region_tau0995.json")
        OUT_COMPARE = os.path.join(WORK, "gene_capture_tau_compare.json")
        OUT_META = os.path.join(WORK, "run_meta.json")
    os.makedirs(REPRO, exist_ok=True)
    if args.merge_supplement:
        merge_supplement()
        return
    t0 = time.time()
    meta = {"started": time.strftime("%Y-%m-%dT%H:%M:%S%z"), "taus": TAUS,
            "select_mode": "legacy" if S.LEGACY_SELECT else "default_union",
            "inputs": {p: os.path.getsize(p) for p in (PATHS.selection_out, PATHS.hits_json,
                                                        PATHS.pgfam_cache)}}

    log("streaming selection.json")
    selection = load_selection_slim()
    log("%d selection records" % len(selection))

    log("enumerating the genomes the include and both exclude re-selections can consult")
    need, n_hits = needed_genomes(selection)
    log("%d hit records; %d genomes needed" % (n_hits, len(need)))
    meta.update(n_selection=len(selection), n_hit_records=n_hits, n_needed=len(need))

    prov, cache_keys = Provider(), set()
    log("streaming the shared PGFam cache for the needed subset")
    n_kept = fill_provider(prov, need, cache_keys)
    supp = load_supplement()
    for g, fams in supp.items():
        if g not in cache_keys:
            prov.add(g, fams)
            cache_keys.add(g)
    log("shared cache %d genomes; kept %d; supplement %d; interned PGFams %d"
        % (len(cache_keys) - len(supp), n_kept, len(supp), len(prov.intern)))

    missing = sorted(g for g in need if g not in prov.sets)
    rounds = []
    for rnd in range(1, 6):
        if missing:
            log("round %d: fetching %d genomes absent from cache" % (rnd, len(missing)))
            got, err = fetch_missing(missing)
            for g, fams in got.items():
                supp[g] = fams
                prov.add(g, fams)
                cache_keys.add(g)
            save_supplement(supp)
            rounds.append({"round": rnd, "fetched": len(got), "errors": err})
            log("round %d fetched %d (err=%d)" % (rnd, len(got), err))
        prov.misses = set()
        log("scoring pass %d" % rnd)
        cells, best_nonself = scoring_pass(selection, prov, cache_keys)
        misses = sorted(prov.misses)
        log("scoring pass %d: %d provider misses" % (rnd, len(misses)))
        rounds.append({"round": rnd, "provider_misses": len(misses)})
        if not misses:
            break
        # a miss that the shared cache already holds just needs loading
        in_cache = set(misses) & cache_keys
        if in_cache:
            fill_provider(prov, in_cache, set())
        missing = sorted(set(misses) - set(prov.sets))
    else:
        raise SystemExit("no fixed point after 5 rounds")
    meta["rounds"] = rounds
    meta["n_supplement"] = len(supp)

    log("writing the 0.987 reproduction and the 0.995 outputs")
    s87 = write_set(cells, best_nonself, 1, MODE87,
                    os.path.join(REPRO, "gene_capture_cells.csv"),
                    os.path.join(REPRO, "gene_capture_long.csv"),
                    os.path.join(REPRO, "gene_capture_by_region.json"))
    repro = {}
    if LIMIT_KEYS is not None:
        import csv
        ref_rows = [r for r in csv.reader(open(PATHS.genecap_cells_csv))]
        head, ref_rows = ref_rows[0], [r for r in ref_rows[1:] if r[0] in LIMIT_KEYS]
        new_rows = list(csv.reader(open(os.path.join(REPRO, "gene_capture_cells.csv"))))
        repro["cells_rows_match_existing"] = (new_rows[0] == head and new_rows[1:] == ref_rows)
        log("smoke: %d reproduced rows vs %d existing rows" % (len(new_rows) - 1, len(ref_rows)))
    for name, ref in (() if LIMIT_KEYS is not None else (("gene_capture_cells.csv", PATHS.genecap_cells_csv),
                      ("gene_capture_long.csv", PATHS.gene_capture_long_csv),
                      ("gene_capture_by_region.json", PATHS.genecap_summary_json))):
        repro[name] = filecmp.cmp(os.path.join(REPRO, name), ref, shallow=False)
    meta["repro_0987_byte_identical"] = repro
    log("0.987 reproduction byte-identical: %s" % repro)

    s95 = write_set(cells, best_nonself, 2, MODE95, OUT_CELLS, OUT_LONG, OUT_SUMMARY)
    compare = {"_meta": {
        "note": ("Gene capture under self-exclusion at tau 0.987 and 0.995 on the depth-100 "
                 "original panel (exact-tie union). 'resolved' = OK cells only; 'population' "
                 "= OK + abstained + pred_empty with recall 0 for abstentions (the honest "
                 "novel-organism number). The include mode is identical under both."),
        "repro_0987_byte_identical": repro,
        "outcome_counts_total": {
            "include": s87["_meta"]["outcome_counts_total"][S.MODE_INCLUDE],
            MODE87: s87["_meta"]["outcome_counts_total"][MODE87],
            MODE95: s95["_meta"]["outcome_counts_total"][MODE95]},
        "headline_macro_recall_pct_resolved": {
            MODE87: s87["_meta"]["headline_macro_recall_pct_resolved"][MODE87],
            MODE95: s95["_meta"]["headline_macro_recall_pct_resolved"][MODE95]},
        "headline_macro_recall_pct_population": {
            MODE87: s87["_meta"]["headline_macro_recall_pct_population"][MODE87],
            MODE95: s95["_meta"]["headline_macro_recall_pct_population"][MODE95]},
    }, "regions": compare_block(s87, s95)}
    with open(OUT_COMPARE, "w") as fh:
        json.dump(compare, fh, indent=1)

    meta["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S%z")
    meta["elapsed_s"] = round(time.time() - t0)
    meta["maxrss_gb"] = round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1048576, 2)
    json.dump(meta, open(OUT_META, "w"), indent=1)
    h = compare["_meta"]
    log("population recall %% 0.987: %s" % json.dumps(h["headline_macro_recall_pct_population"][MODE87]))
    log("population recall %% 0.995: %s" % json.dumps(h["headline_macro_recall_pct_population"][MODE95]))
    log("done in %.0f s" % (time.time() - t0))
    if not all(repro.values()):
        sys.exit(3)


if __name__ == "__main__":
    main()
