#!/usr/bin/env python
"""Close the PGFam cache gap on the k-means run's EXCLUDE-mode reps.

Why there is a gap at all: the S8 candidate prefetch covered the genomes named in
the as-is (include-mode) selection -- 39,060 of them, all cached. But exclude mode
does not reuse that selection: score_genecap rebuilds each amplicon's hit record
with self-hits removed and re-runs select_representatives over it, which reaches
deeper into the top-20 and can name genomes the include-mode selection never did.
Those were never prefetched, so the cell scored `pgfam_missing` -- 813 of them for
FullLength16S alone, against 0 in the original panel's run, whose prefetch covered
a wider universe. Left uncorrected it depresses k-means coverage by ~9 points for a
reason that has nothing to do with the selection scheme.

select_representatives consults the gene provider while scoring candidates, so
filling the cache can change which reps it names. The enumeration is therefore
iterated to a fixed point rather than run once.

Writes by value and never persists a failure: a genome whose fetch raises is left
absent (retryable) rather than cached as [], which is the poisoning mode that cost
this project five separate repairs.

Usage: close_exclude_rep_gap.py --tree <tree dir name under ~/Documents/prFBA> [--legacy]
  --legacy enumerates the exclude-mode reps under the LEGACY ordered reducer
  (score_genecap.LEGACY_SELECT = True, i.e. knobs {"select_all": False}) -- use it for a
  tree scored with `score_genecap.py --legacy`; the default enumerates the exact-tie union.
"""
import json, os, sys, time
from concurrent.futures import ThreadPoolExecutor

R = os.path.expanduser("~/Documents/prFBA")
# --tree picks which validation tree's scripts (and so which BASE/selection/hits)
# are imported. Both trees must be scored against the SAME cache state for the
# cross-panel comparison to mean anything, so the gap is closed on each in turn.
TREE = "region_validation_kmeans"
LEGACY = "--legacy" in sys.argv
for i, arg in enumerate(sys.argv):
    if arg == "--tree":
        TREE = sys.argv[i + 1]
if not os.path.isdir(os.path.join(R, TREE, "scripts")):
    sys.exit("[gap] --tree %r: no %s/%s/scripts" % (TREE, R, TREE))
sys.path.insert(0, R)
sys.path.insert(0, R + "/" + TREE + "/scripts")
import prefetch_gene_families as PF
import score_genecap as S
import select_references as SR

# score_genecap._select_knobs() reads LEGACY_SELECT at call time, so setting it here
# makes enumerate_exclude_reps use the same reducer as `score_genecap.py --legacy`.
S.LEGACY_SELECT = LEGACY
CACHE = S.PATHS.pgfam_cache
WORKERS = 16


def log(m): print("[gap:%s] %s" % (TREE.replace("region_validation", "rv"), m), flush=True)


def enumerate_exclude_reps(selection, hits, provider, cache_keys):
    reps = set()
    for key, inc in selection.items():
        hr = hits.get(key)
        if hr is None or "top20" not in hr:
            continue
        src_gid, _ = S.split_key(key)
        src_tid = S._src_taxon_id(src_gid)
        excl_rec, _ = S.exclude_self_record(hr, src_gid, src_tid, S.TAU)
        sel = SR.select_representatives(excl_rec, gene_provider=provider,
                                        knobs=S._select_knobs())
        reps.update(S.selected_rep_gids(sel or {}))
    return reps


def main():
    log("select mode: %s (knobs %s); data %s" % ("legacy" if S.LEGACY_SELECT else "exact-tie union",
                                                  S._select_knobs(), os.path.dirname(S.PATHS.selection_out)))
    log("loading selection + hits ...")
    selection = S.load_selection(S.PATHS.selection_out)
    hits = S.load_hits(S.PATHS.hits_json)
    log("%d amplicons, %d hit records" % (len(selection), len(hits)))

    for rnd in range(1, 5):
        cache_keys = S.load_cache_keys(CACHE)
        provider = SR.bvbrc_gene_provider(CACHE)
        t0 = time.time()
        reps = enumerate_exclude_reps(selection, hits, provider, cache_keys)
        missing = sorted(reps - cache_keys)
        log("round %d: %d exclude-mode reps, %d missing from cache (%.0fs)"
            % (rnd, len(reps), len(missing), time.time() - t0))
        if not missing:
            log("fixed point reached -- no missing exclude-mode reps")
            break

        got, err, empty = {}, 0, 0
        def one(g):
            try:
                return g, PF.fetch_pgfams(g)
            except Exception:
                return g, None
        with ThreadPoolExecutor(WORKERS) as ex:
            for i, (g, fams) in enumerate(ex.map(one, missing), 1):
                if fams is None:
                    err += 1                      # leave ABSENT, never cache []
                else:
                    got[g] = fams
                    if not fams:
                        empty += 1
                if i % 200 == 0:
                    log("  fetched %d/%d (err=%d empty=%d)" % (i, len(missing), err, empty))
        log("round %d fetched %d (err=%d genuinely-empty=%d)" % (rnd, len(got), err, empty))

        if got:
            cache = json.load(open(CACHE))
            added = 0
            for g, fams in got.items():
                if g not in cache:                # never clobber an existing value
                    cache[g] = fams; added += 1
            tmp = CACHE + ".tmp"
            json.dump(cache, open(tmp, "w"))
            os.replace(tmp, CACHE)
            log("round %d merged %d new entries -> cache now %d" % (rnd, added, len(cache)))
        if err and not got:
            log("every fetch failed this round; stopping rather than looping")
            break


if __name__ == "__main__":
    main()
