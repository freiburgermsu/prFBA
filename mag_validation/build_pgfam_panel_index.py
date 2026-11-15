#!/usr/bin/env python
"""Rebuild a MAG-arm panel PGFam index from the shared BV-BRC PGFam cache.

    python mag_validation/build_pgfam_panel_index.py --old-index OLD --out NEW [--restore-old-values]

The panel indices (`mag_validation/pgfam_panel.json`, `kmeans_panel/pgfam_panel_kmeans.json`,
`random_panel/pgfam_panel_random.json`; `{"n_vocab": N, "genomes": {gid: [int, ...]}}`, read by
`score_mag.load_pgfams` and `af_gate_counterfactual.py`) were extracted from the shared cache on
2026-10-08/09/18 by a builder that was never committed. This script is its replacement
(re-score plan W1; spec: pipelines_mag_deploy.md section A.4):

* the genome set and the key order come from OLD (the original panel's order matches neither
  selection_10k.json nor benchmark_truth.json, so it cannot be re-derived);
* the cache is streamed with ijson (read only); every panel genome must be in it;
* PGFam strings are sorted within a genome and int-encoded by first appearance over the genomes
  in OLD's order; each list is the sorted unique ints (score_mag only ever takes sets of them);
* --restore-old-values puts back the pre-2026-11-13 values of the repaired genomes
  (`bvbrc_cache/CORRECTIONS_2026-11-13.json` `backup_old_values`): the control index, which must
  score the committed MAG tables byte for byte;
* the output is written to NEW.tmp and renamed; a sidecar NEW.stamp.json records the cache stamp
  before and after (must be equal), the corrections-log sha256, counts, and every genome whose
  list length differs from OLD.

--check runs the W1 checks against OLD and the cache and exits non-zero on any failure:
same key set and order; len(new[g]) == len(set(cache[g])) for every g; lengths equal to OLD
outside the repaired genomes; for the repaired panel genomes, OLD length == the log's old_size and
new length == new_size (or the reverse with --restore-old-values); |A & B| equal between OLD and
NEW on random pairs of unrepaired genomes (the encodings differ only by relabelling).
"""
import argparse
import hashlib
import json
import os
import random
import sys
import time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CACHE = os.path.join(ROOT, "bvbrc_cache", "genome_gene_families.json")
CORR = os.path.join(ROOT, "bvbrc_cache", "CORRECTIONS_2026-11-13.json")


def sha256_file(path, chunk=64 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def stamp(path, with_sha=False):
    st = os.stat(path)
    out = {"path": os.path.abspath(path), "size": st.st_size, "mtime": st.st_mtime,
           "mtime_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(st.st_mtime))}
    if with_sha:
        out["sha256"] = sha256_file(path)
    return out


def read_cache_values(cache_path, want):
    """{gid: list} for the genomes in `want`, streamed from the cache (ijson, read only)."""
    import ijson
    vals = {}
    with open(cache_path, "rb") as fh:
        for g, v in ijson.kvitems(fh, ""):
            if g in want:
                vals[g] = list(v)
    return vals


def load_corrections(corr_path, cache_key):
    c = json.load(open(corr_path))
    backup = c["backup_old_values"].get(cache_key, {})
    sizes = {e["genome"]: (e["old_size"], e["new_size"], e["reason"]) for e in c["entries"]
             if e.get("file") == cache_key and e.get("status") == "applied"}
    return backup, sizes


def encode(keys, vals):
    """First-appearance int encoding over genomes in `keys` order, strings sorted per genome."""
    vocab, genomes = {}, {}
    for g in keys:
        genomes[g] = sorted({vocab.setdefault(f, len(vocab)) for f in sorted(set(vals[g]))})
    return len(vocab), genomes


def build(old_index, out, cache_path=CACHE, corr_path=CORR, restore_old=False, cache_key=None):
    t0 = time.time()
    old = json.load(open(old_index))
    keys = list(old["genomes"])
    want = set(keys)
    cache_key = cache_key or os.path.realpath(cache_path)
    st0 = stamp(cache_path, with_sha=True)
    vals = read_cache_values(cache_path, want)
    missing = sorted(want - set(vals))
    if missing:
        raise SystemExit(f"[abort] {len(missing)} panel genomes missing from the cache: {missing[:5]}")
    backup, sizes = load_corrections(corr_path, cache_key) if os.path.exists(corr_path) else ({}, {})
    restored = []
    if restore_old:
        if not backup:
            raise SystemExit(f"[abort] no backup_old_values for {cache_key} in {corr_path}")
        for g in sorted(want & set(backup)):
            vals[g] = list(backup[g])
            restored.append(g)
    n_vocab, genomes = encode(keys, vals)
    st1 = stamp(cache_path, with_sha=True)
    if (st0["size"], st0["mtime"], st0["sha256"]) != (st1["size"], st1["mtime"], st1["sha256"]):
        raise SystemExit("[abort] the cache changed while the index was built")
    tmp = out + ".tmp"
    with open(tmp, "w") as fh:
        json.dump({"n_vocab": n_vocab, "genomes": genomes}, fh)
    os.replace(tmp, out)
    len_diff = sorted(g for g in keys if len(genomes[g]) != len(old["genomes"][g]))
    side = {
        "what": "MAG-arm panel PGFam index rebuilt from the shared cache "
                + ("with the pre-2026-11-13 values of repaired genomes restored (control)"
                   if restore_old else "(current values)"),
        "tool": os.path.abspath(__file__), "built": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "old_index": stamp(old_index, with_sha=True), "out": stamp(out, with_sha=True),
        "cache_before": st0, "cache_after": st1,
        "corrections": {"path": os.path.abspath(corr_path),
                        "sha256": sha256_file(corr_path) if os.path.exists(corr_path) else None},
        "restore_old_values": restore_old, "n_restored": len(restored),
        "n_genomes": len(keys), "n_vocab": n_vocab, "old_n_vocab": old.get("n_vocab"),
        "n_empty": sum(1 for g in keys if not genomes[g]),
        "repaired_in_panel": sorted(want & set(sizes)),
        "len_differs_from_old": len_diff,
        "seconds": round(time.time() - t0, 1),
    }
    with open(out + ".stamp.json", "w") as fh:
        json.dump(side, fh, indent=1)
    return side


def check(old_index, new_index, cache_path=CACHE, corr_path=CORR, restore_old=False,
          n_pairs=1000, seed=1729, cache_key=None):
    """The W1 checks; returns (ok, report)."""
    old = json.load(open(old_index))["genomes"]
    new = json.load(open(new_index))["genomes"]
    cache_key = cache_key or os.path.realpath(cache_path)
    keys_old, keys_new = list(old), list(new)
    rep = {"same_keys_and_order": keys_old == keys_new}
    want = set(keys_new)
    vals = read_cache_values(cache_path, want)
    backup, sizes = load_corrections(corr_path, cache_key)
    if restore_old:
        for g in want & set(backup):
            vals[g] = backup[g]
    bad_len = [g for g in keys_new if len(new[g]) != len(set(vals.get(g, [])))]
    rep["len_eq_cache_set"] = not bad_len
    rep["len_ne_cache_set_examples"] = bad_len[:10]
    repaired = want & set(sizes)
    rep["n_repaired_in_panel"] = len(repaired)
    outside = [g for g in keys_new if g not in repaired and len(new[g]) != len(old[g])]
    rep["len_eq_old_outside_repaired"] = not outside
    rep["len_ne_old_outside_repaired_examples"] = outside[:10]
    # the repaired genomes: OLD holds old_size; NEW holds new_size (or old_size when restored)
    wrong = []
    for g in sorted(repaired):
        o, n, _ = sizes[g]
        exp_new = o if restore_old else n
        if len(old[g]) != o or len(new[g]) != exp_new:
            wrong.append({"genome": g, "old_len": len(old[g]), "log_old_size": o,
                          "new_len": len(new[g]), "expected_new": exp_new})
    rep["repaired_sizes_eq_log"] = not wrong
    rep["repaired_sizes_mismatch"] = wrong[:10]
    rep["n_len_eq_old_among_repaired"] = sum(1 for g in repaired if len(new[g]) == len(old[g]))
    # relabelling invariance: |A & B| on random unrepaired pairs
    rng = random.Random(seed)
    pool = [g for g in keys_new if g not in repaired]
    bad_pairs = 0
    for _ in range(n_pairs):
        a, b = rng.sample(pool, 2)
        if len(set(old[a]) & set(old[b])) != len(set(new[a]) & set(new[b])):
            bad_pairs += 1
    rep["n_pairs"] = n_pairs
    rep["intersection_invariant"] = bad_pairs == 0
    rep["n_pairs_differing"] = bad_pairs
    ok = all(rep[k] for k in ("same_keys_and_order", "len_eq_cache_set", "len_eq_old_outside_repaired",
                              "repaired_sizes_eq_log", "intersection_invariant"))
    return ok, rep


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--old-index", required=True, help="the existing index: genome set and order")
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache", default=CACHE)
    ap.add_argument("--corrections", default=CORR)
    ap.add_argument("--cache-key", help="the cache's key in the corrections log (default: realpath of --cache)")
    ap.add_argument("--restore-old-values", action="store_true",
                    help="put back the pre-2026-11-13 values of repaired genomes (control index)")
    ap.add_argument("--check", action="store_true", help="run the W1 checks after building")
    ap.add_argument("--check-only", action="store_true", help="only run the checks on an existing --out")
    ap.add_argument("--pairs", type=int, default=1000)
    a = ap.parse_args(argv)
    if os.path.realpath(a.out) == os.path.realpath(a.old_index):
        ap.error("--out must differ from --old-index (move the old index aside after checking)")
    if not a.check_only:
        side = build(a.old_index, a.out, a.cache, a.corrections, a.restore_old_values, a.cache_key)
        print(f"[panel-index] {a.out}: {side['n_genomes']} genomes, n_vocab {side['n_vocab']} "
              f"(old {side['old_n_vocab']}), {side['n_empty']} empty, {side['n_restored']} restored, "
              f"{len(side['len_differs_from_old'])} lengths differ from the old index "
              f"({side['seconds']}s)", flush=True)
    if a.check or a.check_only:
        ok, rep = check(a.old_index, a.out, a.cache, a.corrections, a.restore_old_values, a.pairs,
                        cache_key=a.cache_key)
        with open(a.out + ".check.json", "w") as fh:
            json.dump(dict(rep, ok=ok), fh, indent=1)
        print(f"[panel-index] checks {'PASS' if ok else 'FAIL'}: {json.dumps(rep)}", flush=True)
        return 0 if ok else 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
