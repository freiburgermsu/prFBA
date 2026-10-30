#!/usr/bin/env python
"""
build_d100.py -- run build_synthetic_genomes_parallel.py unchanged in its selection and
merge logic, with four run-time patches the depth-100 deployment build needs:

 1. BV-BRC User-Agent. build_synthetic_genomes._api_json (the 'fast' genome builder used
    by the cache warm step) calls urllib with the bare Python User-Agent, which BV-BRC now
    answers with HTTP 403 (re-verified 2026-10-29: bare UA -> 403, the header below -> 200;
    prefetch_gene_families.fetch_pgfams already sends it). Without this every selected genome
    missing from bvbrc_cache/kbase_genome_cache fails to fetch and is silently dropped from
    its ASV's merge. The patch sends the same header as prefetch_gene_families.
 2. --tieset SIDE.json (from tieset_gene_cache.py): wraps gene_count_fn and counts queries
    for genomes outside the tie set; if any occur the small cache could differ from the
    shared one, so the run ABORTS before writing the selection.
 3. gc.freeze() before the merge ProcessPool forks, so the workers' cyclic GC does not touch
    (and copy-on-write duplicate) the parent's hits/audit objects -- at depth 100 those are
    ~5x the 20-hit size. Output is unaffected.
 4. Atomic genome-cache writes. build_synthetic_genomes.build_genome writes
    bvbrc_cache/kbase_genome_cache/{gid}.json in place; a kill mid-write (e.g. MemoryMax) would
    leave a truncated file in that SHARED cache, which every later build would fail to parse and
    silently drop from its merge. The patch writes {gid}.json.tmp.<pid>.<thread> and os.replace()s
    it; reads and the fetched content are unchanged.

Every other argument goes to build_synthetic_genomes_parallel.py as-is (see its --help).

    python deployment_d100/build_d100.py --tieset GC.tieset.json \
        --hits H --gene-cache GC --out-dir OUT --selection-out SEL --procs 8 --io-workers 8

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import gc
import json
import os
import subprocess
import sys
import threading
import time
import urllib.request

PRFBA = "/home/freiburger/Documents/prFBA"
sys.path.insert(0, PRFBA)

import build_synthetic_genomes as BSG            # noqa: E402
import build_synthetic_genomes_parallel as BSGP  # noqa: E402
import select_references as S                    # noqa: E402

UA = {"User-Agent": "prFBA/1.0 (academic; metabolic-model benchmarking)"}  # = prefetch_gene_families


def _api_json_ua(core, query, timeout=60, retries=3):
    url = f"{BSG.API}/{core}/?{query}&http_accept=application/json"
    last = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as resp:
                return json.load(resp)
        except Exception as exc:
            last = exc
            time.sleep(1.5 * (attempt + 1))
    raise last


def _build_genome_atomic(gid, cache_dir, builder):
    """build_synthetic_genomes.build_genome with an atomic cache write (same reads/fetches)."""
    os.makedirs(cache_dir, exist_ok=True)
    path = os.path.join(cache_dir, f"{gid}.json")
    if os.path.exists(path):
        return BSG._load(path)
    g = BSG.fast_build_genome(gid) if builder == "fast" else \
        BSG.KB.BVBRCToKBaseConverter(gid, base_url=BSG.API).build_kbase_genome()
    g.setdefault("id", str(gid))
    tmp = f"{path}.tmp.{os.getpid()}.{threading.get_ident()}"
    with open(tmp, "w") as f:
        json.dump(g, f)
    os.replace(tmp, path)
    return g


def _pop_arg(argv, name):
    if name in argv:
        i = argv.index(name)
        val = argv[i + 1]
        del argv[i:i + 2]
        return val
    return None


def main():
    argv = sys.argv[1:]
    tieset_path = _pop_arg(argv, "--tieset")
    sys.argv = [os.path.join(PRFBA, "build_synthetic_genomes_parallel.py")] + argv

    # provenance of the code that runs
    try:
        head = subprocess.run(["git", "-C", PRFBA, "rev-parse", "HEAD"], capture_output=True,
                              text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", PRFBA, "status", "--porcelain", "--",
                                "build_synthetic_genomes.py", "build_synthetic_genomes_parallel.py",
                                "select_references.py", "prefetch_gene_families.py"],
                               capture_output=True, text=True).stdout.strip()
    except Exception:
        head, dirty = "?", "?"
    print(f"[build_d100] prFBA HEAD {head} | pipeline files dirty: {dirty or 'no'} | "
          f"merge module {BSG.KB.__file__}", flush=True)

    # 1. BV-BRC User-Agent for the fast genome builder
    BSG._api_json = _api_json_ua
    print("[build_d100] patched build_synthetic_genomes._api_json with the BV-BRC User-Agent", flush=True)
    BSG.build_genome = _build_genome_atomic
    print("[build_d100] patched build_synthetic_genomes.build_genome: atomic genome-cache writes", flush=True)

    # 2. tie-set guard on gene_count_fn
    stats = {"queries": 0, "outside": 0, "examples": []}
    if tieset_path:
        tie = set(json.load(open(tieset_path))["tieset"])
        orig = S.gene_count_from_provider

        def guarded(provider):
            f = orig(provider)

            def count(gid):
                stats["queries"] += 1
                if str(gid) not in tie:
                    stats["outside"] += 1
                    if len(stats["examples"]) < 10:
                        stats["examples"].append(str(gid))
                return f(gid)
            return count
        S.gene_count_from_provider = guarded
        orig_exact = BSGP.exact_select

        def exact_select_checked(*a, **k):
            audit, summary = orig_exact(*a, **k)
            print(f"[build_d100] tie-set guard: {stats['queries']} gene_count_fn queries, "
                  f"{stats['outside']} outside the {len(tie)}-genome tie set", flush=True)
            if stats["outside"]:
                raise SystemExit(f"[abort] gene_count_fn queried {stats['outside']} genomes outside the "
                                 f"tie set (e.g. {stats['examples']}): the small cache is not "
                                 f"equivalent to the shared one; nothing was written")
            return audit, summary
        BSGP.exact_select = exact_select_checked
        print(f"[build_d100] tie-set guard on ({len(tie)} genomes from {tieset_path})", flush=True)

    # 3. freeze the parent's objects before the merge pool forks
    _PPE = BSGP.ProcessPoolExecutor

    def ppe_frozen(*a, **k):
        gc.collect()
        gc.freeze()
        print(f"[build_d100] gc.freeze() before forking the merge pool "
              f"({gc.get_freeze_count()} objects frozen)", flush=True)
        return _PPE(*a, **k)
    BSGP.ProcessPoolExecutor = ppe_frozen

    rc = BSGP.main()
    print(f"[build_d100] modelseedpy imported: {'modelseedpy' in sys.modules}", flush=True)
    return rc


if __name__ == "__main__":
    sys.exit(main())
