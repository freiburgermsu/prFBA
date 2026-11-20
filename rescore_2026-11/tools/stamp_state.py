#!/usr/bin/env python
"""State stamp for the 2026-11 re-score (plan R1.0 step 2 and R1.6 step 1). Read-only.

    stamp_state.py --out $W/stamps/R1.0.json [--keys-out K.txt] [--repaired-values-out V.json]
                   [--compare $W/stamps/R1.0.json]

Records: the shared PGFam cache's size, mtime, sha256, n_entries and key-set sha256 (sha256 over
each key's bytes + NUL in file order, the method of build_snapshot.stream_restore); optionally the
full key list (one per line, file order) and the current values of the 1,006 repaired genomes
(so R1.6 can prove no repaired entry changed); the E3 catalogue's _meta (empty-genome list as a
count + sha256) and whether it is current against the cache; sha256 of CORRECTIONS_2026-11-13.json;
the prFBA and paper git HEADs; `sha256sum -c region_validation_d100/d20_baseline.sha256`.
--compare prints which of size/mtime/sha256/keys differ from an earlier stamp (exit 0 either way).
"""
import argparse
import hashlib
import json
import os
import re
import subprocess
import sys
import time

R = "/home/freiburger/Documents/prFBA"
P = "/home/freiburger/Documents/prFBA_paper"
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
CACHE = os.path.join(R, "bvbrc_cache", "genome_gene_families.json")
CATALOGUE = os.path.join(R, "bvbrc_cache", "e3_clade_catalogue.json")
CORR = os.path.join(R, "bvbrc_cache", "CORRECTIONS_2026-11-13.json")
ENTRY = re.compile(rb'"([^"]+)": \[([^\]]*)\]')


def sha256_file(path, chunk=64 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(chunk), b""):
            h.update(b)
    return h.hexdigest()


def stream_cache(path, want_values, chunk=256 << 20):
    h_file, h_keys = hashlib.sha256(), hashlib.sha256()
    keys, values, n = [], {}, 0
    carry = b""
    with open(path, "rb") as f:
        while True:
            buf = f.read(chunk)
            h_file.update(buf)
            data = carry + buf
            last = 0
            for m in ENTRY.finditer(data):
                k = m.group(1)
                h_keys.update(k + b"\0")
                ks = k.decode()
                keys.append(ks)
                if ks in want_values:
                    values[ks] = json.loads(b"[" + m.group(2) + b"]")
                n += 1
                last = m.end()
            carry = data[last:]
            if not buf:
                break
    return {"sha256": h_file.hexdigest(), "keys_sha256": h_keys.hexdigest(), "n_entries": n}, keys, values


def git_head(repo):
    try:
        return subprocess.run(["git", "-C", repo, "rev-parse", "HEAD"], capture_output=True, text=True,
                              env=dict(os.environ, GIT_OPTIONAL_LOCKS="0")).stdout.strip()
    except Exception as e:  # pragma: no cover
        return f"error: {e}"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", required=True)
    ap.add_argument("--keys-out")
    ap.add_argument("--repaired-values-out")
    ap.add_argument("--compare")
    a = ap.parse_args()
    t0 = time.time()
    st0 = os.stat(CACHE)
    rep = json.load(open(os.path.join(W, "expected", "repaired_ids.json")))
    repaired = set(rep["repaired_2026_11_13"])
    want = repaired if a.repaired_values_out else set()
    cache_info, keys, values = stream_cache(CACHE, want)
    st1 = os.stat(CACHE)
    if (st0.st_size, st0.st_mtime) != (st1.st_size, st1.st_mtime):
        raise SystemExit("the cache changed while it was being stamped; re-run")
    cache_info.update(path=CACHE, size=st1.st_size, mtime=st1.st_mtime,
                      mtime_iso=time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(st1.st_mtime)),
                      tmp_exists=os.path.exists(CACHE + ".tmp"))
    sys.path.insert(0, R)
    import pfunc_estimator as PE
    meta = PE.catalogue_meta(CATALOGUE)
    empties = meta.pop("empty_pgfam_genomes", None)
    if empties is not None:
        meta["empty_pgfam_genomes_n"] = len(empties)
        meta["empty_pgfam_genomes_sha256"] = hashlib.sha256("\n".join(empties).encode()).hexdigest()
    cst = os.stat(CATALOGUE)
    base = subprocess.run(["sha256sum", "-c", "region_validation_d100/d20_baseline.sha256"], cwd=R,
                          capture_output=True, text=True)
    lines = [ln for ln in base.stdout.splitlines() if ln.strip()]
    out = {
        "stamped": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "cache": cache_info,
        "catalogue": {"path": CATALOGUE, "size": cst.st_size, "mtime": cst.st_mtime,
                      "sha256": sha256_file(CATALOGUE), "meta": meta,
                      "catalogue_is_current": PE.catalogue_meta_is_current(meta, CACHE)},
        "corrections_sha256": sha256_file(CORR),
        "git": {"prFBA": git_head(R), "paper": git_head(P)},
        "d20_baseline": {"rc": base.returncode, "n_ok": sum(ln.endswith(": OK") for ln in lines),
                         "n_lines": len(lines), "failed": [ln for ln in lines if not ln.endswith(": OK")]},
        "n_repaired_values": len(values),
        "seconds": None,
    }
    if a.keys_out:
        with open(a.keys_out, "w") as fh:
            fh.write("\n".join(keys) + "\n")
        out["keys_file"] = os.path.abspath(a.keys_out)
    if a.repaired_values_out:
        missing = sorted(repaired - set(values))
        if missing:
            raise SystemExit(f"{len(missing)} repaired genomes are not in the cache: {missing[:5]}")
        with open(a.repaired_values_out, "w") as fh:
            json.dump(values, fh, sort_keys=True)
        out["repaired_values_file"] = os.path.abspath(a.repaired_values_out)
        out["repaired_values_sha256"] = sha256_file(a.repaired_values_out)
    out["seconds"] = round(time.time() - t0, 1)
    with open(a.out, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({k: out[k] for k in ("cache", "corrections_sha256", "git", "d20_baseline")}, indent=1))
    print("catalogue_is_current:", out["catalogue"]["catalogue_is_current"])
    if a.compare:
        old = json.load(open(a.compare))
        diff = {k: (old["cache"][k], out["cache"][k]) for k in ("size", "mtime", "sha256", "keys_sha256", "n_entries")
                if old["cache"][k] != out["cache"][k]}
        print("CACHE UNCHANGED vs", a.compare if not diff else "")
        if diff:
            print("CACHE CHANGED vs", a.compare, json.dumps(diff, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
