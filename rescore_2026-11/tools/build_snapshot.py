#!/usr/bin/env python
"""Rebuild the pre-repair (2026-11-11) state of prFBA's shared PGFam cache from the repaired cache
and the 2026-11-13 corrections log (rescore plan W2 / R1.1).

The 11-13 repair (`r4_apply.stream_rewrite`) copied the cache entry by entry, byte for byte, and
replaced the value of 1,006 keys with `json.dumps(new)[1:-1]` inside the same `"key": [` ... `]`
frame. It logged every old value (`backup_old_values[<cache path>]`, original order) and the old
file's size, mtime and sha256. This script inverts that rewrite: it streams the repaired cache with
the same entry regex and puts `json.dumps(old)[1:-1]` back, so the output should be byte-identical
to the pre-repair file (self-test: its sha256 equals the log's `before.sha256`). If the bytes ever
differ, it falls back to an entry-level check (ijson lockstep: every key in the same order, the
1,006 equal to the backup, every other value equal to the repaired cache) and says so.

`--reset-keys FILE` additionally sets the listed keys to `[]` (the d20 / 09-17-era variant
`$SNAP20`: the 20 keys of `bvbrc_cache/poisoned_keys.json` that the 09-17 d20 chain saw as empty
lists and the 10-09 poisoned-key repair later filled). That output has no logged hash, so it is
always verified at entry level.

The output mtime is set to the logged pre-repair mtime (`--mtime`, default the log's
`before.mtime`), so `pfunc_estimator.catalogue_is_current` (size + local-time mtime to the second)
accepts the pre-repair E3 catalogue against it in control C-E3. A sidecar `<out>.stamp.json`
records the inputs' stamps, the replaced and reset keys, and the verification result.

Read-only on the cache: its size, mtime and sha256 are stamped before and after and must agree.

    build_snapshot.py --cache $CACHE --corrections $R/bvbrc_cache/CORRECTIONS_2026-11-13.json --out $SNAP
    build_snapshot.py --cache $CACHE --corrections ... --out $SNAP20 --d20 \
        --reset-keys $W/expected/d20_poison_20.json
"""
import argparse
import hashlib
import itertools
import json
import os
import re
import sys
import time

# the regex of r4_apply / deployment_d100/tieset_gene_cache._ENTRY: one cache entry
ENTRY = re.compile(rb'"([^"]+)": \[([^\]]*)\]')


def sha256_file(path, chunk=64 << 20):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def stamp(path, with_sha=True):
    st = os.stat(path)
    out = {"path": os.path.abspath(path), "size": st.st_size, "mtime": st.st_mtime,
           "mtime_iso": time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(st.st_mtime))}
    if with_sha:
        out["sha256"] = sha256_file(path)
    return out


def stream_restore(src, dst, repl, chunk=256 << 20):
    """Copy src to dst entry by entry; the value of each key in `repl` becomes json.dumps(repl[k]).
    Returns counts and hashes. Aborts on unexpected bytes between entries or a duplicate key."""
    h_out = hashlib.sha256()
    keys_h = hashlib.sha256()
    n, seen = 0, set()
    with open(src, "rb") as f, open(dst, "wb") as out:
        first = f.read(1)
        if first != b"{":
            raise SystemExit(f"[abort] {src}: not a JSON object")
        out.write(first)
        h_out.update(first)
        carry, expect_sep = b"", b""
        while True:
            buf = f.read(chunk)
            data = carry + buf
            last = 0
            for m in ENTRY.finditer(data):
                sep = data[last:m.start()]
                if sep != expect_sep:
                    raise SystemExit(f"[abort] {src}: unexpected bytes {sep[:60]!r} before entry {n}")
                expect_sep = b", "
                gid = m.group(1).decode()
                keys_h.update(m.group(1) + b"\0")
                n += 1
                if gid in repl:
                    if gid in seen:
                        raise SystemExit(f"[abort] {src}: duplicate key {gid}")
                    seen.add(gid)
                    piece = sep + data[m.start():m.start(2)] + json.dumps(repl[gid])[1:-1].encode() + b"]"
                else:
                    piece = data[last:m.end()]
                out.write(piece)
                h_out.update(piece)
                last = m.end()
            carry = data[last:]
            if not buf:
                break
        if carry != b"}":
            raise SystemExit(f"[abort] {src}: trailing bytes {carry[:80]!r}")
        out.write(carry)
        h_out.update(carry)
        out.flush()
        os.fsync(out.fileno())
    missing = sorted(set(repl) - seen)
    if missing:
        raise SystemExit(f"[abort] {len(missing)} keys to restore are not in {src}: {missing[:5]}")
    return {"n_entries": n, "sha256": h_out.hexdigest(), "keys_sha256": keys_h.hexdigest()}


def lockstep_check(snap, cache, repl):
    """Entry-level equality: same keys in the same order; keys in repl hold repl[k] in snap; every
    other value is equal in both files."""
    import ijson
    end = object()
    n = n_repl = 0
    with open(snap, "rb") as fs, open(cache, "rb") as fc:
        it = itertools.zip_longest(ijson.kvitems(fs, "", use_float=True),
                                   ijson.kvitems(fc, "", use_float=True), fillvalue=(end, None))
        for (ks, vs), (kc, vc) in it:
            if ks is end or kc is end or ks != kc:
                return {"ok": False, "error": f"key mismatch at entry {n}: {ks!r} vs {kc!r}"}
            if ks in repl:
                if vs != repl[ks]:
                    return {"ok": False, "error": f"restored key {ks} != backup value"}
                n_repl += 1
            elif vs != vc:
                return {"ok": False, "error": f"unchanged key {ks} differs"}
            n += 1
    return {"ok": True, "n_entries": n, "n_restored_or_reset": n_repl}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cache", required=True, help="the repaired shared PGFam cache (read only)")
    ap.add_argument("--corrections", required=True, help="bvbrc_cache/CORRECTIONS_2026-11-13.json")
    ap.add_argument("--out", required=True)
    ap.add_argument("--reset-keys", help="JSON list of keys to set to [] after restoring (d20 variant)")
    ap.add_argument("--d20", action="store_true",
                    help="the 09-17-era variant: requires --reset-keys; no logged sha256 target")
    ap.add_argument("--expect-sha256", help="target sha256 (default: the log's before.sha256; "
                                            "none with --reset-keys)")
    ap.add_argument("--mtime", type=float, help="output mtime (default: the log's before.mtime)")
    ap.add_argument("--log-key", help="the backup_old_values / files key (default: realpath of --cache)")
    ap.add_argument("--skip-entry-check", action="store_true",
                    help="skip the ijson lockstep check when the sha256 already matches")
    a = ap.parse_args(argv)
    if a.d20 and not a.reset_keys:
        ap.error("--d20 needs --reset-keys (the 20 poisoned keys the 09-17 chain saw as [])")

    t0 = time.time()
    cache = os.path.realpath(a.cache)
    log_key = a.log_key or cache
    corr = json.load(open(a.corrections))
    finfo = corr["files"][log_key]
    backup = corr["backup_old_values"][log_key]
    before = finfo["before"]
    n_by_reason = {}
    for e in corr["entries"]:
        if e.get("file") == log_key and e.get("status") == "applied":
            n_by_reason[e["reason"]] = n_by_reason.get(e["reason"], 0) + 1
    if len(backup) != finfo["n_entries_replaced"]:
        raise SystemExit(f"[abort] backup has {len(backup)} keys, log says {finfo['n_entries_replaced']}")

    repl = {g: list(v) for g, v in backup.items()}
    reset = []
    if a.reset_keys:
        reset = list(json.load(open(a.reset_keys)))
        clash = sorted(set(reset) & set(repl))
        if clash:
            raise SystemExit(f"[abort] reset keys overlap the 11-13 backup: {clash[:5]}")
        for g in reset:
            repl[g] = []
    expect = a.expect_sha256 or (None if reset else before["sha256"])
    mtime = a.mtime if a.mtime is not None else before["mtime"]

    print(f"[snapshot] stamping {cache}", flush=True)
    st_cache0 = stamp(cache)
    after = finfo["after"]
    cache_is_logged_after = (st_cache0["size"], st_cache0["sha256"]) == (after["size"], after["sha256"])
    if not cache_is_logged_after:
        print(f"[warn] the cache is not the logged post-repair file (size {st_cache0['size']} sha "
              f"{st_cache0['sha256'][:12]} vs log {after['size']} {after['sha256'][:12]}); keys added "
              f"since 11-13 will be carried into the snapshot", flush=True)

    tmp = a.out + ".tmp"
    print(f"[snapshot] restoring {len(backup)} backup values"
          + (f" and resetting {len(reset)} keys to []" if reset else "") + f" -> {a.out}", flush=True)
    info = stream_restore(cache, tmp, repl)
    st_cache1 = stamp(cache)
    if (st_cache1["size"], st_cache1["mtime"], st_cache1["sha256"]) != \
            (st_cache0["size"], st_cache0["mtime"], st_cache0["sha256"]):
        os.remove(tmp)
        raise SystemExit("[abort] the cache changed while the snapshot was written")

    sha_ok = (info["sha256"] == expect) if expect else None
    entry = None
    if sha_ok is not True or not a.skip_entry_check:
        print("[snapshot] entry-level lockstep check (ijson) ...", flush=True)
        entry = lockstep_check(tmp, cache, repl)
        if not entry["ok"]:
            os.remove(tmp)
            raise SystemExit(f"[abort] entry-level check failed: {entry['error']}")
    if expect and not sha_ok:
        print(f"[warn] sha256 {info['sha256']} != target {expect}; falling back to entry-level "
              f"equality (passed: {entry})", flush=True)

    os.replace(tmp, a.out)
    ns = int(round(mtime * 1e9))
    os.utime(a.out, ns=(ns, ns))
    st_out = os.stat(a.out)
    mt_iso = time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(st_out.st_mtime))
    want_iso = time.strftime("%Y-%m-%dT%H:%M:%S%z", time.localtime(mtime))
    if mt_iso != want_iso:
        raise SystemExit(f"[abort] mtime {mt_iso} != {want_iso}")

    side = {
        "what": "pre-repair snapshot of the shared PGFam cache, rebuilt by inverting the 2026-11-13 repair"
                + (" plus the d20 (09-17-era) reset of poisoned keys" if reset else ""),
        "tool": os.path.abspath(__file__),
        "argv": sys.argv if argv is None else ["build_snapshot.py"] + list(argv),
        "built": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "cache_before": st_cache0, "cache_after": st_cache1,
        "cache_is_logged_post_repair_file": cache_is_logged_after,
        "corrections": {"path": os.path.abspath(a.corrections), "sha256": sha256_file(a.corrections)},
        "n_restored": len(backup), "restored_by_reason": n_by_reason,
        "reset_keys": sorted(reset), "reset_keys_file": os.path.abspath(a.reset_keys) if a.reset_keys else None,
        "out": {"path": os.path.abspath(a.out), "size": st_out.st_size, "sha256": info["sha256"],
                "mtime": st_out.st_mtime, "mtime_iso": mt_iso, "n_entries": info["n_entries"],
                "keys_sha256": info["keys_sha256"]},
        "target_sha256": expect, "sha256_equals_target": sha_ok,
        "logged_before": before,
        "size_equals_logged_before": (st_out.st_size == before["size"]) if not reset else None,
        "entry_level_check": entry,
        "seconds": round(time.time() - t0, 1),
    }
    with open(a.out + ".stamp.json", "w") as fh:
        json.dump(side, fh, indent=1)
    verdict = ("BYTE-IDENTICAL to the logged pre-repair file" if sha_ok else
               "entry-level equal (no byte target)" if expect is None else
               "entry-level equal only (bytes differ from the logged file)")
    print(f"[snapshot] {verdict}: {a.out} size {st_out.st_size} sha256 {info['sha256']} "
          f"mtime {mt_iso} ({side['seconds']}s)", flush=True)
    return 0 if (sha_ok or expect is None or (entry and entry["ok"])) else 1


if __name__ == "__main__":
    sys.exit(main())
