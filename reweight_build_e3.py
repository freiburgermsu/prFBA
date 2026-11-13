#!/usr/bin/env python
"""
reweight_build_e3.py -- apply the E3 feature-probability estimator to an EXISTING synthetic-
genome build, copy-on-write: the source build is only read; every output goes to a NEW
directory, laid out like a build (so upload_kbase_genomes.py / write_mapping_provenance.py /
the model loaders can use it directly).

What it checks first (run level, e3 mode, before any audit streaming or PGFam work):
  * the E3 clade catalogue must exist and match the PGFam cache (else it stops and prints the
    rebuild command; --pfunc-allow-stale-catalogue overrides);
  * the source build must be a k/n build (a manifest _meta.pfunc naming E3 is refused);
  * the source build must be BV-BRC PATRIC-only (manifest _meta.annotation ==
    build_synthetic_genomes.ANNOTATION_STAMP, in every mode). Builds made at or before git
    3ce12d0 -- the shipped synthetic_genomes_d100 deployments and their synthetic_genomes_d100_e3
    re-weightings among them -- merged PATRIC and the imported RefSeq CDSs and are refused
    (nothing written): rebuild them with the current builder instead. The source genome dicts
    are read BV-BRC PATRIC-only (load_sources), as the builders read them. Every real build
    without the selection stamp below also predates the annotation stamp, so the 'legacy build'
    rules below (and --allow-audit-newer / --allow-cache-newer) no longer apply to any existing
    build; they are kept for a stamped manifest that lost its input stamps;
  * SELECTION AUDIT. A build that stamped its audit (manifest _meta.selection_stamp, written by
    every builder since 2026-11-12, in every --pfunc mode): the sha256 of the audit now must
    equal the stamp. On a mismatch EVERY built ASV in scope (n = 1 included) is copied byte for
    byte (k/n kept) with reason 'inputs drifted since build: selection audit', no audit record
    is streamed and no PGFam store is built. A legacy build (no stamp): the run STOPS before
    writing anything when the audit's mtime is newer than the source manifest's (it may have
    been regenerated since the build), unless --allow-audit-newer (then a warning, and
    _meta.pfunc.input_checks.selection_audit records the override). The audit's content is
    covered only by these two checks; the per-ASV checks below compare the selected genomes,
    not decisions[] (the neighbour list).
Then, per built ASV of the source manifest (in manifest order):
  * n = the manifest's n_genomes (= built source genomes, the n of k/n; cross-checked against
    provenance/{asv}.json 'genomes').
  * n < --pfunc-min-n (default 2), or no estimate possible -> {asv}.json, genome_objects/,
    provenance/ and taxonomies/ files are COPIED byte for byte.
  * per-ASV input-drift checks for n >= min_n (else COPIED, k/n kept, reason 'inputs drifted
    since build: <which>', counted in _meta.pfunc.n_inputs_drifted), in this order:
      - 'selected genomes': the audit record's selected genome ids must equal the provenance
        'genomes' plus the entry's failed[] ids (as multisets);
      - every built source genome must be in --genome-cache-dir (else reason 'source genome(s)
        missing from the genome cache: [...]');
      - 'function_prevalence': build_synthetic_genomes.function_prevalence of those cached
        dicts must equal the provenance's;
      - 'source digest': when the build has pfunc/source_digests.json (every build since
        2026-11-12, every mode), pfunc_estimator.source_digest of each cached dict must equal
        the recorded one (a genome missing from the file, or recorded as null, counts as
        drift); the drifted ids are in the entry's pfunc 'drift_genomes';
      - 'source genome cache newer than build': a legacy build without that file instead
        fails an ASV when any of its source-genome cache files has an mtime newer than the
        source manifest's, unless --allow-cache-newer (the ids are then listed in the entry's
        pfunc 'cache_newer_allowed_genomes' and counted in _meta.pfunc.input_checks).
    A record absent from the audit: reason 'record not in selection audit'
    (_meta.pfunc.n_record_not_in_audit).
  * otherwise E3 (pfunc_estimator.reweight_asv, the same call build_synthetic_genomes.py
    makes at build time; clade prior anchored on the first built rep) replaces k/n on the
    features with a PGFam, in BOTH copies ({asv}.json and genome_objects/{asv}.json, features
    and their '{fid}.CDS' entries; nothing else changes), provenance/{asv}.json gains
    'probability_kn' (the replaced k/n) and 'pfunc' (diagnostics), the manifest entry gains a
    'pfunc' summary.
  * an ASV that raises is COPIED unchanged (k/n kept), reason 'error: ...', counted in
    _meta.pfunc.n_errors. After every ASV (raised or not) the worker parses its {asv}.json,
    genome_objects/ and provenance/ output files; an ASV with one missing or unparseable (a
    source build that was already broken) becomes status 'failed' with the reason in its
    'error' (_meta.pfunc.n_built_demoted_to_failed; n_built_demoted_unreadable_output counts
    those that did not raise). The exit code is 1 when any ASV raised or was demoted.
  * the new manifest _meta keeps the source _meta (selection, selection_stamp, gene_data,
    builder, source_digests) and adds 'pfunc' (estimator, parameters, gating rule, code sha256
    + git state, catalogue/store provenance, input_checks, counts) and 'reweighted_from'; the
    source build's pfunc/source_digests.json is copied along.
Inputs: the selection audit (default: the source manifest's _meta.selection), streamed with
ijson, only the in-scope records kept; the source genome dicts (--genome-cache-dir) for the
feature -> PGFam product-string route; the PGFam sets of every neighbour / clade member,
stream-extracted ONCE from the shared cache (never loaded whole, never written) into
OUT/pfunc/pgfam_store.json, missing genomes fetched from BV-BRC into OUT/pfunc/pgfam_fetched.json
(--pfunc-no-fetch to skip). Not copied: kbase_upload_log.json (it records uploads of the
source objects), mapping_provenance.json and the abstention reports (re-run
write_mapping_provenance.py on the new directory).

    ~/Documents/py_venv/bin/python reweight_build_e3.py \\
        --build-dir BUILD_KN \\
        --out-dir   BUILD_KN_e3 \\
        [--selection AUDIT] [--limit 300] [--asvs LIST.txt] [--procs 8]
        [--allow-audit-newer] [--allow-cache-newer]      # legacy builds only, see above

(BUILD_KN: a --pfunc kn build made by this code; the shipped d100 builds are refused, above.)

Memory (measured 2026-11-12 on the mixed-annotation EmilyKin d100 build, before the BV-BRC-only
rule: 1,410 n >= 2 ASVs, 20,697 of 21,003
genomes in the store, store file 0.69 GB): parent RSS 3.1 GB, each worker's PSS 0.6-0.8 GB,
354 s at --procs 8, no fetching (before the per-ASV output parse was added: ~94 ms per built
ASV, measured serially, i.e. ~45 s more at --procs 8 for EmilyKin's 3,755 built genomes). The
cgroup memory.peak reaches a 40G MemoryMax through page cache (the 5 GB cache stream + ~21 GB
written), not process memory. Run full builds under systemd-run. With --procs > 1 the workers are fork-started and share the store copy-on-write.

Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations

import argparse
import json
import multiprocessing
import os
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import build_synthetic_genomes as BSG  # noqa: E402
import pfunc_estimator as P  # noqa: E402

MANIFEST = "synthetic_genome_manifest.json"
SUBDIRS = ("genome_objects", "provenance", "taxonomies", "ASVset_taxonomies")


# --------------------------------------------------------------------------- inputs
def stream_records(path, wanted, log=print):
    """{asv: record} for the ASVs in `wanted`, plus the audit's _meta, streaming the audit
    with ijson (never json.load: a d100 audit is ~200 MB) and stopping once all are found."""
    import ijson
    wanted = set(wanted)
    out, meta = {}, {}
    t0 = time.time()
    with open(path, "rb") as fh:
        for k, v in ijson.kvitems(fh, "", use_float=True):
            if k == "_meta":
                meta = v
            elif k in wanted:
                out[k] = v
                if len(out) == len(wanted):
                    break
    if log:
        log(f"[reweight] selection: {len(out)}/{len(wanted)} records streamed from {path} "
            f"({time.time() - t0:.1f}s)")
    return out, meta


def load_sources(gids, cache_dir):
    """The built source genome dicts (fast-builder cache files), in order, BV-BRC PATRIC
    annotation only (build_synthetic_genomes.load_source_genome: the same filter the builders
    apply on read). Returns (sources, missing_ids); never fetches (the shared genome cache is
    not written)."""
    sources, missing = [], []
    for g in gids:
        p = os.path.join(cache_dir, f"{g}.json")
        if not os.path.exists(p):
            missing.append(g)
            continue
        sources.append(BSG.load_source_genome(p))
    return sources, missing


def check_annotation(src_meta):
    """Run-level: the source build must be BV-BRC PATRIC-only (manifest _meta.annotation ==
    build_synthetic_genomes.ANNOTATION_STAMP). Every build made before that stamp existed (git
    3ce12d0 and earlier, e.g. the shipped synthetic_genomes_d100 deployments) merged PATRIC and
    the imported RefSeq CDSs; reweighting one would carry its RefSeq-derived features into the
    new build, so it is refused (SystemExit, nothing written)."""
    got = src_meta.get("annotation")
    if got == BSG.ANNOTATION_STAMP:
        return
    raise SystemExit(
        f"[reweight] refusing: the source build is not BV-BRC PATRIC-only (manifest "
        f"_meta.annotation = {got!r}, expected {BSG.ANNOTATION_STAMP!r}). Builds made at or "
        f"before git {BSG.LAST_MIXED_COMMIT} merged PATRIC and RefSeq features; prFBA now uses "
        f"BV-BRC (PATRIC) annotations only. Rebuild it with the current builder (e.g. "
        f"deployment_d100/d100_pipeline.sh into a NEW SYN, or build_synthetic_genomes_parallel.py "
        f"--selection AUDIT --out-dir NEW); the rebuild applies E3 by default. Nothing was "
        f"written.")


# --------------------------------------------------------------------------- per ASV
FILE_KINDS = ("root", "genome_objects", "provenance", "taxonomies")
REQUIRED_KINDS = ("root", "genome_objects", "provenance")   # what 'built' promises downstream


def asv_paths(asv, base_dir):
    """{kind: path} of one ASV's four per-ASV files under a build directory."""
    return {k: (os.path.join(base_dir, k, f"{asv}.json") if k != "root"
                else os.path.join(base_dir, f"{asv}.json")) for k in FILE_KINDS}


def copy_asv(asv, build_dir, out_dir):
    """Copy one ASV's files byte for byte (the k/n build kept). A file the source lacks is
    removed from out_dir too (no half-reweighted leftover). Returns the kinds copied."""
    src, dst = asv_paths(asv, build_dir), asv_paths(asv, out_dir)
    done = []
    for k in FILE_KINDS:
        if os.path.exists(src[k]):
            shutil.copyfile(src[k], dst[k])
            done.append(k)
        elif os.path.exists(dst[k]):
            os.remove(dst[k])
    return done


def unusable_files(asv, out_dir, parse=False):
    """The REQUIRED_KINDS files of one ASV that are missing from out_dir (parse=True: or do
    not parse as JSON). A 'built' manifest entry promises these files."""
    bad = []
    for k in REQUIRED_KINDS:
        p = asv_paths(asv, out_dir)[k]
        if not os.path.exists(p):
            bad.append(k)
        elif parse:
            try:
                BSG._load(p)
            except Exception:
                bad.append(k)
    return bad


def drift_reason(which):
    return f"{BSG.REASON_DRIFT}: {which}"


# the <which> of 'inputs drifted since build: <which>' (fixed strings: they are counted as is)
DRIFT_AUDIT = "selection audit"
DRIFT_SELECTED = "selected genomes"
DRIFT_PREVALENCE = "function_prevalence"
DRIFT_DIGEST = "source digest"
DRIFT_CACHE_NEWER = "source genome cache newer than build"


def check_selection_audit(src_meta, sel_path, manifest_mtime_ns, allow_newer=False, log=print):
    """Run-level check that `sel_path` is the selection audit the build used (see the module
    docstring). Returns (record for _meta.pfunc.input_checks, forced reason | None): a build
    with _meta.selection_stamp whose sha256 differs from the audit's now -> the drift reason
    every built ASV then gets (k/n copy). A legacy build (no stamp) whose audit is newer than
    the source manifest -> SystemExit, unless allow_newer (warning; recorded)."""
    now = BSG.file_stamp(sel_path)
    stamp = src_meta.get("selection_stamp")
    rec = {"audit": sel_path, "audit_now": now}
    if isinstance(stamp, dict) and stamp.get("sha256"):
        rec.update(method="sha256 vs the build's manifest _meta.selection_stamp",
                   stamp_in_build=stamp)
        if now["sha256"] == stamp["sha256"]:
            rec["result"] = "match"
            return rec, None
        rec["result"] = "MISMATCH: every built ASV kept k/n"
        if log:
            log(f"[reweight] WARNING: the selection audit {sel_path} is not the one the build "
                f"used (sha256 {now['sha256'][:12]}... vs the build's stamp "
                f"{stamp['sha256'][:12]}...): every built ASV is copied unchanged (k/n) with "
                f"reason '{drift_reason(DRIFT_AUDIT)}'. Pass --selection with the audit the "
                f"build used.")
        return rec, drift_reason(DRIFT_AUDIT)
    rec.update(method="legacy build (no _meta.selection_stamp): audit mtime vs the source "
                      "manifest's mtime", manifest_mtime_ns=manifest_mtime_ns,
               audit_newer_than_build=now["mtime_ns"] > manifest_mtime_ns,
               allow_audit_newer=bool(allow_newer))
    if not rec["audit_newer_than_build"]:
        rec["result"] = "audit not newer than the build"
        return rec, None
    fmt = lambda ns: time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(ns / 1e9))  # noqa: E731
    msg = (f"the build predates the selection-audit stamp, and its audit {sel_path} "
           f"(mtime {fmt(now['mtime_ns'])}) is NEWER than the build's manifest "
           f"(mtime {fmt(manifest_mtime_ns)}): it may have been regenerated since the build, "
           f"and E3 would then use a neighbour list the build never saw")
    if not allow_newer:
        raise SystemExit(f"[reweight] refusing: {msg}. Nothing was written. If the audit is "
                         f"the one the build used (e.g. an unmodified copy), re-run with "
                         f"--allow-audit-newer; otherwise pass --selection with that audit.")
    rec["result"] = "audit newer than the build; allowed by --allow-audit-newer"
    if log:
        log(f"[reweight] WARNING: {msg}; continuing because --allow-audit-newer was given "
            f"(recorded in _meta.pfunc.input_checks.selection_audit)")
    return rec, None


def source_genome_check(build_dir, src_meta, manifest_mtime_ns, allow_cache_newer=False):
    """How check_inputs compares the cached source genomes with the build (see the module
    docstring): {'mode': 'digests', 'digests': {gid: digest}} when the build wrote
    pfunc/source_digests.json, else {'mode': 'mtime', 'manifest_mtime_ns',
    'allow_cache_newer'} (legacy build)."""
    p = os.path.join(build_dir, BSG.PFUNC_DIR, BSG.SOURCE_DIGESTS)
    if os.path.isfile(p):
        return {"mode": "digests", "file": p, "digests": BSG._load(p)}
    return {"mode": "mtime", "manifest_mtime_ns": manifest_mtime_ns,
            "allow_cache_newer": bool(allow_cache_newer),
            "file_named_in_meta_but_missing": bool(src_meta.get("source_digests"))}


def check_inputs(entry, rec, prov, cache_dir, src_check=None):
    """The retro tool's per-ASV input checks for one n >= min_n built ASV (see the module
    docstring, in that order). Returns (reason | None, sources, info): reason None when E3
    may be applied, with `sources` the built source genome dicts (in provenance order) loaded
    from cache_dir; info: {'digests': {gid: digest}} when computed, 'drift_genomes' (the ids
    that failed the digest / cache-mtime check), 'cache_newer_allowed_genomes' (legacy build,
    --allow-cache-newer)."""
    gids = [str(g) for g in (prov.get("genomes") or [])]
    info = {}
    if rec is None:
        return BSG.REASON_RECORD_MISSING, None, info
    if not rec.get("selected"):
        return "no selected reps", None, info
    sel_ids = [str(s["genome_id"]) for s in rec["selected"]]
    failed_ids = [str((x or {}).get("genome_id")) for x in entry.get("failed") or []]
    if sorted(sel_ids) != sorted(gids + failed_ids):
        return drift_reason(DRIFT_SELECTED), None, info
    sources, missing = load_sources(gids, cache_dir)
    if missing:
        return f"source genome(s) missing from the genome cache: {missing[:5]}", None, info
    if dict(BSG.function_prevalence(sources)) != (prov.get("function_prevalence") or {}):
        return drift_reason(DRIFT_PREVALENCE), None, info
    sc = src_check or {}
    if sc.get("mode") == "digests":
        now = P.source_digests(gids, sources)
        info["digests"] = now
        bad = [g for g in gids if sc["digests"].get(g) != now[g]]
        if bad:
            info["drift_genomes"] = bad
            return drift_reason(DRIFT_DIGEST), None, info
    elif sc.get("mode") == "mtime":
        newer = [g for g in gids if os.stat(os.path.join(cache_dir, f"{g}.json")).st_mtime_ns
                 > sc["manifest_mtime_ns"]]
        if newer and not sc.get("allow_cache_newer"):
            info["drift_genomes"] = newer
            return drift_reason(DRIFT_CACHE_NEWER), None, info
        if newer:
            info["cache_newer_allowed_genomes"] = newer
    return None, sources, info


def reweight_one(asv, entry, rec, build_dir, out_dir, cache_dir, ctx, forced_reason=None,
                 src_check=None):
    """Copy or reweight one built ASV into out_dir. Returns the manifest 'pfunc' summary.

    forced_reason: the run-level selection-audit drift (check_selection_audit): the ASV is
    copied unchanged (k/n) with that reason, whatever its n. Otherwise n < min_n (or --pfunc
    kn) -> copied unchanged; n >= min_n -> E3 only when check_inputs passes, i.e. exactly:
      * selected genomes: the audit record's selected ids == the provenance 'genomes' (built)
        plus the manifest entry's failed[] ids (as multisets). This does NOT compare the
        record's decisions[] (the neighbour list): a changed audit is caught only by the
        run-level sha256 stamp (or, for a legacy build, the audit-newer-than-build refusal);
      * the built source genomes are all in --genome-cache-dir;
      * function_prevalence of those cached dicts == the provenance's (re-fetched product
        strings, which the product-match route reads);
      * source digest: with the build's pfunc/source_digests.json, source_digest of every
        cached dict == the recorded one (also catches a changed PGFAM alias alone); a legacy
        build without it: no cached source file newer than the source manifest, unless
        --allow-cache-newer.
    A failed check -> copied unchanged with reason 'inputs drifted since build: <which>' (or
    'record not in selection audit' / 'source genome(s) missing from the genome cache')."""
    paths, dst = asv_paths(asv, build_dir), asv_paths(asv, out_dir)

    def copy_all():
        copy_asv(asv, build_dir, out_dir)

    n = int(entry.get("n_genomes") or 0)
    base = {"estimator": "E3" if ctx.e3 else "kn", "mode": ctx.mode, "n_built": n,
            "min_n": ctx.min_n}
    if forced_reason:
        copy_all()
        return {**base, "applied": False, "reason": forced_reason}
    with open(paths["provenance"]) as f:
        prov = json.load(f)
    gids = [str(g) for g in (prov.get("genomes") or [])]
    if len(gids) != n or [str(g) for g in entry.get("genomes") or []] != gids:
        copy_all()
        return {**base, "applied": False,
                "reason": f"n mismatch: manifest n_genomes {n}, provenance genomes {len(gids)}"}
    if not ctx.e3 or n < ctx.min_n:
        diag, _ = ctx.reweight(None, rec, None, n)       # returns before touching the genome
        copy_all()
        return BSG.pfunc_summary(diag)
    reason, sources, info = check_inputs(entry, rec, prov, cache_dir, src_check)
    if reason:
        copy_all()
        out = {**base, "applied": False, "reason": reason}
        if info.get("drift_genomes"):
            out["drift_genomes"] = info["drift_genomes"]
        return out
    extra = ({"cache_newer_allowed_genomes": info["cache_newer_allowed_genomes"]}
             if info.get("cache_newer_allowed_genomes") else {})
    root = BSG._load(paths["root"])
    diag, kn = ctx.reweight(root, rec, sources, n, built_ids=gids,
                            source_digests=info.get("digests"))
    if not diag.get("applied"):
        copy_all()
        return {**BSG.pfunc_summary(diag), **extra}
    q = BSG.reweighted_q(root, kn)
    go = BSG._load(paths["genome_objects"])
    old_go = P.apply_to_genome(go, q)
    if old_go != kn:
        copy_all()
        return {**base, "applied": False,
                "reason": "{asv}.json and genome_objects/{asv}.json carry different k/n"}
    BSG._dump_json(root, dst["root"])
    BSG._dump_json(go, dst["genome_objects"])
    prov["probability_kn"] = kn
    prov["pfunc"] = diag
    BSG._dump_json(prov, dst["provenance"])
    if os.path.exists(paths["taxonomies"]):
        shutil.copyfile(paths["taxonomies"], dst["taxonomies"])
    return {**BSG.pfunc_summary(diag), **extra}


_W = {}


def _winit(ctx, build_dir, out_dir, cache_dir, forced_reason=None, src_check=None):
    _W.update(ctx=ctx, build_dir=build_dir, out_dir=out_dir, cache_dir=cache_dir,
              forced_reason=forced_reason, src_check=src_check)


def _wone(task):
    """reweight_one, never raising: an ASV that raises is copied unchanged (the k/n build
    kept, overwriting anything half written) and reported with reason 'error: ...'. Then, in
    either case, its required output files are parsed. Returns (asv, summary | None,
    error | None, required files missing or unparseable in out_dir)."""
    asv, entry, rec = task
    try:
        summ = reweight_one(asv, entry, rec, _W["build_dir"], _W["out_dir"], _W["cache_dir"],
                            _W["ctx"], forced_reason=_W.get("forced_reason"),
                            src_check=_W.get("src_check"))
    except Exception as exc:          # one ASV never kills the run; it keeps k/n
        err = f"{type(exc).__name__}: {str(exc)[:300]}"
        try:
            copy_asv(asv, _W["build_dir"], _W["out_dir"])
        except Exception as exc2:
            err += f"; fallback copy failed: {type(exc2).__name__}: {str(exc2)[:200]}"
        return asv, None, err, unusable_files(asv, _W["out_dir"], parse=True)
    # a copied file can be unreadable without anything raising (a corrupt source build)
    return asv, summ, None, unusable_files(asv, _W["out_dir"], parse=True)


# --------------------------------------------------------------------------- main
def _check_dirs(build_dir, out_dir):
    b, o = os.path.realpath(build_dir), os.path.realpath(out_dir)
    if b == o or o.startswith(b + os.sep) or b.startswith(o + os.sep):
        raise SystemExit(f"[reweight] --out-dir must be a new directory outside --build-dir "
                         f"({out_dir} vs {build_dir})")
    if os.path.exists(o) and os.listdir(o):
        raise SystemExit(f"[reweight] --out-dir {out_dir} exists and is not empty; give a NEW "
                         f"directory (copy-on-write: the source build is never modified)")
    if not os.path.exists(os.path.join(b, MANIFEST)):
        raise SystemExit(f"[reweight] no {MANIFEST} in --build-dir {build_dir}")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--build-dir", required=True, help="existing build (manifest, {asv}.json, "
                    "genome_objects/, provenance/, taxonomies/); only read")
    ap.add_argument("--out-dir", required=True, help="NEW directory for the reweighted build")
    ap.add_argument("--selection", default=None,
                    help="selection audit (default: the source manifest's _meta.selection)")
    ap.add_argument("--genome-cache-dir", default=os.path.join(HERE, "bvbrc_cache", "kbase_genome_cache"),
                    help="the source genome dicts the build merged (fast-builder cache)")
    ap.add_argument("--limit", type=int, default=0,
                    help="only the first N BUILT ASVs of the manifest (abstain entries are then "
                         "left out); 0 = all entries")
    ap.add_argument("--asvs", default=None,
                    help="file with one ASV id per line: only these (built) ASVs, in manifest "
                         "order (overrides --limit)")
    ap.add_argument("--procs", type=int, default=1, help="worker processes (fork; default 1)")
    ap.add_argument("--allow-audit-newer", action="store_true",
                    help="legacy builds (no _meta.selection_stamp) only: proceed although the "
                         "selection audit's mtime is newer than the source manifest's (else the "
                         "run stops); recorded in _meta.pfunc.input_checks")
    ap.add_argument("--allow-cache-newer", action="store_true",
                    help="legacy builds (no pfunc/source_digests.json) only: apply E3 although "
                         "an ASV's source-genome cache file is newer than the source manifest "
                         "(else that ASV keeps k/n, 'inputs drifted since build: source genome "
                         "cache newer than build'); recorded per ASV and in _meta.pfunc")
    BSG.add_pfunc_args(ap, default="e3")
    args = ap.parse_args(argv)

    t0 = time.time()
    build_dir = os.path.abspath(args.build_dir)
    out_dir = os.path.abspath(args.out_dir)
    cache_dir = os.path.abspath(args.genome_cache_dir)
    _check_dirs(build_dir, out_dir)
    preflight = BSG.pfunc_preflight(args.pfunc, args.pfunc_catalogue, args.pfunc_pgfam_cache,
                                    allow_stale=args.pfunc_allow_stale_catalogue)
    manifest_path = os.path.join(build_dir, MANIFEST)
    manifest_mtime_ns = os.stat(manifest_path).st_mtime_ns
    with open(manifest_path) as f:
        src_manifest = json.load(f)
    src_meta = src_manifest.pop("_meta", {}) or {}
    if (src_meta.get("pfunc") or {}).get("estimator") == P.ESTIMATOR_NAME:
        raise SystemExit("[reweight] the source build already carries E3 probabilities "
                         "(manifest _meta.pfunc); reweight a k/n build")
    check_annotation(src_meta)
    sel_path = os.path.abspath(args.selection or src_meta.get("selection") or "")
    if not os.path.isfile(sel_path):
        raise SystemExit(f"[reweight] selection audit not found: {sel_path!r} (pass --selection)")
    if src_meta.get("selection") and os.path.abspath(src_meta["selection"]) != sel_path:
        print(f"[reweight] WARNING: --selection {sel_path} differs from the build's "
              f"{src_meta['selection']}", flush=True)

    only = None
    if args.asvs:
        with open(args.asvs) as f:
            only = {x.strip() for x in f if x.strip()}
    scope, n_built_scope = [], 0
    for a, e in src_manifest.items():
        built = e.get("status") == "built"
        if only is not None:
            if a in only and built:
                scope.append(a)
            continue
        if args.limit:
            if not built:
                continue
            if n_built_scope >= args.limit:
                break
        scope.append(a)
        n_built_scope += built
    built_scope = [a for a in scope if src_manifest[a].get("status") == "built"]
    n_by_asv = {a: int(src_manifest[a].get("n_genomes") or 0) for a in built_scope}
    built_by_asv = {a: [str(g) for g in src_manifest[a].get("genomes") or []] for a in built_scope}
    print(f"[reweight] {build_dir} -> {out_dir} | {len(scope)} entries ({len(built_scope)} built; "
          f"{sum(1 for a in built_scope if n_by_asv[a] >= args.pfunc_min_n)} with n >= "
          f"{args.pfunc_min_n}) | pfunc={args.pfunc}", flush=True)

    want = [a for a in built_scope if n_by_asv[a] >= args.pfunc_min_n] if args.pfunc == "e3" else []
    # run-level input checks, before anything is written (a legacy build whose audit is newer
    # than its manifest stops here unless --allow-audit-newer)
    input_checks, forced, src_check = {}, None, None
    if want:
        input_checks["selection_audit"], forced = check_selection_audit(
            src_meta, sel_path, manifest_mtime_ns, allow_newer=args.allow_audit_newer)
        src_check = source_genome_check(build_dir, src_meta, manifest_mtime_ns,
                                        allow_cache_newer=args.allow_cache_newer)
        sg = {k: v for k, v in src_check.items() if k != "digests"}
        if src_check["mode"] == "digests":
            sg["n_genomes_recorded"] = len(src_check["digests"])
        elif args.allow_cache_newer:
            print("[reweight] WARNING: legacy build (no pfunc/source_digests.json) and "
                  "--allow-cache-newer: source-genome cache files newer than the build are "
                  "accepted (listed per ASV in pfunc 'cache_newer_allowed_genomes')", flush=True)
        input_checks["source_genomes"] = sg
        if forced:
            want = []          # every built ASV is copied: no audit streaming, no PGFam store

    os.makedirs(out_dir, exist_ok=True)
    for d in SUBDIRS:
        os.makedirs(os.path.join(out_dir, d), exist_ok=True)
    recs, sel_meta = stream_records(sel_path, want) if want else ({}, {})
    ctx, pfunc_meta = BSG.prepare_pfunc(
        recs, out_dir, mode=args.pfunc, min_n=args.pfunc_min_n, dedup=args.pfunc_dedup,
        pgfam_cache=os.path.abspath(args.pfunc_pgfam_cache),
        catalogue=os.path.abspath(args.pfunc_catalogue), fetch=not args.pfunc_no_fetch,
        fetch_workers=args.pfunc_fetch_workers, n_by_asv=n_by_asv, built_by_asv=built_by_asv,
        allow_stale=args.pfunc_allow_stale_catalogue, preflight=preflight,
        code_files=BSG.PFUNC_PROVENANCE_FILES + [os.path.abspath(__file__)])

    tasks = [(a, src_manifest[a], recs.get(a)) for a in built_scope]
    results, errors, unusable = {}, {}, {}
    t1 = time.time()
    if args.procs > 1 and len(tasks) > 1:
        with ProcessPoolExecutor(max_workers=args.procs, initializer=_winit,
                                 initargs=(ctx, build_dir, out_dir, cache_dir, forced,
                                           src_check),
                                 mp_context=multiprocessing.get_context("fork")) as ex:
            it = ex.map(_wone, tasks, chunksize=2)
            for i, (a, summ, err, bad) in enumerate(it, 1):
                results[a], errors[a], unusable[a] = summ, err, bad
                if i % 100 == 0 or i == len(tasks):
                    print(f"[reweight] {i}/{len(tasks)} ({time.time() - t1:.0f}s)", flush=True)
    else:
        _winit(ctx, build_dir, out_dir, cache_dir, forced, src_check)
        for i, t in enumerate(tasks, 1):
            a, summ, err, bad = _wone(t)
            results[a], errors[a], unusable[a] = summ, err, bad
            if i % 100 == 0 or i == len(tasks):
                print(f"[reweight] {i}/{len(tasks)} ({time.time() - t1:.0f}s)", flush=True)

    # manifest: a 'built' entry always has its files (REQUIRED_KINDS) in out_dir, parseable
    # (checked in the worker after every ASV). An ASV that raised was copied unchanged (k/n
    # kept, reason 'error: ...'); if a required file is still missing or unparseable (a broken
    # source build, with or without an exception), the entry becomes 'failed'.
    manifest, n_err, n_err_failed, n_unreadable = {}, 0, 0, 0
    n_cache_newer_allowed = 0
    for a in scope:
        e = dict(src_manifest[a])
        if a in results:
            if errors.get(a):
                n_err += 1
                e["pfunc"] = {"estimator": "E3" if ctx.e3 else "kn", "mode": ctx.mode,
                              "n_built": int(e.get("n_genomes") or 0), "min_n": ctx.min_n,
                              "applied": False, "reason": f"error: {errors[a]}"}
            elif ctx.e3:
                e["pfunc"] = results[a]
                n_cache_newer_allowed += bool((results[a] or {}).get("cache_newer_allowed_genomes"))
            bad = unusable.get(a) or []
            if bad:
                n_err_failed += 1
                n_unreadable += not errors.get(a)
                e["status"] = "failed"
                e["error"] = (f"reweight_build_e3: {', '.join(bad)} missing or unreadable in "
                              f"{out_dir} "
                              + (f"after the k/n fallback copy ({errors[a]})" if errors.get(a)
                                 else "(copied from the source build, whose file is missing or "
                                      "unreadable)"))[:500]
        manifest[a] = e
    if n_err or n_err_failed:
        print(f"[reweight] {n_err} ASVs raised errors and were copied unchanged (k/n kept; see "
              f"manifest pfunc.reason); {n_err_failed} of the scope's built ASVs have missing or "
              f"unreadable files and are marked 'failed' ({n_unreadable} of them without an "
              f"error)", flush=True)
    if ctx.e3:
        pfunc_meta.update(BSG.pfunc_counts(manifest.values()))
    if input_checks:
        if src_check and src_check["mode"] == "mtime":
            input_checks["source_genomes"]["n_asvs_cache_newer_allowed"] = n_cache_newer_allowed
        pfunc_meta["input_checks"] = input_checks
    pfunc_meta["n_errors"] = n_err
    pfunc_meta["n_errors_kept_kn"] = n_err - sum(1 for a in scope if errors.get(a)
                                                  and manifest[a].get("status") == "failed")
    pfunc_meta["n_built_demoted_to_failed"] = n_err_failed
    pfunc_meta["n_built_demoted_unreadable_output"] = n_unreadable
    # the source genomes' digests describe this build too (copied ASVs used them; E3 was
    # applied only where the cached genomes still match them)
    dig_src = os.path.join(build_dir, BSG.PFUNC_DIR, BSG.SOURCE_DIGESTS)
    if os.path.isfile(dig_src):
        os.makedirs(os.path.join(out_dir, BSG.PFUNC_DIR), exist_ok=True)
        shutil.copyfile(dig_src, os.path.join(out_dir, BSG.PFUNC_DIR, BSG.SOURCE_DIGESTS))
    statuses = [manifest[a].get("status") for a in scope]
    meta = {**src_meta, "annotation": BSG.ANNOTATION_STAMP,
            "n_asvs": len(scope), "built": statuses.count("built"),
            "abstain": statuses.count("abstain"), "failed": statuses.count("failed"),
            "selection": sel_path, "pfunc": pfunc_meta,
            "reweighted_from": {"build_dir": build_dir, "source_meta": src_meta,
                                "scope": ("asv list " + os.path.abspath(args.asvs)) if args.asvs
                                else (f"first {args.limit} built ASVs" if args.limit else "all"),
                                "genome_cache_dir": cache_dir,
                                "selection_meta_gene_data": (sel_meta or {}).get("gene_data"),
                                "tool": "reweight_build_e3.py",
                                "seconds": round(time.time() - t0, 1)}}
    with open(os.path.join(out_dir, MANIFEST), "w") as f:
        json.dump({"_meta": meta, **manifest}, f, indent=1)
    c = pfunc_meta
    print(f"[reweight] done in {time.time() - t0:.0f}s: applied={c.get('n_applied')} "
          f"not_applied={c.get('not_applied_by_reason')} errors={n_err} "
          f"demoted_to_failed={n_err_failed} -> {os.path.join(out_dir, MANIFEST)}", flush=True)
    if forced:
        print(f"[reweight] WARNING: nothing was reweighted: {forced} (see "
              f"_meta.pfunc.input_checks.selection_audit)", flush=True)
    n_newer = (c.get("not_applied_by_reason") or {}).get(drift_reason(DRIFT_CACHE_NEWER), 0)
    if n_newer:
        print(f"[reweight] WARNING: {n_newer} ASVs kept k/n because a source-genome cache file is "
              f"newer than the build's manifest (legacy build without pfunc/source_digests.json). "
              f"If those files are unchanged copies (e.g. a copied cache), re-run with "
              f"--allow-cache-newer", flush=True)
    return 0 if not (n_err or n_err_failed) else 1


if __name__ == "__main__":
    sys.exit(main())
