#!/usr/bin/env python
"""
upload_kbase_genomes.py — save the per-ASV synthetic genomes of a build directory to a KBase
narrative workspace as KBaseGenomes.Genome objects (+ one KBaseSearch.GenomeSet grouping them),
carrying the full mapping provenance in each object's metadata + workspace provenance record.

  * resumable: objects already present in the workspace (by name) are skipped
  * batched save_objects (by cumulative size), bisecting a failing batch to isolate bad objects
  * each Genome is made typespec-valid for KBaseGenomes.Genome-11.1 without losing information:
      - empty assembly_ref / taxon_ref are dropped (they are optional; '' is not a valid ws ref)
      - legacy ontology_terms {src: [term,...]} -> {src: {term: [event_idx]}} + ontology_events
      - per-feature prFBA fields (probability, selection_provenance, source_genome_count) are
        ALSO written into the spec-sanctioned inference_data list so they survive any consumer
      - notes/taxon_assignments summarise the mapping provenance
  * per-object workspace metadata: asv, tier, best identity, anchor + every source genome,
    flags, prFBA commit; provenance action: prFBA service/commit, selection knobs, external_data
    = the BV-BRC source genomes
  * writes kbase_upload_log.json (asv -> UPA) into the genome dir
  * E3 builds (manifest _meta.pfunc.estimator == 'E3') also record the estimator: provenance
    method_params[0]['feature_probability'] = the manifest's _meta.pfunc stamp, and ws
    metadata 'probability_estimator' ('E3', or 'k/n (E3 not applied: <reason>)' for that
    genome); a k/n build (or one made before E3) uploads exactly as before. Objects are named
    by ASV and an existing name is skipped (resume), so uploading an E3 rebuild of a build
    already in a workspace needs a NEW workspace (--ws-id) or new object names.
  * functions: a BV-BRC PATRIC-only build (manifest _meta.annotation, every build since
    2026-11-13) is described as BV-BRC (PATRIC) annotations in the RAST ontology event, the notes,
    the provenance description + method_params[0]['function_annotation'], ws metadata
    'function_annotation' and the GenomeSet description.
  * REFUSES any other build (since 2026-11-13; only BV-BRC's own annotation is used): a manifest
    without that stamp (every build made at or before prFBA 3ce12d0, which merged PATRIC and
    RefSeq features) exits non-zero before the token is read or KBase is contacted, in every
    mode including --dry-run, so nothing is uploaded and no kbase_upload_log.json is written.
    There is no override; rebuild it. (The mixed-build description text is kept in the helpers
    below but is no longer reachable from main.) The stamp certifies the query, not that every
    source genome was fetched complete: see bvbrc_cache/CORRECTIONS_2026-11-13.json.
  * FEATURE GATE (since 2026-11-13, after the stamp gate, also before any token/network use and
    in every mode): every built genome's {asv}.json (the file uploaded) is streamed, one at a
    time (--gate-workers processes), and the build is REFUSED (exit 1, nothing uploaded or
    logged) if any feature or CDS (also non_coding_features / mrnas) has no PATRIC_id alias, an
    empty one, or one that does not start with 'fig|<g>.' for a g in that ASV's source genomes
    (provenance/{asv}.json 'genomes', which must equal the manifest entry's 'genomes'). The
    stamp records which query a build used, not where each feature came from: a stamped build
    from a cache file holding RefSeq rows with ANOTHER genome's patric_id (294.202's RefSeq CDSs
    carry fig|294.71.peg.N) passed the old non-empty filter and the stamp gate alike. There is
    no override; rebuild the ASVs it names. The gate is per ASV (a merged feature does not record
    which source it came from, so an id of ANY of the ASV's sources passes); the per-genome rule
    is build_synthetic_genomes.bvbrc_only's. On the 2026-11-13 d100 BV-BRC-only builds it took
    28-37 s (8 processes; 3,225 / 3,755 genomes) and refused EmilyKin's on 7ac9047b alone.

    ~/Documents/py_venv/bin/python upload_kbase_genomes.py --genome-dir DIR --ws-id 261948 [--limit 3]
"""
from __future__ import annotations

import argparse
import concurrent.futures
import datetime as dt
import json
import os
import sys
import time
import copy
from collections import Counter

import requests

import pfunc_estimator as P   # probability_description: what Feature.probability means

WS_URL = "https://kbase.us/services/ws"
TOKEN_FILE = os.path.expanduser("~/.kbase/token")
GENOME_TYPE = "KBaseGenomes.Genome"
GENOMESET_TYPE = "KBaseSearch.GenomeSet"
PRFBA_FEATURE_KEYS = ("probability", "selection_provenance", "source_genome_count")
GATE_COLLECTIONS = ("features", "cdss", "non_coding_features", "mrnas")   # every one checked
GATE_EXAMPLES = 3            # example offending PATRIC_ids kept per ASV for the message

try:                         # the feature gate parses every built genome once: orjson if present
    import orjson as _orjson

    def _json_bytes(b):
        return _orjson.loads(b)
except ImportError:          # pragma: no cover - the venv has orjson
    def _json_bytes(b):
        return json.loads(b)


# --------------------------------------------------------------------------- ws client
class WS:
    def __init__(self, token, url=WS_URL):
        self.url, self.h = url, {"Authorization": token}

    def call(self, method, params, timeout=3600, retries=3):
        last = None
        for i in range(retries):
            try:
                r = requests.post(self.url, json={"method": method, "params": params, "version": "1.1", "id": "1"},
                                  headers=self.h, timeout=timeout)
                j = r.json()
                if "error" in j:
                    raise RuntimeError(j["error"].get("message") or j["error"].get("error") or str(j["error"])[:2000])
                return j["result"][0] if j.get("result") else None
            except (requests.ConnectionError, requests.Timeout, ValueError) as exc:  # transient
                last = exc
                time.sleep(5 * (i + 1))
        raise RuntimeError(f"transient failure after {retries} tries: {last}")

    def existing(self, ws_id, type_prefix=None):
        out = {}
        params = {"ids": [ws_id], "limit": 10000}
        if type_prefix:
            params["type"] = type_prefix
        for o in self.call("Workspace.list_objects", [params]) or []:
            out[o[1]] = f"{o[6]}/{o[0]}/{o[4]}"
        return out


def now_kbase():
    """KBase provenance timestamp: YYYY-MM-DDThh:mm:ss+hhmm."""
    return dt.datetime.now().astimezone().strftime("%Y-%m-%dT%H:%M:%S%z")


# --------------------------------------------------------------------------- genome transform
def _fix_ontology(feat, event_idx):
    ot = feat.get("ontology_terms")
    if not isinstance(ot, dict):
        return False
    changed = False
    for src, terms in list(ot.items()):
        if isinstance(terms, list):
            ot[src] = {str(t): [event_idx] for t in terms}
            changed = True
        elif isinstance(terms, dict):
            for t, v in list(terms.items()):
                if not isinstance(v, list):
                    terms[t] = [event_idx]; changed = True
    return changed


def annotation_meta(annotation):
    """The manifest _meta subset pfunc_estimator.functions_description reads."""
    return {"annotation": annotation} if annotation else {}


def to_kbase_genome(g, *, asv, sel, man, commit, keep_extra_keys=True, pfunc_meta=None,
                    annotation=None):
    """pfunc_meta: the build manifest's _meta['pfunc'] (None for builds made before E3, or
    estimator 'kn': the probability text is then the original k/n sentence, unchanged).
    annotation: the manifest's _meta['annotation'] (BV-BRC PATRIC-only builds); None for a
    build made before it, whose functions mix PATRIC and RefSeq and are described so."""
    g = copy.deepcopy(g)
    for k in ("assembly_ref", "taxon_ref"):
        if k in g and not g[k]:
            g.pop(k)
    ts = now_kbase()
    fdesc = P.functions_description(annotation_meta(annotation))
    if not g.get("ontology_events"):
        if P.is_bvbrc_only_build(annotation_meta(annotation)):
            method = ("prFBA synthetic genome: union of the BV-BRC (PATRIC) annotation functions of the "
                      "selected 16S-matched reference genomes")
            desc = ("term keys are the BV-BRC (PATRIC, RASTtk) function strings carried in "
                    "feature.functions; prFBA " + commit)
        else:
            method = ("prFBA synthetic genome: union of the BV-BRC genome_feature CDS products (PATRIC "
                      "and imported RefSeq annotations) of the selected 16S-matched reference genomes")
            desc = ("term keys are the function strings carried in feature.functions (BV-BRC PATRIC "
                    "and RefSeq products; build made before the BV-BRC-only rule); prFBA " + commit)
        g["ontology_events"] = [{"id": "RAST", "method": method,
                                 "method_version": commit, "timestamp": ts,
                                 "description": desc}]
    ev = len(g["ontology_events"]) - 1
    for coll in ("features", "cdss", "non_coding_features", "mrnas"):
        for f in g.get(coll) or []:
            _fix_ontology(f, ev)
            if coll == "features":
                inf = [x for x in (f.get("inference_data") or []) if x.get("category") != "prFBA"]
                for k in PRFBA_FEATURE_KEYS:
                    if k in f:
                        inf.append({"category": "prFBA", "type": k, "evidence": str(f[k])})
                if inf:
                    f["inference_data"] = inf
            if not keep_extra_keys:
                for k in PRFBA_FEATURE_KEYS:
                    f.pop(k, None)
    selected = sel.get("selected") or []
    anchor = selected[0] if selected else {}
    srcs = ", ".join(f"{s.get('genome_id')} ({s.get('organism')}, id={s.get('identity')}, {s.get('role')})" for s in selected)
    g["notes"] = (f"prFBA synthetic genome for 16S ASV {asv} (codiffusion bioreactor, V4-V5). Built by prFBA commit "
                  f"{commit} from the exact-tie de-duplicated union (+P3 tie-break) of BV-BRC 16S-matched reference "
                  f"genomes: {srcs}. Tier={sel.get('tier')}, best identity={sel.get('best_identity')}, "
                  f"flags={sel.get('flags')}. Functions: {fdesc}. Feature.probability = "
                  f"{P.probability_description(pfunc_meta, (man or {}).get('pfunc'))} "
                  f"(also in inference_data category=prFBA). No sequences: functions + RAST ontology + PATRIC/PGFam aliases only.")
    if anchor.get("taxon_id") is not None:
        g["taxon_assignments"] = {"ncbi": str(anchor["taxon_id"])}
    g.setdefault("warnings", [])
    return g


def is_e3_build(pfunc_meta):
    """True for a build whose manifest _meta.pfunc names the E3 estimator. Only those builds
    get the extra 'feature_probability' / 'probability_estimator' provenance, so a k/n build
    (or one made before E3) uploads exactly as before."""
    return bool(pfunc_meta) and pfunc_meta.get("estimator") == P.ESTIMATOR_NAME


def probability_estimator_tag(entry_pfunc):
    """Short ws-metadata value: which estimator set THIS genome's feature probabilities."""
    pf = entry_pfunc or {}
    if pf.get("applied"):
        return "E3"
    return f"k/n (E3 not applied: {pf.get('reason')})"[:200]


def ws_meta(asv, sel, man, commit, rel_ab, pfunc_meta=None, annotation=None):
    selected = sel.get("selected") or []
    anchor = selected[0] if selected else {}
    m = {
        "asv": asv, "prFBA_commit": commit[:12],
        "selection": "exact-tie dedup union + P3 gene tie-break (select_references.DEFAULT_KNOBS)",
        "tier": str(sel.get("tier")), "best_identity": str(sel.get("best_identity")),
        "asv_len": str(sel.get("asv_len")), "rel_abundance": str(rel_ab),
        "n_source_genomes": str(len(selected)),
        "source_genomes": "|".join(str(s.get("genome_id")) for s in selected)[:800],
        "anchor_genome": str(anchor.get("genome_id")), "anchor_organism": str(anchor.get("organism"))[:200],
        "anchor_identity": str(anchor.get("identity")), "anchor_ncbi_taxon_id": str(anchor.get("taxon_id")),
        "consensus_family": str(sel.get("consensus_family")),
        "flags": ",".join(sel.get("flags") or [])[:800],
        "n_features": str(man.get("n_features")),
        "n_equivalent_genomes": str(man.get("n_equivalent_genomes")),
        "midas_taxonomy": (json.dumps(sel.get("midas_taxonomy")) if not isinstance(sel.get("midas_taxonomy"), str)
                           else sel.get("midas_taxonomy"))[:800],
    }
    if is_e3_build(pfunc_meta):           # E3 builds only: k/n uploads keep their metadata
        m["probability_estimator"] = probability_estimator_tag((man or {}).get("pfunc"))
    if P.is_bvbrc_only_build(annotation_meta(annotation)):   # stamped builds only
        m["function_annotation"] = "BV-BRC PATRIC (RASTtk) only"
    return {k: v for k, v in m.items() if v not in (None, "None")}


def ws_provenance(asv, sel, prov, commit, hits_path, sel_path, gdir, pfunc_meta=None,
                  annotation=None):
    """pfunc_meta: the manifest _meta['pfunc']; for an E3 build it is recorded as
    method_params[0]['feature_probability'] (estimator, parameters, gating rule, code sha256
    and git state, counts). A k/n build's provenance is unchanged. annotation: the manifest
    _meta['annotation']; a BV-BRC PATRIC-only build records it as
    method_params[0]['function_annotation'] and says so in the description."""
    selected = sel.get("selected") or []
    knobs = ((prov.get("stages") or {}).get("3_selection") or {}).get("knobs") or {}
    ext = [{"resource_name": "BV-BRC", "resource_url": "https://www.bv-brc.org",
            "data_id": str(s.get("genome_id")), "data_url": f"https://www.bv-brc.org/view/Genome/{s.get('genome_id')}",
            "description": f"source reference genome: {s.get('organism')}; 16S identity={s.get('identity')}, "
                           f"coverage={s.get('coverage')}, role={s.get('role')}, confidence={s.get('confidence')}"}
           for s in selected]
    gitinfo = (prov.get("pipeline") or {}).get("git") or {}
    out = [{
        "time": now_kbase(), "service": "prFBA", "service_ver": commit,
        "method": "gpu_align.py --full -> select_references (exact-tie dedup union + P3) -> build_synthetic_genomes_parallel.py",
        "method_params": [{"asv": asv, "selection_knobs": knobs, "p3_gene_tiebreak": True, "gene_data": "bvbrc",
                           "alignment_scoring": ((prov.get("stages") or {}).get("1_alignment") or {}).get("scoring"),
                           "tier": sel.get("tier"), "best_identity": sel.get("best_identity"),
                           "flags": sel.get("flags"), "selected": [
                               {k: s.get(k) for k in ("genome_id", "organism", "taxon_id", "identity", "coverage",
                                                       "role", "confidence", "est_genes", "gene_set_empty")}
                               for s in selected]}],
        "script": "upload_kbase_genomes.py", "script_ver": commit, "script_command_line": " ".join(sys.argv),
        "description": f"prFBA 16S ASV -> reference-genome -> synthetic-genome mapping; commit {commit} "
                       f"({gitinfo.get('commit_date')}); repo {gitinfo.get('remote')}",
        "external_data": ext,
        "custom": {"alignment_hits": hits_path, "selection_audit": sel_path,
                   "mapping_provenance": os.path.join(gdir, "mapping_provenance.json"),
                   "per_asv_provenance": os.path.join(gdir, "provenance", f"{asv}.json"),
                   "genome_json": os.path.join(gdir, f"{asv}.json")},
    }]
    if is_e3_build(pfunc_meta):
        out[0]["method_params"][0]["feature_probability"] = pfunc_meta
    if P.is_bvbrc_only_build(annotation_meta(annotation)):
        out[0]["method_params"][0]["function_annotation"] = annotation
        out[0]["description"] += "; functions: " + P.FUNCTIONS_BVBRC
    return out


# --------------------------------------------------------------------------- upload
def save_batch(ws, ws_id, objs, log, keep_extra_keys, stats):
    """Save a batch; on failure bisect to isolate the offending object(s)."""
    if not objs:
        return
    t0 = time.perf_counter()
    try:
        infos = ws.call("Workspace.save_objects", [{"id": ws_id, "objects": [o["obj"] for o in objs]}])
        dt_ = time.perf_counter() - t0
        for o, info in zip(objs, infos):
            log[o["asv"]] = {"name": info[1], "upa": f"{info[6]}/{info[0]}/{info[4]}", "type": info[2],
                             "size_bytes": info[9], "saved_at": info[3]}
        stats["n"] += len(objs); stats["bytes"] += sum(o["nbytes"] for o in objs); stats["calls"] += 1
        mb = sum(o["nbytes"] for o in objs) / 1e6
        print(f"[upload] +{len(objs)} ({mb:.0f} MB in {dt_:.0f}s, {mb/max(dt_,1e-6):.1f} MB/s) "
              f"total={stats['n']}", flush=True)
    except Exception as exc:
        msg = str(exc)
        if len(objs) > 1:
            print(f"[upload] batch of {len(objs)} failed ({msg[:160]}); bisecting", flush=True)
            h = len(objs) // 2
            save_batch(ws, ws_id, objs[:h], log, keep_extra_keys, stats)
            save_batch(ws, ws_id, objs[h:], log, keep_extra_keys, stats)
        else:
            o = objs[0]
            if keep_extra_keys and ("dditional" in msg or "not allowed" in msg or "unexpected" in msg.lower()):
                print(f"[upload] {o['asv']}: extra feature keys rejected -> retrying with them stripped", flush=True)
                o2 = dict(o); o2["obj"] = dict(o["obj"]); o2["obj"]["data"] = to_kbase_genome(
                    o["raw"], asv=o["asv"], sel=o["sel"], man=o["man"], commit=o["commit"], keep_extra_keys=False,
                    pfunc_meta=o.get("pfunc_meta"), annotation=o.get("annotation"))
                stats["stripped"] += 1
                return save_batch(ws, ws_id, [o2], log, False, stats)
            log[o["asv"]] = {"error": msg[:2000]}
            stats["failed"] += 1
            print(f"[upload] FAILED {o['asv']}: {msg[:300]}", flush=True)


def require_bvbrc_only_build(man_meta, gdir):
    """Refuse (SystemExit with a message, exit status 1) a build whose manifest _meta lacks the
    BV-BRC PATRIC-only annotation stamp (pfunc_estimator.ANNOTATION_STAMP): every build made at
    or before prFBA 3ce12d0 merged BV-BRC's PATRIC and imported RefSeq features. main() calls it
    before reading the token or contacting KBase, so nothing is uploaded or logged."""
    if P.is_bvbrc_only_build(man_meta):
        return
    raise SystemExit(
        f"[upload] REFUSED: {gdir} is not a BV-BRC PATRIC-only build (manifest _meta.annotation = "
        f"{(man_meta or {}).get('annotation')!r}, expected {P.ANNOTATION_STAMP!r}). It was made at "
        f"or before prFBA 3ce12d0 and its genomes merge PATRIC and RefSeq features; prFBA uploads "
        f"BV-BRC (PATRIC) annotations only. Rebuild it (deployment_d100/README.md section 10). "
        f"Nothing was uploaded.")


def _patric_id_owner(pid):
    """'fig|294.71.peg.2281' -> 'fig|294.71.' (the genome a PATRIC id belongs to, for messages)."""
    if not isinstance(pid, str) or not pid.startswith("fig|"):
        return repr(pid)
    parts = pid[4:].split(".")
    return "fig|" + ".".join(parts[:2]) + "." if len(parts) >= 3 else pid


def patric_id_violations(genome, source_ids):
    """Check one built genome object against the BV-BRC-own rule: every feature (and CDS,
    non-coding feature, mRNA) carries at least one ['PATRIC_id', id] alias and every such id
    starts with 'fig|<g>.' (pfunc_estimator.patric_id_prefix) for a g in source_ids. Returns
    {n_checked, n_bad, missing, empty, foreign: {'fig|<owner>.': n}, examples: [...]}; n_bad = 0
    means it passes."""
    prefixes = tuple(P.patric_id_prefix(g) for g in source_ids)
    n_checked = n_bad = missing = empty = 0
    foreign, examples = Counter(), []
    for coll in GATE_COLLECTIONS:
        for f in genome.get(coll) or ():
            n_checked += 1
            pids = [a[1] for a in (f.get("aliases") or ())
                    if isinstance(a, (list, tuple)) and len(a) >= 2 and a[0] == "PATRIC_id"]
            why = []
            if not pids:
                missing += 1
                why.append("no PATRIC_id")
            for pid in pids:
                if not isinstance(pid, str) or not pid:
                    empty += 1
                    why.append(f"PATRIC_id {pid!r}")
                elif not (prefixes and pid.startswith(prefixes)):
                    foreign[_patric_id_owner(pid)] += 1
                    why.append(pid)
            if why:
                n_bad += 1
                if len(examples) < GATE_EXAMPLES:
                    examples.append(f"{coll}/{f.get('id')}: {', '.join(why)}")
    return {"n_checked": n_checked, "n_bad": n_bad, "missing": missing, "empty": empty,
            "foreign": dict(foreign), "examples": examples}


def _read_source_genomes(gdir, asv):
    """provenance/{asv}.json 'genomes' (the built source genomes of the ASV)."""
    with open(os.path.join(gdir, "provenance", f"{asv}.json"), "rb") as fh:
        prov = _json_bytes(fh.read())
    return prov.get("genomes") if isinstance(prov, dict) else None


def check_built_genome(job):
    """The feature gate for one built ASV: job = (gdir, asv, manifest entry 'genomes' or None).
    Returns (asv, result): result is patric_id_violations' dict plus 'sources', or
    {'error': ...} when the genome cannot be checked (no provenance source list, a list that
    differs from the manifest entry's, or an unreadable genome file); an error refuses the build
    too. Top-level, so a process pool can run it."""
    gdir, asv, man_genomes = job
    try:
        srcs = _read_source_genomes(gdir, asv)
    except (OSError, ValueError) as exc:
        return asv, {"error": f"provenance/{asv}.json unreadable ({type(exc).__name__}: {exc})"}
    if not (isinstance(srcs, list) and srcs and all(isinstance(g, str) and g for g in srcs)):
        return asv, {"error": f"provenance/{asv}.json has no source-genome list ('genomes')"}
    if man_genomes is not None and list(man_genomes) != srcs:
        return asv, {"error": f"provenance/{asv}.json 'genomes' differs from the manifest entry's"}
    try:
        with open(os.path.join(gdir, f"{asv}.json"), "rb") as fh:
            genome = _json_bytes(fh.read())
    except (OSError, ValueError) as exc:
        return asv, {"error": f"{asv}.json unreadable ({type(exc).__name__}: {exc})"}
    res = patric_id_violations(genome, srcs)
    res["sources"] = len(srcs)
    return asv, res


def require_own_patric_ids(gdir, manifest, workers=8, log=print):
    """The FEATURE GATE: refuse (SystemExit, exit status 1) a build in which any built genome
    holds a feature whose PATRIC_id is missing, empty, or not 'fig|<g>.' for one of that ASV's
    source genomes (check_built_genome; every 'built' manifest entry, whatever --limit or the
    workspace already holds). The genomes are streamed one file at a time, by `workers`
    processes (serially for small builds). main() calls it after require_bvbrc_only_build and
    before reading the token or contacting KBase, so a refused build uploads and logs nothing.
    Returns {n_genomes, n_checked} when every genome passes."""
    jobs = [(gdir, a, (e.get("genomes") if isinstance(e.get("genomes"), list) else None))
            for a, e in manifest.items()
            if not a.startswith("_") and isinstance(e, dict) and e.get("status") == "built"]
    t0 = time.perf_counter()
    if workers and workers > 1 and len(jobs) >= 4 * workers:
        with concurrent.futures.ProcessPoolExecutor(max_workers=workers) as ex:
            results = list(ex.map(check_built_genome, jobs, chunksize=8))
    else:
        results = [check_built_genome(j) for j in jobs]
    bad = [(a, r) for a, r in results if r.get("error") or r.get("n_bad")]
    n_checked = sum(r.get("n_checked", 0) for _, r in results)
    if bad:
        lines = []
        for a, r in bad[:20]:
            if r.get("error"):
                lines.append(f"  {a}: cannot be checked: {r['error']}")
                continue
            kinds = []
            if r["foreign"]:
                kinds.append("another genome's id " + ", ".join(
                    f"{k} x{v}" for k, v in sorted(r["foreign"].items(), key=lambda kv: -kv[1])[:5]))
            if r["empty"]:
                kinds.append(f"empty id x{r['empty']}")
            if r["missing"]:
                kinds.append(f"no PATRIC_id x{r['missing']}")
            eg = f"; e.g. {'; '.join(r['examples'])}" if r["examples"] else ""
            lines.append(f"  {a}: {r['n_bad']} of {r['n_checked']} features/CDSs "
                         f"({'; '.join(kinds)}){eg}")
        more = f"\n  ... and {len(bad) - 20} more" if len(bad) > 20 else ""
        raise SystemExit(
            f"[upload] REFUSED: {len(bad)} of {len(jobs)} built genomes in {gdir} hold features "
            f"whose PATRIC_id is not 'fig|<g>.' for one of the ASV's source genomes "
            f"(provenance/<asv>.json 'genomes'), i.e. not BV-BRC's own PATRIC annotation of a "
            f"source genome:\n" + "\n".join(lines) + more + "\nRebuild those ASVs with the "
            f"current builder (build_synthetic_genomes.bvbrc_only keeps a genome's own PATRIC ids "
            f"only, since 2026-11-13). Nothing was uploaded.")
    log(f"[upload] feature gate: {len(jobs)} built genomes, {n_checked} features/CDSs, every "
        f"PATRIC_id is a source genome's own ({time.perf_counter() - t0:.0f}s)")
    return {"n_genomes": len(jobs), "n_checked": n_checked}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome-dir", required=True)
    ap.add_argument("--ws-id", type=int, required=True)
    ap.add_argument("--selection", default=None, help="selection audit JSON (default: manifest _meta.selection)")
    ap.add_argument("--provenance", default=None, help="mapping_provenance.json (default: <genome-dir>/mapping_provenance.json)")
    ap.add_argument("--batch-mb", type=float, default=100.0)
    ap.add_argument("--batch-n", type=int, default=25)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--no-genomeset", action="store_true")
    ap.add_argument("--genomeset-name", default=None)
    ap.add_argument("--genomeset-description", default=None,
                    help="GenomeSet description (default: the codiffusion V4-V5 text)")
    ap.add_argument("--strip-extra-keys", action="store_true", help="drop the raw prFBA per-feature keys "
                    "(they are always duplicated into inference_data)")
    ap.add_argument("--gate-workers", type=int, default=8,
                    help="processes for the feature gate's pass over the build (default 8)")
    args = ap.parse_args()

    gdir = os.path.abspath(args.genome_dir)
    manifest = json.load(open(os.path.join(gdir, "synthetic_genome_manifest.json")))
    man_meta = manifest.get("_meta", {})
    pfunc_meta = man_meta.get("pfunc")       # None for builds made before E3 -> k/n text
    annotation = man_meta.get("annotation")  # None for builds made before the BV-BRC-only rule
    require_bvbrc_only_build(man_meta, gdir)  # exits non-zero, before any token/network use
    require_own_patric_ids(gdir, manifest, workers=args.gate_workers)   # likewise
    sel_path = os.path.abspath(args.selection or man_meta["selection"])
    sel = json.load(open(sel_path)); sel_meta = sel.pop("_meta", {})
    hits_path = sel_meta.get("generated_from")
    prov_path = args.provenance or os.path.join(gdir, "mapping_provenance.json")
    prov = json.load(open(prov_path)) if os.path.exists(prov_path) else {}
    commit = ((prov.get("pipeline") or {}).get("git") or {}).get("commit") or "unknown"

    token = open(TOKEN_FILE).read().strip()
    ws = WS(token)
    info = ws.call("Workspace.get_workspace_info", [{"id": args.ws_id}])
    print(f"[ws] {info[1]} (id {info[0]}, owner {info[2]}, {info[8].get('narrative_nice_name')}) | prFBA {commit[:12]}")
    existing = ws.existing(args.ws_id, GENOME_TYPE)
    print(f"[ws] {len(existing)} Genome objects already present")

    log_path = os.path.join(gdir, "kbase_upload_log.json")
    log = json.load(open(log_path)) if os.path.exists(log_path) else {}
    log["_meta"] = {"workspace_id": args.ws_id, "workspace_name": info[1], "prFBA_commit": commit,
                    "genome_dir": gdir, "selection_audit": sel_path, "started": log.get("_meta", {}).get("started") or now_kbase(),
                    "command": " ".join(sys.argv)}
    for a, upa in existing.items():           # adopt anything already there (resume)
        if a in manifest and a not in log:
            log[a] = {"name": a, "upa": upa, "adopted": True}

    todo = [a for a in manifest if not a.startswith("_") and manifest[a].get("status") == "built"
            and a not in existing and not (a in log and log[a].get("upa"))]
    if args.limit:
        todo = todo[:args.limit]
    print(f"[upload] {len(todo)} genomes to upload (batch <= {args.batch_mb} MB / {args.batch_n} objects)"
          + (" [DRY RUN]" if args.dry_run else ""), flush=True)

    stats = {"n": 0, "bytes": 0, "calls": 0, "failed": 0, "stripped": 0}
    batch, bsize, t_all = [], 0, time.perf_counter()
    keep = not args.strip_extra_keys
    for i, asv in enumerate(todo, 1):
        raw = json.load(open(os.path.join(gdir, f"{asv}.json")))
        s = sel.get(asv, {}); m = manifest[asv]
        rel_ab = s.get("rel_ab")
        data = to_kbase_genome(raw, asv=asv, sel=s, man=m, commit=commit, keep_extra_keys=keep,
                               pfunc_meta=pfunc_meta, annotation=annotation)
        name = asv if not asv.isdigit() else f"ASV_{asv}"
        obj = {"type": GENOME_TYPE, "name": name, "data": data,
               "meta": ws_meta(asv, s, m, commit, rel_ab, pfunc_meta=pfunc_meta,
                               annotation=annotation),
               "provenance": ws_provenance(asv, s, prov, commit, hits_path, sel_path, gdir,
                                           pfunc_meta=pfunc_meta, annotation=annotation)}
        nbytes = len(json.dumps(data))
        if args.dry_run:
            print(f"[dry] {name} {nbytes/1e6:.1f} MB meta={obj['meta']}")
            continue
        if batch and (bsize + nbytes > args.batch_mb * 1e6 or len(batch) >= args.batch_n):
            save_batch(ws, args.ws_id, batch, log, keep, stats)
            json.dump(log, open(log_path, "w"), indent=1)
            batch, bsize = [], 0
            el = time.perf_counter() - t_all
            print(f"[upload] {stats['n']}/{len(todo)} | {stats['bytes']/1e9:.2f} GB | {el/60:.1f} min | "
                  f"eta {el/max(stats['n'],1)*(len(todo)-stats['n'])/60:.0f} min | failed={stats['failed']}", flush=True)
        batch.append({"asv": asv, "obj": obj, "nbytes": nbytes, "raw": raw, "sel": s, "man": m, "commit": commit,
                      "pfunc_meta": pfunc_meta, "annotation": annotation})
        bsize += nbytes
    if batch and not args.dry_run:
        save_batch(ws, args.ws_id, batch, log, keep, stats)
        json.dump(log, open(log_path, "w"), indent=1)
    print(f"[upload] done: saved={stats['n']} failed={stats['failed']} stripped_extra_keys={stats['stripped']} "
          f"calls={stats['calls']} {stats['bytes']/1e9:.2f} GB in {(time.perf_counter()-t_all)/60:.1f} min", flush=True)

    # ---- GenomeSet grouping every uploaded genome ----
    if args.no_genomeset or args.dry_run:
        return
    elements = {}
    for a, e in log.items():
        if a.startswith("_") or not e.get("upa"):
            continue
        s = sel.get(a, {}); sel_ = s.get("selected") or []
        elements[a] = {"ref": e["upa"], "metadata": {
            "tier": str(s.get("tier")), "best_identity": str(s.get("best_identity")),
            "anchor_genome": str(sel_[0].get("genome_id")) if sel_ else "",
            "source_genomes": "|".join(str(x.get("genome_id")) for x in sel_)[:800],
            "rel_abundance": str(s.get("rel_ab"))}}
    if not elements:
        return
    gs_name = args.genomeset_name or f"codiffusion_synthetic_genomes_prFBA_{commit[:7]}"
    n_abst = man_meta.get("abstain")
    gs = {"description": args.genomeset_description or (f"prFBA synthetic genomes for the codiffusion bioreactor 16S V4-V5 ASVs: {len(elements)} of "
                          f"{man_meta.get('n_asvs')} ASVs ({n_abst} abstained: no BV-BRC reference passed the selection "
                          f"gates). prFBA commit {commit}; selection = exact-tie de-duplicated union + P3 gene tie-break; "
                          f"functions = {P.functions_description(man_meta)}; "
                          f"see mapping_provenance.json / abstained_asvs.md in {gdir}."),
          "elements": elements}
    gs_prov = [{"time": now_kbase(), "service": "prFBA", "service_ver": commit, "method": "upload_kbase_genomes.py",
                "script": "upload_kbase_genomes.py", "script_ver": commit, "script_command_line": " ".join(sys.argv),
                "description": gs["description"], "input_ws_objects": [e["ref"] for e in elements.values()],
                "custom": {"mapping_provenance": prov_path, "selection_audit": sel_path, "alignment_hits": str(hits_path)}}]
    gs_meta = {"n_genomes": str(len(elements)), "n_asvs": str(man_meta.get("n_asvs")), "n_abstained": str(n_abst),
               "prFBA_commit": commit[:12], "genome_dir": gdir[:800]}
    info = ws.call("Workspace.save_objects", [{"id": args.ws_id, "objects": [
        {"type": GENOMESET_TYPE, "name": gs_name, "data": gs, "meta": gs_meta, "provenance": gs_prov}]}])[0]
    log["_genomeset"] = {"name": info[1], "upa": f"{info[6]}/{info[0]}/{info[4]}", "n_elements": len(elements)}
    json.dump(log, open(log_path, "w"), indent=1)
    print(f"[upload] GenomeSet {info[1]} -> {log['_genomeset']['upa']} ({len(elements)} genomes)", flush=True)


if __name__ == "__main__":
    main()
