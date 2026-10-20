#!/usr/bin/env python
"""
build_gtdb_bvbrc_map.py — map every GTDB r226 representative accession (the
genomes in the skani GTDB database) to a BV-BRC genome_id.

Why
---
GTDB skani hits (`mag_skani_hits.py search --db-kind gtdb`) name each reference
by its NCBI assembly accession (GCA_/GCF_). The PGFam gene provider, the
prefetch and the synthetic-genome merge need numeric BV-BRC genome_ids. This
script builds that lookup once, over the whole representative universe, so a
rerun with a different top-k or new MAGs needs no new API calls. It does NOT
change any selection: nothing reads the map until the integration steps listed
in the stats sidecar are made (see `mag_skani_hits.load_acc_map`).

Construction (design: scout_accmap.md, 2026-10-20)
--------------------------------------------------
1. universe  every GC[AF]_#########.# accession in skani's index.db
2. metadata  GTDB r226 bac120 + ar53 metadata, downloaded and MD5-checked
             against the release MD5SUM.txt; if that fails, the local r232
             files (`--fallback-meta-dir`) -- the stats record which was used.
             Per accession: GenBank accession, ncbi_taxid, ncbi_genome_category,
             ncbi_biosample, genome_size.
3. pass 1    in(assembly_accession,(...)) over {the accession, its prefix-swapped
             form (same version), the metadata GenBank accession}
4. pass 2    only for accessions still without a usable record (usable = not
             Deprecated/Plasmid and passing the length gate below):
             (a) version-free wildcard on both prefixes,
                 or(eq(assembly_accession,GCA_<digits>.*),eq(...,GCF_<digits>.*))
             (b) then in(biosample_accession,(...)) on the metadata BioSample
   pass 3    (added after verification; --no-pass3 disables) only for accessions with
             no usable record after the guard and the gate below:
             or(and(eq(taxon_id,T),ge(genome_length,lo),le(genome_length,hi)),...) on the
             GTDB ncbi_taxid (and its merged.dmp target) with a +-1% window on the GTDB
             genome_size (match_type taxid_length). A record is accepted only if it is
             the ONLY live (non-Deprecated/Plasmid) record in that taxid+length window,
             is neither chosen nor accession-claimed by another representative, and
             carries no different assembly accession (a different assembly core is a
             different NCBI assembly); a record accepted for two representatives is
             rejected for both. "Only survivor after exclusions" is not enough: on a
             catch-all taxid (e.g. 1869227) the window holds dozens of other MAGs and
             the one without an assembly accession is just another MAG. Weakest
             evidence: no shared accession or BioSample.
5. choose    among a representative's candidates, in order:
             (1) drop genome_status Deprecated and Plasmid records;
             (2) match strength exact > prefix_swap > metadata_gca >
                 wildcard_core > biosample > taxid_length;
             (3) genome_id taxid prefix == GTDB ncbi_taxid (directly or through
                 NCBI merged.dmp);
             (4) genome_quality Good > missing > Poor;
             (5) non-empty entry in the PGFam cache;
             (6) genome_status Complete > WGS > other;
             (7) highest cds;  (8) lowest genome_id.
             Every candidate (Deprecated/Plasmid included) is kept in the row.
   guard     (added after verifying pass 2; --no-biosample-guard disables) a
             biosample-matched record is ineligible if it carries ANOTHER
             representative's assembly accession (one GTDB representative is one
             assembly; a metagenome BioSample can hold many MAGs), or if its taxid
             is not consistent with the GTDB ncbi_taxid -- unless it has the same
             BioSample and a genome_length within 1% of GTDB's genome_size (an NCBI
             taxon rename, e.g. GCF_037478005.1 -> 233316.43, identical 5,012,782 bp).
   length    (added after verification; --no-length-gate disables) every
   gate      non-Deprecated/Plasmid candidate is compared with the GTDB genome_size,
             d = |genome_length - genome_size| / genome_size, and is ineligible when
             d > 0.05 for wildcard_core, metadata_gca and biosample matches (version
             moves and shared BioSamples can point at another sequence), or d > 0.10
             for exact/prefix_swap matches whose genome_id taxid prefix is a
             'mismatch' (crossed BV-BRC records, e.g. GCF_001182745.1 Nocardia ->
             562.7382, a 1.8 Mb E. coli record). A record without genome_length (or a
             representative without genome_size) cannot be tested and is kept.

API etiquette: bulk POSTs with a descriptive User-Agent (BV-BRC answers 403 to
the default Python one), at most 4 concurrent requests, retries with
exponential backoff, and a batch that keeps failing aborts the run. Only a
successful answer is written, so a failed batch is never cached as empty. A
batch whose answer is truncated (Content-Range total > rows, or rows == limit)
is split in halves and re-queried. Rerunning resumes from the raw batch files.

Outputs (--outdir, default reference_maps/)
-------------------------------------------
gtdb_r226_acc_to_bvbrc.tsv.gz      one row per universe accession, columns:
    gtdb_accession            key (as emitted by mag_skani_hits.ref_id)
    gtdb_genbank_accession    GTDB metadata ncbi_genbank_assembly_accession
    acc_core                  the 9 accession digits (prefix and version dropped)
    gtdb_ncbi_taxid           GTDB metadata ncbi_taxid
    ncbi_genome_category      none (isolate) | derived from metagenome |
                              derived from single cell | derived from environmental sample
    bvbrc_genome_id           chosen BV-BRC genome_id; '' when none is usable
    bvbrc_assembly_accession  that record's assembly_accession
    match_type                exact | prefix_swap | metadata_gca | wildcard_core |
                              biosample | taxid_length (pass 3) for a chosen record;
                              when nothing is usable, in this precedence:
                              deprecated_only (a Deprecated/Plasmid record is an
                              accession match within 5% of the GTDB genome_size, or
                              only Deprecated/Plasmid records were found) |
                              length_rejected (the length gate rejected the rest) |
                              biosample_rejected (the biosample guard rejected the
                              rest) | unmapped (no record by accession or BioSample,
                              and no unique taxid+length record).
                              bvbrc_genome_id is '' for all four.
    n_candidates              number of BV-BRC records found (all statuses)
    candidates                every record, 'genome_id:genome_status;...', chosen first,
                              then the other usable ones in tie-break order, then those
                              the guard / length gate rejected (see rejected_candidates),
                              then the Deprecated/Plasmid ones.
                              Membership checks must use bvbrc_genome_id, not this.
    genome_status, genome_quality, bvbrc_taxon_id, cds   of the chosen record
    taxid_consistent          exact | merged (equal through merged.dmp) | mismatch |
                              missing (no GTDB taxid), comparing the chosen genome_id's taxid prefix with
                              gtdb_ncbi_taxid ('' when unmapped)
    in_pgfam_cache            nonempty | empty | absent (chosen genome_id in the
                              PGFam cache at build time; '' when unmapped)
    queried_at                UTC time of the API answer that returned the chosen record
    gtdb_genome_size          GTDB metadata genome_size (bp)
    bvbrc_genome_length       the chosen record's genome_length ('' when unmapped/absent)
    length_ratio              bvbrc_genome_length / gtdb_genome_size, 4 decimals
                              ('' when either is missing)
    bvbrc_label_taxid_consistent  exact | merged | mismatch | missing, comparing the
                              chosen record's BV-BRC taxon_id (its label) with
                              gtdb_ncbi_taxid ('' when unmapped). Some sequence-correct
                              records carry another strain's label, so integration must
                              take taxon_id from gtdb_ncbi_taxid, not from BV-BRC.
    rejected_candidates       'genome_id:reason;...' for every candidate made ineligible:
                              claimed_by_another_rep | taxid_mismatch (biosample guard),
                              length (length gate). Pass-3 records that were not accepted
                              are records of other genomes and are in neither column (see
                              stats.pass3_taxid_length: taxid_length_taken |
                              taxid_length_other_assembly | taxid_length_shared |
                              taxid_length_ambiguous)
    (the last five columns were appended in v2; the first 17 keep their order)
gtdb_r226_acc_to_bvbrc.stats.json  provenance, coverage, match types, one-to-many,
                                   taxid consistency, runtime, integration to-do
raw/gtdb_<release>_reps_meta.tsv.gz  the metadata fields used, per accession
raw/{pass1_exact,pass2_wildcard,pass2_biosample,pass3_taxid_length}/bNNNNN.json.gz
                                   every API answer, with its select field list; a batch
                                   is reused only if both its values and its fields match

Usage
-----
  python build_gtdb_bvbrc_map.py --meta-dir <dir for the r226 metadata> \\
      --report-set real_selected=<accs.txt> --report-mag-selection real=<selection.json>
"""

import argparse
import collections
import csv
import datetime
import gzip
import hashlib
import http.client
import json
import math
import mmap
import os
import random
import re
import sys
import threading
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
API = "https://www.bv-brc.org/api"
UA = ("prFBA-gtdb-bvbrc-map/1.0 (academic research; GTDB-to-BV-BRC accession map; "
      "+https://github.com/freiburgermsu/prfba)")
GTDB_URL = "https://data.gtdb.ecogenomic.org/releases/release226/226.0/"
META_FILES = ("bac120_metadata_r226.tsv.gz", "ar53_metadata_r226.tsv.gz")
FALLBACK_FILES = ("bac120_metadata_r232.tsv.gz", "ar53_metadata_r232.tsv.gz")
META_COLS = ["accession", "gtdb_representative", "ncbi_genbank_assembly_accession",
             "ncbi_taxid", "ncbi_species_taxid", "ncbi_genome_category", "ncbi_biosample",
             "genome_size"]
FIELDS = ["genome_id", "assembly_accession", "genome_status", "genome_quality", "taxon_id",
          "cds", "genome_name", "checkm_completeness", "checkm_contamination",
          "biosample_accession", "date_inserted", "genome_length"]
ACC_FIND = re.compile(rb"GC[AF]_[0-9]{9}\.[0-9]+")
GID_RE = re.compile(r"^\d+\.\d+$")

MATCH_RANK = {"exact": 0, "prefix_swap": 1, "metadata_gca": 2, "wildcard_core": 3,
              "biosample": 4, "taxid_length": 5}
ACCESSION_MATCHES = ("exact", "prefix_swap", "metadata_gca", "wildcard_core")
LEN_WEAK_TYPES = ("wildcard_core", "metadata_gca", "biosample")  # gated at --len-tol-weak
GUARD_REASONS = ("claimed_by_another_rep", "taxid_mismatch")
UNRESOLVED_TYPES = ("deprecated_only", "length_rejected", "biosample_rejected", "unmapped")
DROP_STATUS = {"Deprecated", "Plasmid"}
QUALITY_RANK = {"Good": 0, None: 1, "": 1, "Poor": 3}  # any other value ranks 2
STATUS_RANK = {"Complete": 0, "WGS": 1}                # any other value ranks 2
TIE_STEPS = ["match_strength", "taxid_prefix", "genome_quality", "pgfam_cache",
             "complete_before_wgs", "highest_cds", "lowest_genome_id"]
COLUMNS = ["gtdb_accession", "gtdb_genbank_accession", "acc_core", "gtdb_ncbi_taxid",
           "ncbi_genome_category", "bvbrc_genome_id", "bvbrc_assembly_accession", "match_type",
           "n_candidates", "candidates", "genome_status", "genome_quality", "bvbrc_taxon_id",
           "taxid_consistent", "cds", "in_pgfam_cache", "queried_at",
           # appended in v2 (2026-10-20, after verification); the 17 above keep their order
           "gtdb_genome_size", "bvbrc_genome_length", "length_ratio",
           "bvbrc_label_taxid_consistent", "rejected_candidates"]

INTEGRATION_TODO = [
    "mag_skani_hits.py: add --acc-map to the search subparser and load it after the taxonomy; "
    "in each top-k entry set genome_id to the mapped BV-BRC id (fall back to the accession), "
    "taxon_id to int(gtdb_ncbi_taxid) -- NOT the genome_id prefix and NOT bvbrc_taxon_id, "
    "because crossed-label BV-BRC records (bvbrc_label_taxid_consistent == 'mismatch' on a "
    "length-matched accession match, e.g. GCF_002020805.1 -> 1918949.3) carry another "
    "strain's taxid -- and add gtdb_accession, bvbrc_mapped, match_type and "
    "is_mag = (ncbi_genome_category != 'none'). Keep ref_id/best_ref as the accession "
    "(GTDB-Tk concordance and same_ref_as_gtdbtk depend on it); add best_ref_bvbrc and "
    "map-coverage counts to run_stats; leave the lineage as GTDB.",
    "select_references.py L570: is_low = r['uncultured'] or bool(r['_h'].get('is_mag')) or "
    "_is_mag(...). Numeric mapped ids escape the _MAG regex, so without this the real-MAG "
    "handfuls grow past max_low_conf=1 and the reported GTDB selections change.",
    "Guard against non-numeric ids: skip any id not matching ^\\d+\\.\\d+$ in "
    "bvbrc_gene_provider.provider (select_references.py) and in "
    "prefetch_gene_families.candidate_genome_ids, so fetch_pgfams('GCA_...') can never "
    "cache [] in the shared PGFam cache.",
    "Unmapped-hit policy as an explicit flag: take tier and taxonomy from the full GTDB hit "
    "list, build the gene union from mapped hits only, and record unmapped_skipped (today "
    "an unmapped anchor is excluded at stage 7 as no_gene_set and the next hit becomes the "
    "anchor, which changes the taxonomic call). Every row with bvbrc_genome_id == '' is "
    "unmapped: match_type unmapped, deprecated_only, biosample_rejected and length_rejected "
    "alike. Decide whether taxid_length rows (taxid + length only, no shared accession or "
    "BioSample) are used or treated as unmapped.",
    "build_synthetic_genomes.py: add --knobs-json and thread it into exact_selection and "
    "PF.populate_cache (prefetch_gene_families.py L92 calls candidate_genome_ids with the 16S "
    "default knobs); add the User-Agent to _api_json (bare urllib is refused with 403).",
    "PGFam prefetch for the mapped best-hit union (and the top-k union if needed) before any "
    "rerun that uses the gene provider.",
    "GTDB arm scoring: translate each hit's ref to the mapped id and call "
    "score_mag.score_query (NCBI-lineage rank cells, PGFam capture, same-taxon exclusion); "
    "needs an int-encoded PGFam index of the mapped ids and truth files for the k-means and "
    "random arms. Use bvbrc_genome_id for any `src_gid ==` check; the candidates column "
    "also lists rejected (rejected_candidates) and Deprecated/Plasmid records, so a bare "
    "`src_gid in candidates` test is wrong.",
    "Manuscript: update prFBA.tex L368-369 and prFBA_andrew_voice.tex L358-359 once the "
    "integrated results exist (user sign-off needed: the integration changes reported "
    "results).",
]

INTEGRATION_NOTES = [
    "taxon_id must come from gtdb_ncbi_taxid. The chosen record's own label (bvbrc_taxon_id; "
    "compared in bvbrc_label_taxid_consistent) and its genome_id prefix (taxid_consistent) can "
    "belong to another strain: crossed-label BV-BRC records keep the right sequence (length "
    "equals GTDB's genome_size) but another strain's taxon_id and BioSample.",
    "Unmapped = bvbrc_genome_id == '' (match_type unmapped, deprecated_only, biosample_rejected, "
    "length_rejected).",
    "candidates lists every record found, including those in rejected_candidates and "
    "Deprecated/Plasmid ones; membership checks must use bvbrc_genome_id.",
    "taxid_length (pass 3) rows rest on taxid + genome length within 1% only; flag or "
    "exclude them where a same-assembly guarantee matters.",
    "length_ratio is bvbrc_genome_length / gtdb_genome_size of the chosen record; the gate "
    "thresholds are in stats.length_gate.",
]


# ----------------------------------------------------------------------- utils
def log(msg):
    print(f"[gtdb_bvbrc_map {time.strftime('%H:%M:%S')}] {msg}", file=sys.stderr, flush=True)


def utcnow():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def acc_core(acc):
    """GCA_003722315.1 / GCF_003722315.1 -> 003722315 (same assembly, either prefix)."""
    if not acc or "_" not in acc:
        return None
    return acc.split("_", 1)[1].split(".")[0]


def swap_prefix(acc):
    return ("GCF_" if acc.startswith("GCA_") else "GCA_") + acc[4:]


def gid_key(gid):
    a, b = gid.split(".", 1)
    return int(a), int(b)


def wilson(k, n, z=1.959964):
    if n == 0:
        return None
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(100 * (c - h), 4), round(100 * (c + h), 4)]


def pct(k, n):
    return round(100 * k / n, 4) if n else None


def write_json_gz(path, obj):
    tmp = path + ".tmp"
    with open(tmp, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(json.dumps(obj, separators=(",", ":")).encode())
    os.replace(tmp, path)


def read_json_gz(path):
    with gzip.open(path, "rb") as fh:
        return json.loads(fh.read())


def md5sum(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


# ----------------------------------------------------------------- BV-BRC API
class Transient(Exception):
    """A response that should be retried (server error, suspicious empty answer...)."""


def rql_post(collection, query, retries=6, timeout=300, empty_is_error=False):
    """POST an RQL query; return (rows, content_range_total or None). Retries transient
    failures with exponential backoff; 4xx other than 408/429 raise at once."""
    req = urllib.request.Request(
        f"{API}/{collection}/", data=query.encode(),
        headers={"User-Agent": UA, "Accept": "application/json",
                 "Content-Type": "application/rqlquery+x-www-form-urlencoded"})
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                total = None
                m = re.search(r"/(\d+)\s*$", resp.headers.get("Content-Range") or "")
                if m:
                    total = int(m.group(1))
                rows = json.loads(resp.read())
            if not isinstance(rows, list):
                raise Transient(f"non-list answer: {str(rows)[:200]}")
            if empty_is_error and not rows:
                raise Transient("empty answer for a batch that cannot plausibly be empty")
            return rows, total
        except urllib.error.HTTPError as exc:
            if exc.code not in (408, 429) and exc.code < 500:
                raise
            err = exc
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError,
                http.client.HTTPException, ValueError, Transient) as exc:
            err = exc
        if attempt == retries - 1:
            raise RuntimeError(f"{collection} query failed after {retries} attempts: {err}")
        wait = 2 ** (attempt + 1) + random.uniform(0, 1)
        log(f"{collection} batch failed ({err}); retry in {wait:.0f}s")
        time.sleep(wait)


def fetch_values(build_query, values, limit, empty_is_error, depth=0,
                 single_truncation_ok=False):
    """Run build_query(values, limit); if the answer is truncated, split the values in
    halves and recurse. Returns (rows, n_splits, truncated_values). A single value that
    still returns >= limit rows aborts the run, unless single_truncation_ok (pass 3, where
    such a value is recorded as truncated and its representative treated as ambiguous)."""
    q = build_query(values, limit)
    rows, total = rql_post("genome", q, empty_is_error=empty_is_error)
    truncated = len(rows) >= limit or (total is not None and total > len(rows))
    if not truncated:
        return rows, 0, []
    if len(values) == 1:
        if single_truncation_ok:
            log(f"single value returns >= {limit} rows (total={total}); kept as truncated: "
                f"{values[0]}")
            return rows, 0, [values[0]]
        raise RuntimeError(f"single value returns >= {limit} rows: {values[0]}")
    log(f"truncated answer ({len(rows)} rows, total={total}) for {len(values)} values; "
        f"splitting")
    mid = len(values) // 2
    r1, s1, t1 = fetch_values(build_query, values[:mid], limit, False, depth + 1,
                              single_truncation_ok)
    r2, s2, t2 = fetch_values(build_query, values[mid:], limit, False, depth + 1,
                              single_truncation_ok)
    seen, out = set(), []
    for r in r1 + r2:
        if r.get("genome_id") not in seen:
            seen.add(r.get("genome_id"))
            out.append(r)
    return out, 1 + s1 + s2, t1 + t2


def run_pass(name, values, batch, build_query, raw_dir, workers, limit, pause,
             empty_is_error, fields, single_truncation_ok=False):
    """Query `values` in batches; each successful answer is written to
    raw_dir/name/bNNNNN.json.gz together with the select field list, and reused on a rerun
    only if both its values and its fields are identical (a batch saved before a field was
    added would otherwise be reused with that field missing from every row)."""
    pdir = os.path.join(raw_dir, name)
    os.makedirs(pdir, exist_ok=True)
    chunks = [values[i:i + batch] for i in range(0, len(values), batch)]
    t0 = time.time()
    stats = dict(n_values=len(values), n_batches=len(chunks), reused=0, queried=0,
                 n_rows=0, truncation_splits=0, stale_fields_requeried=0)
    results = [None] * len(chunks)
    todo = []
    for i, chunk in enumerate(chunks):
        path = os.path.join(pdir, f"b{i:05d}.json.gz")
        if os.path.exists(path):
            rec = read_json_gz(path)
            if rec.get("values") == chunk and rec.get("fields") == list(fields):
                results[i] = rec
                stats["reused"] += 1
                continue
            if rec.get("values") == chunk:
                stats["stale_fields_requeried"] += 1
        todo.append(i)
    log(f"{name}: {len(values)} values in {len(chunks)} batches; reused {stats['reused']}, "
        f"to query {len(todo)}")
    lock = threading.Lock()
    done = [0]

    def work(i):
        chunk = chunks[i]
        t = time.time()
        rows, splits, trunc = fetch_values(build_query, chunk, limit,
                                           empty_is_error and len(chunk) >= 20,
                                           single_truncation_ok=single_truncation_ok)
        rec = dict(pass_name=name, batch=i, values=chunk, fields=list(fields),
                   queried_at=utcnow(), secs=round(time.time() - t, 3), n_rows=len(rows),
                   splits=splits, truncated_values=trunc, rows=rows)
        write_json_gz(os.path.join(pdir, f"b{i:05d}.json.gz"), rec)
        if pause:
            time.sleep(pause)
        with lock:
            results[i] = rec
            done[0] += 1
            if done[0] % 25 == 0 or done[0] == len(todo):
                log(f"{name}: {done[0]}/{len(todo)} batches ({time.time() - t0:.0f}s)")
        return i

    if todo:
        ex = ThreadPoolExecutor(max_workers=workers)
        futs = [ex.submit(work, i) for i in todo]
        try:
            for fut in futs:
                fut.result()  # re-raises: a failed batch aborts the run, nothing cached
        finally:
            ex.shutdown(wait=True, cancel_futures=True)
    stats["queried"] = len(todo)
    # a batch file beyond the current batch count is left over from an earlier, larger
    # value list; remove it so raw/ holds exactly the answers this build used
    stale = [f for f in os.listdir(pdir)
             if re.fullmatch(r"b\d{5}\.json\.gz", f) and int(f[1:6]) >= len(chunks)]
    for f in stale:
        os.remove(os.path.join(pdir, f))
    stats["stale_batch_files_removed"] = len(stale)
    stats["n_rows"] = sum(r["n_rows"] for r in results)
    stats["truncation_splits"] = sum(r.get("splits", 0) for r in results)
    stats["truncated_single_values"] = sum(len(r.get("truncated_values") or [])
                                           for r in results)
    stats["secs"] = round(time.time() - t0, 1)
    stats["api_secs_sum"] = round(sum(r["secs"] for r in results), 1)
    stamps = sorted(r["queried_at"] for r in results)
    stats["first_queried_at"] = stamps[0] if stamps else None
    stats["last_queried_at"] = stamps[-1] if stamps else None
    return results, stats


# ------------------------------------------------------------------ inputs
def load_universe(index_db):
    data = open(index_db, "rb").read()
    found = ACC_FIND.findall(data)
    accs = sorted({a.decode() for a in found})
    return accs, len(found)


def fetch_url(url, dest):
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    tmp = dest + ".part"
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=600) as resp, open(tmp, "wb") as out:
                for block in iter(lambda: resp.read(1 << 20), b""):
                    out.write(block)
            os.replace(tmp, dest)
            return
        except (urllib.error.URLError, TimeoutError, ConnectionError, OSError) as exc:
            if attempt == 3:
                raise
            log(f"download {url} failed ({exc}); retry")
            time.sleep(2 ** (attempt + 2))


def get_metadata_files(meta_dir, fallback_dir, offline):
    """Return (release, {file: info}, error). r226 is downloaded into meta_dir (or reused
    if present) and MD5-checked; on any failure fall back to the local r232 files."""
    info, err = {}, None
    try:
        if offline:
            raise RuntimeError("--offline given")
        os.makedirs(meta_dir, exist_ok=True)
        md5_path = os.path.join(meta_dir, "MD5SUM.txt")
        if not os.path.exists(md5_path):
            fetch_url(GTDB_URL + "MD5SUM.txt", md5_path)
        expected = {}
        for line in open(md5_path):
            parts = line.split()
            if len(parts) == 2:
                expected[os.path.basename(parts[1])] = parts[0]
        for fname in META_FILES:
            path = os.path.join(meta_dir, fname)
            downloaded = False
            if not (os.path.exists(path) and md5sum(path) == expected.get(fname)):
                log(f"downloading {GTDB_URL + fname}")
                fetch_url(GTDB_URL + fname, path)
                downloaded = True
            got = md5sum(path)
            if got != expected.get(fname):
                raise RuntimeError(f"MD5 mismatch for {fname}: {got} != {expected.get(fname)}")
            info[fname] = dict(path=os.path.abspath(path), url=GTDB_URL + fname, md5=got,
                               md5_verified=True, bytes=os.path.getsize(path),
                               downloaded_this_run=downloaded)
        return "r226", info, None
    except Exception as exc:  # noqa: BLE001 -- recorded, then the r232 fallback is used
        err = f"{type(exc).__name__}: {exc}"
        log(f"r226 metadata unavailable ({err}); falling back to {fallback_dir}")
    info = {}
    for fname in FALLBACK_FILES:
        path = os.path.join(fallback_dir, fname)
        info[fname] = dict(path=os.path.abspath(path), md5=md5sum(path), md5_verified=False,
                           bytes=os.path.getsize(path))
    return "r232", info, err


def load_metadata(release, files, universe):
    """{accession: meta dict} for the universe; join exact, else (r232 fallback) by digits."""
    import pandas as pd
    frames = [pd.read_csv(f["path"], sep="\t", usecols=META_COLS, dtype=str,
                          compression="gzip", keep_default_na=False)
              for f in files.values()]
    m = pd.concat(frames, ignore_index=True)
    n_rows = len(m)
    m["bare"] = m["accession"].str.replace(r"^(GB|RS)_", "", regex=True)
    recs = {}
    by_core = {}
    for row in m.itertuples(index=False):
        d = row._asdict()
        recs[d["bare"]] = d
        by_core.setdefault(acc_core(d["bare"]), []).append(d)
    meta, join = {}, collections.Counter()
    for acc in universe:
        d = recs.get(acc)
        how = "exact"
        if d is None:
            cands = by_core.get(acc_core(acc)) or []
            same = [c for c in cands if c["bare"][:3] == acc[:3]]
            d = (same or cands or [None])[0]
            how = "core" if d is not None else "none"
        join[how] += 1
        if d is None:
            meta[acc] = dict(meta_accession="", meta_join="none",
                             ncbi_genbank_assembly_accession="", ncbi_taxid="",
                             ncbi_species_taxid="", ncbi_genome_category="", ncbi_biosample="",
                             genome_size="", gtdb_representative="")
            continue
        meta[acc] = dict(meta_accession=d["accession"], meta_join=how,
                         **{k: (d[k] or "").strip() for k in META_COLS if k != "accession"})
    n_rep_flag = sum(1 for a in universe if meta[a].get("gtdb_representative") == "t")
    return meta, dict(release=release, n_metadata_rows=n_rows, join=dict(join),
                      universe_flagged_representative=n_rep_flag)


def scan_pgfam_cache(path):
    """Genome keys of the PGFam cache JSON, split into non-empty and empty, by scanning
    the file with mmap (it is several GB; json.load would need tens of GB)."""
    nonempty, empty, odd = set(), set(), 0
    with open(path, "rb") as fh, mmap.mmap(fh.fileno(), 0, access=mmap.ACCESS_READ) as mm:
        for pat in (b'": [', b'":['):
            pos = 0
            while True:
                j = mm.find(pat, pos)
                if j < 0:
                    break
                k = mm.rfind(b'"', max(0, j - 64), j)
                key = mm[k + 1:j].decode("ascii", "replace")
                nxt = mm[j + len(pat):j + len(pat) + 1]
                if GID_RE.match(key):
                    (empty if nxt == b"]" else nonempty).add(key)
                else:
                    odd += 1
                pos = j + len(pat)
            if nonempty or empty:
                break
    return nonempty, empty, odd


def load_merged(path):
    merged = {}
    with open(path) as fh:
        for line in fh:
            parts = [p.strip() for p in line.split("|")]
            if len(parts) >= 2 and parts[0].isdigit() and parts[1].isdigit():
                merged[int(parts[0])] = int(parts[1])
    return merged


def taxid_relation(a, b, merged):
    """exact | merged | mismatch | missing for two NCBI taxids (ints or None)."""
    if a is None or b is None:
        return "missing"
    if a == b:
        return "exact"
    ma, mb = merged.get(a, a), merged.get(b, b)
    if ma == b or mb == a or ma == mb:
        return "merged"
    return "mismatch"


def to_int(x):
    try:
        return int(str(x).strip())
    except (TypeError, ValueError):
        return None


def rel_len_dev(row, gsize):
    """Signed (genome_length - genome_size) / genome_size, or None if either is missing."""
    gl = to_int(row.get("genome_length"))
    if not gl or not gsize:
        return None
    return (gl - gsize) / gsize


def len_bin(d):
    if d is None:
        return "missing"
    d = abs(d)
    return ("0 identical" if d == 0 else "1 <=1%" if d <= 0.01 else "2 1-5%" if d <= 0.05
            else "3 5-10%" if d <= 0.10 else "4 >10%")


# ------------------------------------------------------------------- choose
def candidate_key(c, gtdb_taxid, merged, pg_nonempty):
    row = c["row"]
    gid = row["genome_id"]
    rel = taxid_relation(gid_key(gid)[0], gtdb_taxid, merged)
    cds = row.get("cds")
    cds = cds if isinstance(cds, int) else (to_int(cds) or -1)
    return (MATCH_RANK[c["match_type"]],
            0 if rel in ("exact", "merged") else 1,
            QUALITY_RANK.get(row.get("genome_quality"), 2),
            0 if gid in pg_nonempty else 1,
            STATUS_RANK.get(row.get("genome_status"), 2),
            -cds,
            gid_key(gid))


def choose(cands, gtdb_taxid, merged, pg_nonempty, rejected=None):
    """(chosen candidate or None, ordered candidate list, decisive step or None).
    `rejected` maps genome_id -> reason for the candidates the biosample guard, the length
    gate or the pass-3 filter made ineligible for this accession."""
    rejected = rejected or {}
    def dead(c):
        return c["row"].get("genome_status") in DROP_STATUS
    usable = [c for c in cands.values() if not dead(c) and c["row"]["genome_id"] not in rejected]
    rej = sorted((c for c in cands.values()
                  if not dead(c) and c["row"]["genome_id"] in rejected),
                 key=lambda c: gid_key(c["row"]["genome_id"]))
    dropped = sorted((c for c in cands.values() if dead(c)),
                     key=lambda c: gid_key(c["row"]["genome_id"]))
    keyed = sorted(((candidate_key(c, gtdb_taxid, merged, pg_nonempty), c) for c in usable),
                   key=lambda kc: kc[0])
    ordered = [c for _, c in keyed] + rej + dropped
    if not keyed:
        return None, ordered, None
    if len(keyed) == 1:
        reasons = {rejected[c["row"]["genome_id"]] for c in rej}
        step = ("length_gate" if "length" in reasons else
                "biosample_guard" if reasons else
                "drop_deprecated_plasmid" if dropped else None)
        return keyed[0][1], ordered, step
    k0, k1 = keyed[0][0], keyed[1][0]
    step = next(TIE_STEPS[i] for i in range(len(k0)) if k0[i] != k1[i])
    return keyed[0][1], ordered, step


def add_candidate(cands, row, match_type, queried_at, pass_name):
    gid = row.get("genome_id")
    if not gid or not GID_RE.match(str(gid)):
        return False
    old = cands.get(gid)
    if old is None or MATCH_RANK[match_type] < MATCH_RANK[old["match_type"]]:
        cands[gid] = dict(row=row, match_type=match_type, queried_at=queried_at,
                          pass_name=pass_name)
    return True


# ------------------------------------------------------------------- report
def group_coverage(accs, out_rows, universe_set):
    accs = list(dict.fromkeys(a.strip() for a in accs if a.strip()))
    inu = [a for a in accs if a in universe_set]
    rows = [out_rows[a] for a in inu]
    res = [r for r in rows if r["bvbrc_genome_id"]]
    gids = {r["bvbrc_genome_id"] for r in res}
    pg = {g: r["in_pgfam_cache"] for r in res for g in [r["bvbrc_genome_id"]]}
    return dict(
        n=len(accs), n_in_universe=len(inu),
        distinct_chosen_genome_ids=len(gids),
        distinct_chosen_pgfam=dict(collections.Counter(pg.values()).most_common()),
        not_in_universe=[a for a in accs if a not in universe_set][:50],
        resolved=len(res), pct_resolved=pct(len(res), len(inu)),
        wilson95=wilson(len(res), len(inu)),
        match_type=dict(collections.Counter(r["match_type"] for r in rows).most_common()),
        chosen_in_pgfam_cache=dict(collections.Counter(r["in_pgfam_cache"]
                                                       for r in res).most_common()),
        taxid_consistent=dict(collections.Counter(r["taxid_consistent"]
                                                  for r in res).most_common()))


def mag_selection_coverage(path, out_rows):
    sel = json.load(open(path))
    n_mags = n_sel = n_all = n_any = n_anchor = 0
    accs = set()
    for key, rec in sel.items():
        if key == "_meta" or not isinstance(rec, dict):
            continue
        n_mags += 1
        gids = [s.get("genome_id") for s in rec.get("selected") or []]
        if not gids:
            continue
        n_sel += 1
        accs.update(gids)
        mapped = [bool(out_rows.get(g, {}).get("bvbrc_genome_id")) for g in gids]
        n_all += all(mapped)
        n_any += any(mapped)
        n_anchor += mapped[0]
    return dict(path=os.path.abspath(path), n_mags=n_mags, n_with_selection=n_sel,
                n_distinct_selected_accessions=len(accs),
                mags_anchor_mapped=n_anchor, mags_all_selected_mapped=n_all,
                mags_any_selected_mapped=n_any,
                pct_of_all_mags_anchor_mapped=pct(n_anchor, n_mags),
                pct_of_selected_mags_anchor_mapped=pct(n_anchor, n_sel))


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--index-db", default=os.path.join(HERE, "skani_db",
                                                       "skani_gtdb_r226-v0.3", "index.db"))
    ap.add_argument("--meta-dir", required=True,
                    help="directory holding (or to receive) the GTDB r226 metadata files")
    ap.add_argument("--fallback-meta-dir", default=os.path.join(
        HERE, "region_validation", "picrust2_comparison", "data", "mag_arm", "gtdb"),
        help="local r232 metadata, used only if r226 cannot be downloaded/verified")
    ap.add_argument("--offline", action="store_true",
                    help="skip the r226 download and use the fallback metadata")
    ap.add_argument("--pgfam-cache", default=os.path.join(HERE, "bvbrc_cache",
                                                          "genome_gene_families.json"))
    ap.add_argument("--merged-dmp", default=os.path.join(HERE, "merged.dmp"))
    ap.add_argument("--outdir", default=os.path.join(HERE, "reference_maps"))
    ap.add_argument("--expect-n", type=int, default=143614,
                    help="expected universe size (abort if different; 0 disables)")
    ap.add_argument("--check-universe", help="optional accession list that must equal "
                                             "the index.db universe")
    ap.add_argument("--subset", type=int, default=0,
                    help="testing only: a seeded random subset of N accessions "
                         "(use a scratch --outdir and --expect-n 0)")
    ap.add_argument("--workers", type=int, default=4, help="concurrent requests (max 4)")
    ap.add_argument("--batch", type=int, default=1000, help="pass-1 strings per POST")
    ap.add_argument("--wild-batch", type=int, default=100, help="pass-2 wildcard terms/POST")
    ap.add_argument("--bios-batch", type=int, default=200, help="pass-2 biosamples per POST")
    ap.add_argument("--limit", type=int, default=3000, help="RQL limit per POST")
    ap.add_argument("--pause", type=float, default=0.2, help="seconds each worker waits "
                                                             "after a request")
    ap.add_argument("--report-set", action="append", default=[], metavar="NAME=PATH",
                    help="accession list (one per line) to report coverage for")
    ap.add_argument("--report-union", action="append", default=[],
                    metavar="NAME=SET1+SET2", help="coverage/demand of a union of --report-set "
                                                   "names")
    ap.add_argument("--no-biosample-guard", action="store_true",
                    help="keep biosample matches that carry another representative's "
                         "assembly or a mismatching taxid (off by default)")
    ap.add_argument("--no-length-gate", action="store_true",
                    help="keep candidates whose genome_length disagrees with the GTDB "
                         "genome_size (off by default)")
    ap.add_argument("--len-tol-weak", type=float, default=0.05,
                    help="length gate for wildcard_core/metadata_gca/biosample matches")
    ap.add_argument("--len-tol-accession", type=float, default=0.10,
                    help="length gate for exact/prefix_swap matches whose genome_id taxid "
                         "prefix is a mismatch")
    ap.add_argument("--len-tol-rename", type=float, default=0.01,
                    help="a same-BioSample biosample match within this length is kept "
                         "despite a taxid mismatch (NCBI taxon rename)")
    ap.add_argument("--len-tol-taxid", type=float, default=0.01,
                    help="pass-3 length window around the GTDB genome_size")
    ap.add_argument("--no-pass3", action="store_true",
                    help="skip pass 3 (GTDB ncbi_taxid + genome length)")
    ap.add_argument("--p3-batch", type=int, default=50, help="pass-3 terms per POST")
    ap.add_argument("--report-mag-selection", action="append", default=[],
                    metavar="NAME=PATH", help="mag_reference_selection.*.json to report "
                                              "MAG-weighted coverage for")
    a = ap.parse_args()
    workers = max(1, min(a.workers, 4))
    t_start = time.time()
    timings = {}
    os.makedirs(a.outdir, exist_ok=True)
    raw_dir = os.path.join(a.outdir, "raw")
    os.makedirs(raw_dir, exist_ok=True)

    # 1. universe
    t = time.time()
    universe, n_occ = load_universe(a.index_db)
    if a.subset:
        universe = sorted(random.Random(226).sample(universe, a.subset))
    universe_set = set(universe)
    log(f"universe: {len(universe)} accessions ({n_occ} occurrences) in {a.index_db}")
    if a.expect_n and len(universe) != a.expect_n:
        sys.exit(f"universe has {len(universe)} accessions, expected {a.expect_n}")
    uni_check = None
    if a.check_universe:
        other = {l.strip() for l in open(a.check_universe) if l.strip()}
        uni_check = dict(path=os.path.abspath(a.check_universe), n=len(other),
                         equal=(other == universe_set),
                         only_in_index=len(universe_set - other),
                         only_in_list=len(other - universe_set))
        log(f"universe cross-check vs {a.check_universe}: equal={uni_check['equal']}")
    timings["universe"] = round(time.time() - t, 1)

    # 2. metadata
    t = time.time()
    release, mfiles, meta_err = get_metadata_files(a.meta_dir, a.fallback_meta_dir, a.offline)
    meta, meta_stats = load_metadata(release, mfiles, universe)
    log(f"metadata {release}: {meta_stats}")
    meta_out = os.path.join(raw_dir, f"gtdb_{release}_reps_meta.tsv.gz")
    with open(meta_out + ".tmp", "wb") as raw, \
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        cols = ["gtdb_accession", "meta_accession", "meta_join",
                "ncbi_genbank_assembly_accession", "ncbi_taxid", "ncbi_species_taxid",
                "ncbi_genome_category", "ncbi_biosample", "genome_size",
                "gtdb_representative"]
        gz.write(("\t".join(cols) + "\n").encode())
        for acc in universe:
            m = meta[acc]
            gz.write(("\t".join([acc] + [m.get(c, "") for c in cols[1:]]) + "\n").encode())
    os.replace(meta_out + ".tmp", meta_out)
    timings["metadata"] = round(time.time() - t, 1)

    # PGFam cache keys + merged.dmp (tie-break inputs)
    t = time.time()
    pg_nonempty, pg_empty, pg_odd = scan_pgfam_cache(a.pgfam_cache)
    st = os.stat(a.pgfam_cache)
    pg_stats = dict(path=os.path.abspath(a.pgfam_cache), bytes=st.st_size,
                    mtime=datetime.datetime.fromtimestamp(
                        st.st_mtime, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                    n_keys=len(pg_nonempty) + len(pg_empty), n_empty=len(pg_empty),
                    n_non_genome_keys=pg_odd,
                    n_gca_gcf_keys=sum(1 for k in pg_nonempty | pg_empty
                                       if k.startswith(("GCA_", "GCF_"))))
    log(f"PGFam cache: {pg_stats['n_keys']} genome keys, {len(pg_empty)} empty "
        f"({time.time() - t:.0f}s)")
    merged = load_merged(a.merged_dmp)
    timings["pgfam_cache_scan_and_merged"] = round(time.time() - t, 1)

    # 3. pass 1: exact / prefix-swap / metadata GCA
    t = time.time()
    s_index = collections.defaultdict(list)   # string -> [(rep, match_type)]
    rep_strings = {}
    for acc in universe:
        strs = {acc: "exact"}
        strs.setdefault(swap_prefix(acc), "prefix_swap")
        gca = meta[acc]["ncbi_genbank_assembly_accession"]
        if gca.startswith(("GCA_", "GCF_")):
            strs.setdefault(gca, "metadata_gca")
        rep_strings[acc] = strs
        for s, mt in strs.items():
            s_index[s].append((acc, mt))
    # batch by representative (each one's strings stay together), not by sorted string:
    # a batch of only GCF swaps of GCA MAGs can legitimately be empty, a batch of ~490
    # whole representatives cannot -- which is what makes the empty-answer guard valid
    strings = list(dict.fromkeys(s for acc in universe for s in rep_strings[acc]))
    fields = ",".join(FIELDS)

    def q_in(field):
        return lambda vals, limit: (f"in({field},({','.join(vals)}))&select({fields})"
                                    f"&limit({limit})")

    gsize = {acc: to_int(meta[acc]["genome_size"]) for acc in universe}
    gtax = {acc: to_int(meta[acc]["ncbi_taxid"]) for acc in universe}
    gate_on = not a.no_length_gate
    cands = {acc: {} for acc in universe}

    def dead(x):
        return x["row"].get("genome_status") in DROP_STATUS

    def gate_fails(acc, x):
        """True when candidate x fails the genome-length gate (see the module docstring);
        a record or representative without a length cannot be tested and passes."""
        if not gate_on:
            return False
        d = rel_len_dev(x["row"], gsize[acc])
        if d is None:
            return False
        mt = x["match_type"]
        if mt in LEN_WEAK_TYPES:
            return abs(d) > a.len_tol_weak
        if mt in ("exact", "prefix_swap"):
            return (abs(d) > a.len_tol_accession and taxid_relation(
                gid_key(x["row"]["genome_id"])[0], gtax[acc], merged) == "mismatch")
        if mt == "taxid_length":
            return abs(d) > a.len_tol_taxid
        return False

    def is_usable(acc):
        """Has a non-Deprecated/Plasmid candidate that passes the length gate (the
        biosample guard needs every pass's claims, so it is applied after pass 2b)."""
        return any(not dead(x) and not gate_fails(acc, x) for x in cands[acc].values())

    def any_live(acc):
        return any(not dead(x) for x in cands[acc].values())

    def dead_lookalike(acc, x):
        """A Deprecated/Plasmid accession match that looks like this genome: its length is
        within --len-tol-weak of the GTDB genome_size (or cannot be tested). A Plasmid
        record that holds only a plasmid of the assembly does not qualify."""
        if not dead(x) or x["match_type"] not in ACCESSION_MATCHES:
            return False
        d = rel_len_dev(x["row"], gsize[acc])
        return d is None or abs(d) <= a.len_tol_weak

    p1, p1_stats = run_pass("pass1_exact", strings, a.batch, q_in("assembly_accession"),
                            raw_dir, workers, a.limit, a.pause, empty_is_error=True,
                            fields=FIELDS)
    unattributed = 0
    for rec in p1:
        for row in rec["rows"]:
            hits = s_index.get(row.get("assembly_accession"))
            if not hits:
                unattributed += 1
                continue
            for acc, mt in hits:
                add_candidate(cands[acc], row, mt, rec["queried_at"], "pass1_exact")
    p1_stats.update(n_strings=len(strings), rows_unattributed=unattributed,
                    reps_with_any_record=sum(1 for c in cands.values() if c),
                    reps_usable_before_length_gate=sum(1 for acc in universe if any_live(acc)),
                    reps_usable=sum(1 for acc in universe if is_usable(acc)))
    log(f"pass 1: {p1_stats}")
    timings["pass1"] = round(time.time() - t, 1)

    # 4a. pass 2: version-free wildcard on both prefixes
    t = time.time()
    unres = [acc for acc in universe if not is_usable(acc)]
    unres_gated = sum(1 for acc in unres if any_live(acc))
    core_reps = collections.defaultdict(list)
    for acc in unres:
        core_reps[acc_core(acc)].append(acc)
    terms = sorted(f"eq(assembly_accession,{p}_{core}.*)" for core in core_reps
                   for p in ("GCA", "GCF"))
    p2w, p2w_stats = run_pass(
        "pass2_wildcard", terms, a.wild_batch,
        lambda vals, limit: f"or({','.join(vals)})&select({fields})&limit({limit})",
        raw_dir, workers, a.limit, a.pause, empty_is_error=False, fields=FIELDS)
    off_core = 0
    for rec in p2w:
        for row in rec["rows"]:
            reps = core_reps.get(acc_core(row.get("assembly_accession")))
            if not reps:
                off_core += 1
                continue
            for acc in reps:
                mt = rep_strings[acc].get(row.get("assembly_accession"), "wildcard_core")
                add_candidate(cands[acc], row, mt, rec["queried_at"], "pass2_wildcard")
    after_w = [acc for acc in unres if not is_usable(acc)]
    p2w_stats.update(n_reps_queried=len(unres),
                     n_reps_queried_only_because_the_length_gate_failed=unres_gated,
                     rows_off_core=off_core, reps_newly_usable=len(unres) - len(after_w))
    log(f"pass 2 wildcard: {p2w_stats}")
    timings["pass2_wildcard"] = round(time.time() - t, 1)

    # 4b. pass 2: BioSample
    t = time.time()
    bios_reps = collections.defaultdict(list)
    after_w_gated = sum(1 for acc in after_w if any_live(acc))
    no_bios = 0
    for acc in after_w:
        b = meta[acc]["ncbi_biosample"]
        if re.match(r"^SAM[A-Z]+\d+$", b or ""):
            bios_reps[b].append(acc)
        else:
            no_bios += 1
    bios = sorted(bios_reps)
    p2b, p2b_stats = run_pass("pass2_biosample", bios, a.bios_batch,
                              q_in("biosample_accession"), raw_dir, workers, a.limit,
                              a.pause, empty_is_error=False, fields=FIELDS)
    bios_rows_per_rep = collections.Counter()
    for rec in p2b:
        for row in rec["rows"]:
            for acc in bios_reps.get(row.get("biosample_accession"), ()):
                ra = row.get("assembly_accession")
                if ra in rep_strings[acc]:
                    mt = rep_strings[acc][ra]
                elif ra and acc_core(ra) == acc_core(acc):
                    mt = "wildcard_core"
                else:
                    mt = "biosample"
                if add_candidate(cands[acc], row, mt, rec["queried_at"], "pass2_biosample"):
                    bios_rows_per_rep[acc] += 1
    after_b = [acc for acc in after_w if not is_usable(acc)]
    p2b_stats.update(n_reps_queried=len(after_w),
                     n_reps_queried_only_because_the_length_gate_failed=after_w_gated,
                     reps_without_biosample=no_bios,
                     reps_newly_usable_before_guard=len(after_w) - len(after_b),
                     biosample_rows_per_rep=dict(sorted(collections.Counter(
                         bios_rows_per_rep.values()).items())))
    log(f"pass 2 biosample: {p2b_stats}")
    timings["pass2_biosample"] = round(time.time() - t, 1)

    # 5. biosample guard: a biosample match is weaker evidence than an assembly accession.
    # A GTDB representative is one assembly, so a biosample-matched genome that already
    # carries ANOTHER representative's assembly accession cannot be this one's genome; and
    # a biosample-matched genome must be taxid-consistent (a metagenome BioSample can hold
    # many MAGs) unless it has the same BioSample and a length within --len-tol-rename (an
    # NCBI taxon rename). Rejected candidates stay in `candidates`; --no-biosample-guard
    # disables.
    t = time.time()
    claims = collections.defaultdict(set)   # genome_id -> reps whose accession it carries
    for acc in universe:
        for gid, x in cands[acc].items():
            if x["match_type"] != "biosample":
                claims[gid].add(acc)
    rejected, guard = {}, collections.Counter()   # acc -> {genome_id: reason}
    rename_kept = []
    if not a.no_biosample_guard:
        for acc in universe:
            rej = {}
            for gid, x in cands[acc].items():
                if x["match_type"] != "biosample" or dead(x):
                    continue
                if claims.get(gid, set()) - {acc}:
                    rej[gid] = "claimed_by_another_rep"
                    guard["candidates_rejected_claimed_by_another_rep"] += 1
                elif taxid_relation(gid_key(gid)[0], gtax[acc], merged) not in ("exact",
                                                                                "merged"):
                    d = rel_len_dev(x["row"], gsize[acc])
                    bs = meta[acc]["ncbi_biosample"]
                    if (bs and x["row"].get("biosample_accession") == bs and d is not None
                            and abs(d) <= a.len_tol_rename):
                        guard["candidates_kept_taxid_differs_same_biosample_length_within_tol"] += 1
                        rename_kept.append(dict(
                            acc=acc, genome_id=gid, gtdb_ncbi_taxid=meta[acc]["ncbi_taxid"],
                            bvbrc_taxon_id=x["row"].get("taxon_id"),
                            genome_name=x["row"].get("genome_name"),
                            genome_length=x["row"].get("genome_length"),
                            gtdb_genome_size=gsize[acc], length_dev_pct=round(100 * d, 3)))
                    else:
                        rej[gid] = "taxid_mismatch"
                        guard["candidates_rejected_taxid_mismatch"] += 1
            if rej:
                rejected[acc] = rej
    guard_only = {acc: dict(r) for acc, r in rejected.items()}

    # 5b. length gate over every non-Deprecated/Plasmid candidate of every representative
    gate = {k: collections.Counter() for k in ("checked", "missing_length", "rejected",
                                               "checked_taxid_mismatch")}
    gate_examples = []
    for acc in universe:
        for gid, x in cands[acc].items():
            if dead(x):
                continue
            mt = x["match_type"]
            d = rel_len_dev(x["row"], gsize[acc])
            gate["missing_length" if d is None else "checked"][mt] += 1
            if d is not None and mt in ("exact", "prefix_swap") and taxid_relation(
                    gid_key(gid)[0], gtax[acc], merged) == "mismatch":
                gate["checked_taxid_mismatch"][mt] += 1
            if gate_fails(acc, x):
                gate["rejected"][mt] += 1
                rejected.setdefault(acc, {}).setdefault(gid, "length")
                if len(gate_examples) < 40:
                    gate_examples.append(dict(
                        acc=acc, genome_id=gid, match_type=mt,
                        bvbrc_assembly_accession=x["row"].get("assembly_accession") or "",
                        genome_length=x["row"].get("genome_length"),
                        gtdb_genome_size=gsize[acc], length_dev_pct=round(100 * d, 2)))

    # 5c. choose over passes 1-2 (before pass 3, so pass 3 can exclude what is taken)
    chosen12, pre_guard_gid, pre_gate_gid = {}, {}, {}
    for acc in universe:
        c = cands[acc]
        ch = choose(c, gtax[acc], merged, pg_nonempty, rejected.get(acc))[0]
        if ch:
            chosen12[acc] = ch["row"]["genome_id"]
        if acc in rejected:
            b0 = choose(c, gtax[acc], merged, pg_nonempty)[0]
            b1 = choose(c, gtax[acc], merged, pg_nonempty, guard_only.get(acc))[0]
        else:
            b0 = b1 = ch
        pre_guard_gid[acc] = b0["row"]["genome_id"] if b0 else ""
        pre_gate_gid[acc] = b1["row"]["genome_id"] if b1 else ""
    timings["guard_and_gate"] = round(time.time() - t, 1)

    # 5d. pass 3: GTDB ncbi_taxid + genome length within --len-tol-taxid, for what is left
    t = time.time()
    unres3 = [acc for acc in universe if acc not in chosen12]
    p3_stats = dict(enabled=not a.no_pass3, n_reps_unresolved_after_pass2=len(unres3))
    p3_reasons, p3_examples, p3_accepted = collections.Counter(), [], {}
    if not a.no_pass3:
        term_reps = collections.defaultdict(list)
        tax_reps = collections.defaultdict(list)
        skipped = collections.Counter()
        for acc in unres3:
            T, S = gtax[acc], gsize[acc]
            if not T or not S:
                skipped["no_gtdb_ncbi_taxid" if not T else "no_gtdb_genome_size"] += 1
                continue
            lo = math.ceil(S * (1 - a.len_tol_taxid))
            hi = math.floor(S * (1 + a.len_tol_taxid))
            for tx in dict.fromkeys([T, merged.get(T, T)]):
                term_reps[f"and(eq(taxon_id,{tx}),ge(genome_length,{lo}),"
                          f"le(genome_length,{hi}))"].append(acc)
                tax_reps[tx].append(acc)
        p3, p3q = run_pass(
            "pass3_taxid_length", sorted(term_reps), a.p3_batch,
            lambda vals, limit: f"or({','.join(vals)})&select({fields})&limit({limit})",
            raw_dir, workers, a.limit, a.pause, empty_is_error=False, fields=FIELDS,
            single_truncation_ok=True)
        truncated_reps = {acc for rec in p3 for v in (rec.get("truncated_values") or [])
                          for acc in term_reps.get(v, ())}
        found = collections.defaultdict(dict)   # acc -> {genome_id: (row, queried_at)}
        for rec in p3:
            for row in rec["rows"]:
                gid = row.get("genome_id")
                if not gid or not GID_RE.match(str(gid)):
                    continue
                for acc in dict.fromkeys(tax_reps.get(to_int(row.get("taxon_id")), ())):
                    d = rel_len_dev(row, gsize[acc])
                    if d is not None and abs(d) <= a.len_tol_taxid:
                        found[acc][gid] = (row, rec["queried_at"])
        taken = set(chosen12.values())
        p3_rej = collections.defaultdict(dict)
        surv = {}
        for acc in unres3:
            s = []
            n_live = sum(1 for row, _ in found.get(acc, {}).values()
                         if row.get("genome_status") not in DROP_STATUS)
            for gid, (row, _) in found.get(acc, {}).items():
                if row.get("genome_status") in DROP_STATUS or gid in cands[acc]:
                    continue   # Deprecated/Plasmid, or already judged in passes 1-2
                ra = row.get("assembly_accession") or ""
                if gid in taken or claims.get(gid, set()) - {acc}:
                    p3_rej[acc][gid] = "taxid_length_taken"
                elif ra and acc_core(ra) != acc_core(acc):
                    p3_rej[acc][gid] = "taxid_length_other_assembly"
                else:
                    s.append(gid)
            if found.get(acc):
                surv[acc] = (s, n_live)
        n_surv = collections.Counter(g for s, _ in surv.values() for g in s)
        p3_window = collections.Counter()
        for acc, (s, n_live) in surv.items():
            p3_window[min(n_live, 5)] += 1
            for g in s:
                if n_surv[g] > 1:
                    p3_rej[acc][g] = "taxid_length_shared"
            s = [g for g in s if n_surv[g] == 1]
            if len(s) == 1 and n_live == 1 and acc not in truncated_reps:
                p3_accepted[acc] = s[0]
            else:
                for g in s:
                    p3_rej[acc][g] = "taxid_length_ambiguous"
        # only the accepted record joins `candidates`: the others are records of OTHER
        # genomes that share the taxid and length window (they are counted here, listed
        # in stats.pass3_taxid_length.examples and kept in raw/pass3_taxid_length)
        for acc in found:
            for gid, why in p3_rej.get(acc, {}).items():
                p3_reasons[why] += 1
        for acc, gid in p3_accepted.items():
            row, qa = found[acc][gid]
            add_candidate(cands[acc], row, "taxid_length", qa, "pass3_taxid_length")
        for acc in unres3:
            if acc in found and len(p3_examples) < 40:
                p3_examples.append(dict(
                    acc=acc, gtdb_ncbi_taxid=meta[acc]["ncbi_taxid"],
                    gtdb_genome_size=gsize[acc], accepted=p3_accepted.get(acc, ""),
                    found={g: dict(status=r.get("genome_status"),
                                   assembly_accession=r.get("assembly_accession") or "",
                                   genome_length=r.get("genome_length"),
                                   name=r.get("genome_name"),
                                   reason=p3_rej.get(acc, {}).get(g, "accepted" if
                                                                  p3_accepted.get(acc) == g
                                                                  else "deprecated_or_known"))
                           for g, (r, _) in sorted(found[acc].items(),
                                                   key=lambda kv: gid_key(kv[0]))[:12]}))
        p3_stats.update(p3q, rule=(
            f"or(and(eq(taxon_id,T),ge(genome_length,lo),le(genome_length,hi))) on the GTDB "
            f"ncbi_taxid and its merged.dmp target, window +-{a.len_tol_taxid} of genome_size; "
            f"accept only when the window holds exactly one live record, not chosen or "
            f"accession-claimed by another representative, with no different assembly "
            f"accession"),
                        live_records_in_window_per_rep_with_any=dict(
                            sorted(((f"{k}+" if k == 5 else str(k)), v)
                                   for k, v in p3_window.items())),
                        n_terms=len(term_reps), reps_skipped=dict(skipped),
                        n_reps_queried=len({x for v in term_reps.values() for x in v}),
                        reps_with_any_record=len(found),
                        reps_truncated=len(truncated_reps),
                        reps_accepted=len(p3_accepted),
                        candidates_rejected=dict(p3_reasons.most_common()))
        log(f"pass 3 taxid+length: {p3_stats}")
    timings["pass3_taxid_length"] = round(time.time() - t, 1)

    # 6. choose + write
    t = time.time()
    out_rows = {}
    guard_examples = []
    bios_cons = collections.defaultdict(collections.Counter)
    decisive = collections.Counter()
    multi_examples = []
    usable_hist, all_hist = collections.Counter(), collections.Counter()
    usable_total = 0
    tax_prefix, tax_bvbrc, tax_species = (collections.Counter(), collections.Counter(),
                                          collections.Counter())
    mismatch_by_mt = collections.Counter()
    bios_diag = collections.Counter()
    len_bins = collections.defaultdict(collections.Counter)
    crossed, crossed_kind = collections.Counter(), collections.Counter()
    crossed_examples = []
    gate_changed = []
    for acc in universe:
        m = meta[acc]
        gtaxid = gtax[acc]
        c = cands[acc]
        rej = rejected.get(acc, {})
        chosen, ordered, step = choose(c, gtaxid, merged, pg_nonempty, rej)
        g_rej = guard_only.get(acc)
        if g_rej:
            guard["rows_with_a_rejected_candidate"] += 1
            guard["rows_left_without_usable_record" if not pre_gate_gid[acc] else
                  "rows_still_resolved_by_another_candidate"] += 1
            if len(guard_examples) < 25:
                guard_examples.append(dict(acc=acc, gtdb_ncbi_taxid=m["ncbi_taxid"],
                                           rejected=sorted(g_rej, key=gid_key),
                                           chosen_after_guard=pre_gate_gid[acc],
                                           chosen_final=chosen["row"]["genome_id"]
                                           if chosen else "",
                                           claimed_by={g: sorted(claims.get(g, set()) - {acc})
                                                       for g in g_rej}))
        if pre_gate_gid[acc] != chosen12.get(acc, ""):
            gate_changed.append(acc)
        n_usable = sum(1 for x in c.values() if not dead(x)
                       and x["row"]["genome_id"] not in rej)
        if c:
            all_hist[min(len(c), 6)] += 1
        rej_listed = [x["row"]["genome_id"] for x in ordered
                      if not dead(x) and x["row"]["genome_id"] in rej]
        row = dict(gtdb_accession=acc, gtdb_genbank_accession=m["ncbi_genbank_assembly_accession"],
                   acc_core=acc_core(acc), gtdb_ncbi_taxid=m["ncbi_taxid"],
                   ncbi_genome_category=m["ncbi_genome_category"],
                   n_candidates=len(c),
                   candidates=";".join(f"{x['row']['genome_id']}:"
                                       f"{x['row'].get('genome_status') or 'NA'}"
                                       for x in ordered),
                   gtdb_genome_size=m["genome_size"],
                   rejected_candidates=";".join(f"{g}:{rej[g]}" for g in rej_listed))
        if chosen is None:
            reasons = {rej[g] for g in rej_listed}
            if any(dead_lookalike(acc, x) for x in c.values()):
                label = "deprecated_only"
            elif "length" in reasons:
                label = "length_rejected"
            elif reasons & set(GUARD_REASONS):
                label = "biosample_rejected"
            elif any(dead(x) and x["match_type"] != "taxid_length" for x in c.values()):
                label = "deprecated_only"
            else:
                label = "unmapped"
            row.update(bvbrc_genome_id="", bvbrc_assembly_accession="", match_type=label,
                       genome_status="", genome_quality="", bvbrc_taxon_id="",
                       taxid_consistent="", cds="", in_pgfam_cache="", queried_at="",
                       bvbrc_genome_length="", length_ratio="",
                       bvbrc_label_taxid_consistent="")
        else:
            r = chosen["row"]
            gid = r["genome_id"]
            rel = taxid_relation(gid_key(gid)[0], gtaxid, merged)
            lrel = taxid_relation(to_int(r.get("taxon_id")), gtaxid, merged)
            gl, gs = to_int(r.get("genome_length")), gsize[acc]
            d = rel_len_dev(r, gs)
            row.update(bvbrc_genome_id=gid, bvbrc_assembly_accession=r.get("assembly_accession")
                       or "", match_type=chosen["match_type"],
                       genome_status=r.get("genome_status") or "",
                       genome_quality=r.get("genome_quality") or "",
                       bvbrc_taxon_id="" if r.get("taxon_id") is None else str(r["taxon_id"]),
                       taxid_consistent=rel,
                       cds="" if r.get("cds") is None else str(r["cds"]),
                       in_pgfam_cache=("nonempty" if gid in pg_nonempty else
                                       "empty" if gid in pg_empty else "absent"),
                       queried_at=chosen["queried_at"],
                       bvbrc_genome_length="" if gl is None else str(gl),
                       length_ratio="" if d is None else f"{gl / gs:.4f}",
                       bvbrc_label_taxid_consistent=lrel)
            len_bins[chosen["match_type"]][len_bin(d)] += 1
            rb, gb = r.get("biosample_accession"), m["ncbi_biosample"]
            if (chosen["match_type"] in ACCESSION_MATCHES and d is not None and abs(d) <= 0.01
                    and (lrel == "mismatch" or rel == "mismatch")):
                rbn = rb if rb not in (None, "", "-") else None   # BV-BRC writes '-' for none
                kind = ("crossed_biosample_differs" if rbn and gb and rbn != gb else
                        "same_biosample_relabelled" if rbn and gb else "biosample_missing")
                crossed[f"{kind}__label_{lrel}__prefix_{rel}"] += 1
                crossed_kind[kind] += 1
                if (len(crossed_examples) < 40 and
                        (kind != "same_biosample_relabelled" or crossed_kind[kind] <= 10)):
                    crossed_examples.append(dict(kind=kind,
                            acc=acc, genome_id=gid, match_type=chosen["match_type"],
                            gtdb_ncbi_taxid=m["ncbi_taxid"], bvbrc_taxon_id=r.get("taxon_id"),
                            genome_name=r.get("genome_name"),
                            gtdb_biosample=m["ncbi_biosample"],
                            bvbrc_biosample=r.get("biosample_accession"),
                            length_ratio=row["length_ratio"]))
            usable_hist[min(n_usable, 6)] += 1
            usable_total += n_usable
            decisive[step or "single_candidate"] += 1
            if n_usable > 1 and len(multi_examples) < 25:
                multi_examples.append(dict(acc=acc, chosen=gid, decided_by=step,
                                           candidates=row["candidates"]))
            tax_prefix[rel] += 1
            if rel == "mismatch":
                mismatch_by_mt[chosen["match_type"]] += 1
                if gid_key(gid)[0] == to_int(m["ncbi_species_taxid"]):
                    tax_species["prefix_equals_gtdb_species_taxid"] += 1
            tax_bvbrc[lrel] += 1
            bios_cons[chosen["match_type"]]["missing" if not rb or not gb else
                                            "equal" if rb == gb else "differ"] += 1
            if chosen["match_type"] == "biosample":
                ra = r.get("assembly_accession")
                bios_diag["chosen_has_no_assembly_accession" if not ra else
                          "chosen_is_a_different_assembly"] += 1
                bios_diag[f"taxid_{rel}"] += 1
                prefixes = {gid_key(x["row"]["genome_id"])[0] for x in c.values()
                            if x["match_type"] == "biosample"}
                if len(prefixes) > 1:
                    bios_diag["candidates_span_several_taxids"] += 1
        out_rows[acc] = row
    shared_final = collections.Counter(r["bvbrc_genome_id"] for r in out_rows.values()
                                       if r["bvbrc_genome_id"])
    shared_final = {g: k for g, k in shared_final.items() if k > 1}
    if shared_final:
        log(f"WARNING: {len(shared_final)} genome_ids are chosen by more than one "
            f"representative: {list(shared_final)[:10]}")

    out_tsv = os.path.join(a.outdir, "gtdb_r226_acc_to_bvbrc.tsv.gz")
    with open(out_tsv + ".tmp", "wb") as raw, \
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        gz.write(("\t".join(COLUMNS) + "\n").encode())
        for acc in universe:
            r = out_rows[acc]
            vals = [str(r[cname]).replace("\t", " ").replace("\n", " ") for cname in COLUMNS]
            gz.write(("\t".join(vals) + "\n").encode())
    os.replace(out_tsv + ".tmp", out_tsv)
    timings["choose_and_write"] = round(time.time() - t, 1)

    # 6. stats
    n = len(universe)
    resolved = [r for r in out_rows.values() if r["bvbrc_genome_id"]]
    def many_to_one(gid_of):
        reps_of = collections.defaultdict(list)
        for acc_, g in gid_of.items():
            if g:
                reps_of[g].append(acc_)
        shared_ = {g: v for g, v in reps_of.items() if len(v) > 1}
        return dict(n_genome_ids_chosen_by_more_than_one_rep=len(shared_),
                    n_reps_involved=sum(len(v) for v in shared_.values()),
                    cases={g: [dict(acc=x, match_type=out_rows[x]["match_type"],
                                    gtdb_ncbi_taxid=out_rows[x]["gtdb_ncbi_taxid"])
                               for x in v] for g, v in list(shared_.items())[:30]})
    by_cat = collections.defaultdict(list)
    by_prefix = collections.defaultdict(list)
    for acc in universe:
        cat = meta[acc]["ncbi_genome_category"] or "missing"
        by_cat[{"none": "isolate (none)"}.get(cat, cat)].append(acc)
        by_prefix[acc[:3]].append(acc)
    report_sets = {}
    for spec in a.report_set:
        name, path = spec.split("=", 1)
        report_sets[name] = dict(path=os.path.abspath(path),
                                 **group_coverage(open(path).read().split(), out_rows,
                                                  universe_set))
    set_accs = {}
    for spec in a.report_set:
        name, path = spec.split("=", 1)
        set_accs[name] = open(path).read().split()
    for spec in a.report_union:
        name, parts = spec.split("=", 1)
        accs = sorted({x for part in parts.split("+") for x in set_accs[part]})
        report_sets[name] = dict(union_of=parts.split("+"),
                                 **group_coverage(accs, out_rows, universe_set))
    mag_sel = {}
    for spec in a.report_mag_selection:
        name, path = spec.split("=", 1)
        mag_sel[name] = mag_selection_coverage(path, out_rows)
    acc_dev = [abs(int(r["bvbrc_genome_length"]) / int(r["gtdb_genome_size"]) - 1)
               for r in resolved
               if r["match_type"] in ("exact", "prefix_swap") and r["bvbrc_genome_length"]
               and r["gtdb_genome_size"] and r["taxid_consistent"] != "mismatch"]
    acc_off_n = sum(1 for x in acc_dev if x > a.len_tol_accession)
    acc_off_weak_n = sum(1 for x in acc_dev if a.len_tol_weak < x <= a.len_tol_accession)
    timings["total"] = round(time.time() - t_start, 1)
    stats = dict(
        built_at=utcnow(), script=os.path.abspath(__file__), argv=sys.argv[1:],
        user_agent=UA,
        outputs=dict(map_tsv=os.path.abspath(out_tsv), raw_dir=os.path.abspath(raw_dir),
                     reps_meta=os.path.abspath(meta_out)),
        universe=dict(index_db=os.path.abspath(a.index_db), n=n, n_occurrences=n_occ,
                      n_gca=len(by_prefix["GCA"]), n_gcf=len(by_prefix["GCF"]),
                      cross_check=uni_check),
        metadata=dict(release_used=release, r226_error=meta_err, files=mfiles,
                      fallback_dir=os.path.abspath(a.fallback_meta_dir),
                      fallback_used=release != "r226", **meta_stats),
        pgfam_cache=pg_stats,
        merged_dmp=dict(path=os.path.abspath(a.merged_dmp), n=len(merged)),
        passes=dict(pass1_exact=p1_stats, pass2_wildcard=p2w_stats,
                    pass2_biosample=p2b_stats, pass3_taxid_length=p3_stats),
        coverage=dict(
            overall=dict(resolved=len(resolved), n=n, pct=pct(len(resolved), n),
                         wilson95=wilson(len(resolved), n),
                         any_record=sum(1 for r in out_rows.values() if r["n_candidates"]),
                         resolved_excluding_taxid_length=sum(
                             1 for r in resolved if r["match_type"] != "taxid_length"),
                         **{lab: sum(1 for r in out_rows.values() if r["match_type"] == lab)
                            for lab in UNRESOLVED_TYPES}),
            by_genome_category={k: group_coverage(v, out_rows, universe_set)
                                for k, v in sorted(by_cat.items())},
            by_prefix={k: group_coverage(v, out_rows, universe_set)
                       for k, v in sorted(by_prefix.items())},
            report_sets=report_sets, mag_selections=mag_sel),
        match_type=dict(collections.Counter(r["match_type"]
                                            for r in out_rows.values()).most_common()),
        one_to_many=dict(
            records_per_rep_with_any_record=dict(sorted(
                ((f"{k}+" if k == 6 else str(k)), v) for k, v in all_hist.items())),
            usable_records_per_resolved_rep=dict(sorted(
                ((f"{k}+" if k == 6 else str(k)), v) for k, v in usable_hist.items())),
            mean_records_per_rep_with_any=round(
                sum(r["n_candidates"] for r in out_rows.values()) /
                max(1, sum(1 for r in out_rows.values() if r["n_candidates"])), 3),
            mean_usable_per_resolved=round(usable_total / max(1, len(resolved)), 3),
            reps_with_multiple_usable=sum(v for k, v in usable_hist.items() if k > 1),
            decided_by=dict(decisive.most_common()),
            examples=multi_examples,
            many_to_one=many_to_one({x: r["bvbrc_genome_id"] for x, r in out_rows.items()}),
            many_to_one_before_biosample_guard_and_length_gate=many_to_one(pre_guard_gid),
            many_to_one_after_guard_before_length_gate=many_to_one(pre_gate_gid)),
        biosample_guard=dict(enabled=not a.no_biosample_guard,
                             taxid_rename_exemption=dict(
                                 rule="a biosample match with a taxid mismatch is kept when its "
                                      "BioSample equals GTDB's and |genome_length - genome_size|"
                                      f" / genome_size <= {a.len_tol_rename}",
                                 kept=rename_kept),
                             rows_resolved_by_biosample_final=sum(
                                 1 for r in out_rows.values() if r["match_type"] == "biosample"),
                             **dict(guard.most_common()),
                             examples=guard_examples),
        length_gate=dict(
            enabled=gate_on,
            rule=(f"reject a non-Deprecated/Plasmid candidate when |genome_length - genome_size|"
                  f" / genome_size > {a.len_tol_weak} for {'/'.join(LEN_WEAK_TYPES)}, or > "
                  f"{a.len_tol_accession} for exact/prefix_swap whose genome_id taxid prefix is "
                  f"a mismatch; taxid_length candidates must be within {a.len_tol_taxid}; a "
                  f"candidate without a length is kept"),
            len_tol_weak=a.len_tol_weak, len_tol_accession=a.len_tol_accession,
            len_tol_taxid=a.len_tol_taxid,
            candidates_checked_by_match_type=dict(gate["checked"].most_common()),
            candidates_checked_exact_prefix_swap_with_taxid_mismatch=dict(
                gate["checked_taxid_mismatch"].most_common()),
            candidates_missing_length_by_match_type=dict(gate["missing_length"].most_common()),
            candidates_failing_gate_by_match_type_incl_guard_rejected=dict(
                gate["rejected"].most_common()),
            reps_whose_pass1_2_choice_the_gate_changed=len(gate_changed),
            reps_left_length_rejected=sum(1 for r in out_rows.values()
                                          if r["match_type"] == "length_rejected"),
            changed=[dict(acc=x, after_guard=pre_gate_gid[x], after_gate=chosen12.get(x, ""),
                          final=out_rows[x]["bvbrc_genome_id"],
                          final_match_type=out_rows[x]["match_type"])
                     for x in gate_changed[:200]],
            examples=gate_examples),
        pass3_taxid_length=dict(p3_stats, accepted=dict(sorted(p3_accepted.items())),
                                examples=p3_examples),
        chosen_length_vs_gtdb_genome_size={k: dict(sorted(v.items()))
                                           for k, v in sorted(len_bins.items())},
        chosen_bvbrc_biosample_vs_gtdb_biosample={k: dict(v.most_common())
                                                  for k, v in sorted(bios_cons.items())},
        taxid=dict(
            genome_id_prefix_vs_gtdb_ncbi_taxid=dict(tax_prefix.most_common()),
            bvbrc_taxon_id_vs_gtdb_ncbi_taxid=dict(tax_bvbrc.most_common()),
            bvbrc_label_taxid_consistent=dict(tax_bvbrc.most_common()),
            label_mismatch_on_sequence_correct_records=dict(
                note="chosen accession matches whose genome_length is within 1% of GTDB's "
                     "genome_size (the right sequence) but whose BV-BRC taxon_id label and/or "
                     "genome_id prefix disagree with gtdb_ncbi_taxid. "
                     "crossed_biosample_differs: the record also carries another BioSample "
                     "(crossed-label BV-BRC records, e.g. GCF_002020805.1 -> 1918949.3); "
                     "same_biosample_relabelled: same BioSample, the taxon was reassigned "
                     "after BV-BRC took the record. Either way integration must take taxon_id "
                     "from gtdb_ncbi_taxid, not from bvbrc_taxon_id or the genome_id prefix.",
                by_kind=dict(crossed_kind.most_common()),
                counts=dict(crossed.most_common()), n=sum(crossed.values()),
                examples=crossed_examples),
            prefix_mismatch_by_match_type=dict(mismatch_by_mt.most_common()),
            biosample_matches=dict(bios_diag.most_common()),
            **dict(tax_species)),
        chosen_in_pgfam_cache=dict(collections.Counter(r["in_pgfam_cache"]
                                                       for r in resolved).most_common()),
        shared_genome_ids_final=len(shared_final),
        runtime_secs=timings,
        integration_notes=INTEGRATION_NOTES + [
            f"{acc_off_n} chosen exact/prefix_swap records whose taxid agrees differ from the "
            f"GTDB genome_size by more than {a.len_tol_accession:.0%} "
            f"({acc_off_weak_n} more by {a.len_tol_weak:.0%}-{a.len_tol_accession:.0%}); the "
            f"length gate tests accession matches only when their taxid disagrees, so filter on "
            f"length_ratio where a same-sequence guarantee matters."],
        integration_todo=INTEGRATION_TODO,
    )
    out_stats = os.path.join(a.outdir, "gtdb_r226_acc_to_bvbrc.stats.json")
    with open(out_stats + ".tmp", "w") as fh:
        json.dump(stats, fh, indent=1)
    os.replace(out_stats + ".tmp", out_stats)
    log(f"resolved {len(resolved)}/{n} ({pct(len(resolved), n)}%); "
        f"match types {stats['match_type']}; total {timings['total']}s")
    log(f"wrote {out_tsv} and {out_stats}")


if __name__ == "__main__":
    main()
