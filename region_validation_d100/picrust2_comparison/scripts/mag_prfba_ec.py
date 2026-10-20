#!/usr/bin/env python
"""Re-derive the prFBA MAG arm's selected reference SETS, so it can be scored in EC space.

``mag_validation/out/mag_cells_all.csv`` keeps only the anchor and the set size, and
its gene capture is in PGFam space; the PICRUSt2 benchmark is scored on ECs. This
re-runs the MAG arm's own selection path -- ``score_mag.build_record`` /
``exclude_hits`` / ``select_references.select_representatives`` with ``ANI_KNOBS``,
tau_ANI 95 -- on the same skani hits (``hits/`` = the RESULTS_v2 basis, depth 100, and
``lowc_hits/``) for every uncontaminated MAG, and checks each re-derived selection
against the scorer's recorded anchor and set size.

Writes ``data/mag_arm/prfba_mag_selections.csv.gz``: mag, genome_id, completeness,
regime, picks (';'-joined), and ``matches_recorded``.
"""
import collections
import csv
import glob
import gzip
import json
import os
import sys
from multiprocessing import Pool

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import mag_common as M  # noqa: E402

sys.path.insert(0, M.MV)
import score_mag as SM  # noqa: E402

PGF = None


def work(job):
    hits_path, mani_path, cells_path = job
    meta = {r["mag"]: r for r in json.load(open(mani_path))["rows"]
            if r["contamination_target"] == 0}
    hits = collections.defaultdict(list)
    with open(hits_path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            q = os.path.basename(row["Query_file"])[:-4]
            if q in meta:
                hits[q].append(dict(ref=os.path.basename(row["Ref_file"])[:-4],
                                    ani=float(row["ANI"]),
                                    af_query=float(row["Align_fraction_query"]),
                                    af_ref=float(row["Align_fraction_ref"])))
    recorded = {}
    with open(cells_path) as fh:
        for row in csv.DictReader(fh):
            recorded[(row["mag"], row["regime"])] = (row["anchor"], int(row["n_selected"]))
    out = []
    for mag, m in sorted(meta.items()):
        gid = m["genome_id"]
        allh = hits.get(mag, [])
        for regime, hh in (("include", allh),
                           ("exclude", SM.exclude_hits(allh, gid, SM.src_taxon(gid)))):
            rec, _ = SM.build_record(mag, hh, m["src_bp"] + m["contam_bp"])
            try:
                sel = SM.SR.select_representatives(
                    rec, knobs=SM.ANI_KNOBS, gene_provider=lambda g: set(PGF.get(g) or ()))
            except Exception:
                sel = {"selected": []}
            picks = [s["genome_id"] for s in (sel.get("selected") or [])]
            got = (picks[0] if picks else "", len(picks))
            out.append(dict(mag=mag, genome_id=gid, completeness=m["completeness_target"],
                            regime=regime, picks=";".join(picks),
                            matches_recorded=int(recorded.get((mag, regime)) == got)))
    return out


def main():
    global PGF
    PGF = SM.load_pgfams()
    jobs = []
    for hdir, mdir, cdir in (("hits", "manifests", "out_batches"),
                             ("lowc_hits", "lowc_manifests", "lowc_out_batches")):
        for hp in sorted(glob.glob(os.path.join(M.MV, hdir, "*.tsv"))):
            b = os.path.basename(hp)[:-4]
            jobs.append((hp, os.path.join(M.MV, mdir, b + ".json"),
                         os.path.join(M.MV, cdir, b, "mag_cells.csv")))
    rows = []
    with Pool(len(jobs)) as pool:
        for part in pool.imap_unordered(work, jobs):
            rows += part
            print(f"[prfba_ec] {len(rows)} selections", flush=True)
    rows.sort(key=lambda r: (r["mag"], r["regime"]))
    path = os.path.join(M.DATA, "prfba_mag_selections.csv.gz")
    with gzip.open(path, "wt", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    ok = sum(r["matches_recorded"] for r in rows)
    print(f"[prfba_ec] {len(rows)} selections; {ok} match the scorer's recorded anchor + set "
          f"size ({len(rows) - ok} do not)")


if __name__ == "__main__":
    sys.exit(main())
