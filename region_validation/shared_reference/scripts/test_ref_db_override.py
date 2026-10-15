#!/usr/bin/env python
"""
test_ref_db_override.py — self-test for component 0 of the shared-reference k-fold
benchmark: the `--ref-db` / `--ref-id-db` override and the `--gpu-cand >= --k2` guard
added to align_hits.py + edlib_biopython_hits.py (SPEC.md D14, §6.2).

WHY THIS EXISTS. Every other component assumes three things that are invisible from the
outside and silently wrong if they break:
  1. the aligner really aligns against ONLY the md5->seq JSON handed to `--ref-db`
     (if it ever falls back to the 112k-genome default, the fold's held-out genome is
     back in the reference and the whole cross-validation is void);
  2. `asv_top20_alignment_hits.json` keeps its exact schema and its exact `identity`
     definition (`n_matches / aligned_len`, local SW), because that number is THE identity
     used by the fold plan's edges, the novelty covariate and select_references' tiers;
  3. the SYNTHETIC `md5 -> header` file that plan_folds.py writes round-trips through
     `edlib_biopython_hits.parse_header` into (organism, genome_id, taxon_id, feature_id).
It builds a 8-genome fixture from the real panel + the real 16S DB, runs the aligner on
both backends with the full DB replaced, then with one genome REMOVED from the replacement
DB, and asserts all three. The removal case is smoke-gate item 7 in miniature: the removed
genome's own amplicon must lose its identity-1.000 self hit.

    ~/Documents/py_venv/bin/python test_ref_db_override.py                  # cpu + cuda
    ~/Documents/py_venv/bin/python test_ref_db_override.py --backends cpu   # no GPU needed
    ~/Documents/py_venv/bin/python test_ref_db_override.py --force          # rebuild fixture

The fixture step loads the 519 MB md5->seq DB (~40 s) and is cached in --workdir, so
re-runs are fast; pass --force to rebuild. Exit 0 = every assertion held.
Interpreter: ~/Documents/py_venv/bin/python
"""
from __future__ import annotations
import argparse, json, os, subprocess, sys

PRFBA = "/home/freiburger/Documents/prFBA"
sys.path.insert(0, PRFBA)
import edlib_biopython_hits as P                      # noqa: E402

PY = os.path.expanduser("~/Documents/py_venv/bin/python")
SELECTION = f"{PRFBA}/region_validation/data/selection_10k.json"
FANOUT = f"{PRFBA}/region_validation/data/amplicon_md5_fanout.json"
AMPLICONS = f"{PRFBA}/region_validation/data/combined_amplicons.unique.fasta"

ENTRY_KEYS = {"asv_len", "midas_taxonomy", "rel_ab", "best_align_score", "best_identity",
              "search", "n_edlib_candidates", "top20"}
HIT_KEYS = {"rank", "align_score", "identity", "n_matches", "aligned_len", "edlib_distance",
            "edlib_identity", "organism", "genome_id", "taxon_id", "feature_id", "md5",
            "ref_seq_len", "ref_aln_start", "ref_aln_end", "lineage"}
PARAM_KEYS = ["prefilter_k", "k2", "min_ref_len", "workers", "gpu_cand", "limit", "fasta",
              "ref_db", "ref_id_db"]
RANKS = {"Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"}


def ok(label, detail=""):
    print(f"  [ok]   {label}{('  — ' + detail) if detail else ''}", flush=True)


# ------------------------------------------------------------------------- fixture
def build_fixture(work: str, n_genomes: int, region: str, force: bool) -> dict:
    meta_p = os.path.join(work, "fixture.json")
    if os.path.exists(meta_p) and not force:
        print(f"[fixture] reusing {meta_p} (use --force to rebuild)", flush=True)
        return json.load(open(meta_p))
    os.makedirs(work, exist_ok=True)

    genomes = json.load(open(SELECTION))["genomes"]
    g2amp: dict[str, dict[str, str]] = {}
    for md5, cells in json.load(open(FANOUT)).items():
        for cell in cells:
            gid, reg = cell.rsplit("__", 1)
            g2amp.setdefault(gid, {})[reg] = md5

    picked = [g for g in genomes
              if 1200 <= g["best_16s_len"] <= 1800 and region in g2amp.get(g["genome_id"], {})]
    picked = picked[:n_genomes]
    if len(picked) < n_genomes:
        raise SystemExit(f"[abort] only {len(picked)} eligible fixture genomes found")
    gids = [g["genome_id"] for g in picked]

    print(f"[fixture] loading {os.path.basename(P.DB_MD5_SEQ)} (519 MB) ...", flush=True)
    md5_seq = json.load(open(P.DB_MD5_SEQ))
    ref_db = {g["best_16s_md5"]: md5_seq[g["best_16s_md5"]].upper() for g in picked}
    del md5_seq
    if len(ref_db) != len(picked):
        raise SystemExit("[abort] md5 collision inside the fixture panel")

    # exactly the synthetic header format plan_folds.py writes (SPEC §4.2)
    ref_id_db = {g["best_16s_md5"]:
                 f"fig|{g['genome_id']}.rna.1|   16S ribosomal RNA   "
                 f"[{g['species']} | {g['genome_id']}]" for g in picked}
    for g in picked:
        got = P.parse_header(ref_id_db[g["best_16s_md5"]])
        want = (g["species"], g["genome_id"], int(g["genome_id"].split(".")[0]),
                f"{g['genome_id']}.rna.1")
        if got != want:
            raise SystemExit(f"[abort] parse_header round trip failed: {got} != {want}")
    ok("synthetic md5->header round-trips through parse_header", f"{len(picked)} entries")

    query_md5_of = {g["genome_id"]: g2amp[g["genome_id"]][region] for g in picked}
    want_md5 = set(query_md5_of.values())
    seqs, cur = {}, None
    with open(AMPLICONS) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line[0] == ">":
                cur = line[1:].split()[0]
            elif cur in want_md5:
                seqs[cur] = seqs.get(cur, "") + line
    if set(seqs) != want_md5:
        raise SystemExit("[abort] missing amplicon sequences for the fixture panel")

    drop_md5 = picked[0]["best_16s_md5"]
    json.dump(ref_db, open(f"{work}/ref_db.json", "w"))
    json.dump(ref_id_db, open(f"{work}/ref_id_db.json", "w"))
    json.dump({k: v for k, v in ref_db.items() if k != drop_md5},
              open(f"{work}/ref_db_minus1.json", "w"))
    json.dump({k: v for k, v in ref_id_db.items() if k != drop_md5},
              open(f"{work}/ref_id_db_minus1.json", "w"))
    with open(f"{work}/queries.fna", "w") as fh:                 # headers = bare md5
        for md5 in sorted(want_md5):
            fh.write(f">{md5}\n{seqs[md5]}\n")
    with open(f"{work}/midas_empty.csv", "w") as fh:             # no external prior
        fh.write("seq,Kingdom,Phylum,Class,Order,Family,Genus,Species,rel_ab\n")

    meta = {"genomes": gids, "region": region, "query_md5_of": query_md5_of,
            "dropped_genome": gids[0], "dropped_md5": drop_md5}
    json.dump(meta, open(meta_p, "w"), indent=2)
    print(f"[fixture] {len(gids)} refs, {len(want_md5)} {region} queries, "
          f"holdout={gids[0]} -> {work}", flush=True)
    return meta


# --------------------------------------------------------------------------- runs
def run_aligner(work, out, ref_db, ref_id_db, backend, *, k2=5, gpu_cand=10, k1=0):
    cmd = [PY, f"{PRFBA}/align_hits.py",
           "--fasta", f"{work}/queries.fna",
           "--ref-db", f"{work}/{ref_db}", "--ref-id-db", f"{work}/{ref_id_db}",
           "--taxonomy-csv", f"{work}/midas_empty.csv",
           "--outdir", out, "--k2", str(k2), "--gpu-cand", str(gpu_cand),
           "--prefilter-k", str(k1 or len(json.load(open(f"{work}/{ref_db}")))),
           "--min-ref-len", "200", "--workers", "4", "--force-backend", backend]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout[-3000:], r.stderr[-3000:])
        raise SystemExit(f"[abort] aligner failed ({backend}, rc={r.returncode})")
    return (json.load(open(f"{out}/asv_top20_alignment_hits.json")),
            json.load(open(f"{out}/method_selection.json")))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--workdir", default="/tmp/ref_db_override_test")
    ap.add_argument("--n-genomes", type=int, default=8)
    ap.add_argument("--region", default="V4")
    ap.add_argument("--backends", default="cpu,cuda")
    ap.add_argument("--k2", type=int, default=5)
    ap.add_argument("--force", action="store_true", help="rebuild the cached fixture")
    args = ap.parse_args()
    work = args.workdir
    backends = [b for b in args.backends.split(",") if b]

    print("== fixture ==", flush=True)
    fx = build_fixture(work, args.n_genomes, args.region, args.force)
    src_of = {v: k for k, v in fx["query_md5_of"].items()}
    refs = set(json.load(open(f"{work}/ref_db.json")))
    refs1 = set(json.load(open(f"{work}/ref_db_minus1.json")))

    print("== guard ==", flush=True)
    rc = subprocess.run([PY, f"{PRFBA}/align_hits.py", "--k2", "20", "--gpu-cand", "10",
                         "--explain", "--outdir", f"{work}/guard"],
                        capture_output=True, text=True).returncode
    assert rc == 2, f"--gpu-cand < --k2 must exit 2, got {rc}"
    ok("--gpu-cand < --k2 exits 2")
    rc = subprocess.run([PY, f"{PRFBA}/align_hits.py", "--k2", "200", "--gpu-cand", "400",
                         "--explain", "--outdir", f"{work}/guard"],
                        capture_output=True, text=True).returncode
    assert rc == 0, f"the spec's --k2 200 --gpu-cand 400 must be accepted, got {rc}"
    sel0 = json.load(open(f"{work}/guard/method_selection.json"))
    assert list(sel0["params"]) == PARAM_KEYS, list(sel0["params"])
    assert sel0["params"]["ref_db"] is None and sel0["params"]["ref_id_db"] is None
    ok("params carry ref_db/ref_id_db (null without the flags)")

    hits = {}
    for backend in backends:
        print(f"== restricted DB, backend={backend} ==", flush=True)
        full, sel = run_aligner(work, f"{work}/hits_{backend}", "ref_db.json",
                                "ref_id_db.json", backend, k2=args.k2)
        hits[backend] = full
        assert sel["n_unique_refs"] == len(refs), (sel["n_unique_refs"], len(refs))
        assert sel["params"]["ref_db"].endswith("ref_db.json")
        ok("aligned against the --ref-db only", f"n_unique_refs={sel['n_unique_refs']}")

        for qid, e in full.items():
            assert set(e) == ENTRY_KEYS, set(e) ^ ENTRY_KEYS
            for h in e["top20"]:
                assert set(h) == HIT_KEYS, set(h) ^ HIT_KEYS
                assert set(h["lineage"]) == RANKS
                assert abs(h["identity"] - h["n_matches"] / h["aligned_len"]) < 5e-5
                assert h["md5"] in refs, h["md5"]
                assert h["organism"] and h["genome_id"] and h["taxon_id"] and h["feature_id"]
        ok("hits schema + identity == n_matches/aligned_len + id-db join",
           f"{len(full)} queries x {len(full[next(iter(full))]['top20'])} hits")

        for qmd5, e in full.items():
            h = e["top20"][0]
            assert h["genome_id"] == src_of[qmd5] and h["identity"] == 1.0, (qmd5, h)
        ok("in-DB source genome is rank-1 at identity 1.000")

        assert all(e["midas_taxonomy"] == "" and e["rel_ab"] == 0.0 for e in full.values())
        ok("header-only MiDAS CSV -> no external lineage prior")

        print(f"== restricted DB minus the holdout, backend={backend} ==", flush=True)
        m1, sel1 = run_aligner(work, f"{work}/hits_minus1_{backend}", "ref_db_minus1.json",
                               "ref_id_db_minus1.json", backend, k2=args.k2)
        assert sel1["n_unique_refs"] == len(refs1)
        assert not any(h["md5"] == fx["dropped_md5"] or h["genome_id"] == fx["dropped_genome"]
                       for e in m1.values() for h in e["top20"])
        qd = fx["query_md5_of"][fx["dropped_genome"]]
        assert full[qd]["best_identity"] == 1.0 and m1[qd]["best_identity"] < 1.0
        ok("holdout genome absent from every hit",
           f"its own amplicon: best_identity 1.0 -> {m1[qd]['best_identity']}")

    if len(hits) > 1:
        a, b = list(hits)[:2]
        assert all(hits[a][q]["top20"][i]["md5"] == hits[b][q]["top20"][i]["md5"]
                   and hits[a][q]["top20"][i]["align_score"] == hits[b][q]["top20"][i]["align_score"]
                   for q in hits[a] for i in range(len(hits[a][q]["top20"])))
        ok(f"{a} and {b} backends return identical top-{args.k2}")

    print("\n[PASS] --ref-db / --ref-id-db restrict the reference, the schema and the "
          "identity definition are unchanged, and the k2/gpu_cand guard exits 2.", flush=True)


if __name__ == "__main__":
    main()
