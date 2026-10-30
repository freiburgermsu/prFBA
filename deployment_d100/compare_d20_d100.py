#!/usr/bin/env python
"""
compare_d20_d100.py -- post-run validation and 20-hit vs 100-hit comparison of the two
real-amplicon deployments (feeds RESULTS.md). Read-only on the deployment directories;
writes only under --out-dir (JSON) and --scratch (control selections).

Per deployment (A codif, B EmilyKin):
  1. validate: hits per record, ranks, ASV coverage (FASTA / hits / selection / manifest /
     per-ASV genome JSONs / provenance / genome_objects / taxonomies), tie-cloud truncation.
  2. hits concordance: the 20-hit record (exhaustive GPU SW) vs the first 20 hits of the
     depth-100 record (edlib top-500 prefilter + CUDA SW).
  3. control selections, same code path as build_synthetic_genomes_parallel.exact_select
     (S.bvbrc_gene_provider + S.gene_count_from_provider [P3 on] + S.build_audit), against
     the depth-100 run's small tie-set gene cache (+ entries streamed read-only from the
     shared cache for any control tie-set genome it lacks), network disabled, with a
     tie-set guard:
       ctl_d20hits_p3   : the 20-hit hits, P3 on          (P3 / cache effect vs published)
       ctl_d100top20_p3 : the depth-100 hits cut to 20, P3 on (aligner-path effect)
     The depth-100 selection itself is the build's (asv_reference_selection_d100.json).
  4. per-ASV selected-set changes between selection pairs (genome ids and species).

    ~/Documents/py_venv/bin/python deployment_d100/compare_d20_d100.py \
        --out-dir deployment_d100/compare --scratch <scratchpad>
"""
from __future__ import annotations

import argparse
import collections
import json
import os
import subprocess
import sys
import time
import urllib.request

PRFBA = "/home/freiburger/Documents/prFBA"
sys.path.insert(0, PRFBA)
sys.path.insert(0, os.path.join(PRFBA, "deployment_d100"))
import select_references as S  # noqa: E402
import tieset_gene_cache as T  # noqa: E402

PY = "/home/freiburger/Documents/py_venv/bin/python"
SHARED_CACHE = os.path.join(PRFBA, "bvbrc_cache/genome_gene_families.json")
C = "/home/freiburger/Documents/codiffusion_bioreactor"
E = "/home/freiburger/Documents/EmilyKin"
DEP = {
    "codif": dict(
        fasta=f"{C}/new_data/dna-sequences_codif_all.fasta",
        hits20=f"{C}/bvbrc_alignment_hits/asv_top20_alignment_hits_gpu.json",
        sel20={"d20_54d1150_p3on(KBase)": f"{C}/bvbrc_alignment_hits/asv_reference_selection_gpu.54d1150.json",
               "d20_p3off": f"{C}/bvbrc_alignment_hits/asv_reference_selection_gpu.json"},
        out100=f"{C}/bvbrc_alignment_hits_d100", syn100=f"{C}/synthetic_genomes_d100"),
    "emilykin": dict(
        fasta=f"{E}/asvs.fasta",
        hits20=f"{E}/bvbrc_alignment_hits/asv_top20_alignment_hits_gpu.json",
        sel20={"d20_published_p3off": f"{E}/bvbrc_alignment_hits/asv_reference_selection_gpu.json",
               "d20_p3on_june": f"{E}/bvbrc_alignment_hits/asv_reference_selection_gpu.p3.json"},
        out100=f"{E}/bvbrc_alignment_hits_d100", syn100=f"{E}/synthetic_genomes_d100"),
}


def recs(d):
    return {k: v for k, v in d.items() if not k.startswith("_")}


def fasta_ids(p):
    return [l[1:].split()[0] for l in open(p) if l.startswith(">")]


def tie_cloud_cut(hits, depth):
    """records whose last (depth-th) hit ties rank 1 on identity AND align_score."""
    n = 0
    for r in hits.values():
        h = r["top20"]
        if len(h) >= depth and h[depth - 1]["identity"] == h[0]["identity"] \
                and h[depth - 1]["align_score"] == h[0]["align_score"]:
            n += 1
    return n


def validate(name, cfg, hits100):
    out = {}
    fids = fasta_ids(cfg["fasta"])
    L = collections.Counter(len(r["top20"]) for r in hits100.values())
    out["fasta_asvs"] = len(fids)
    out["fasta_unique"] = len(set(fids))
    out["hits_records"] = len(hits100)
    out["hits_per_record"] = {str(k): v for k, v in sorted(L.items())}
    out["ranks_contiguous"] = all([h["rank"] for h in r["top20"]] == list(range(1, len(r["top20"]) + 1))
                                  for r in hits100.values())
    out["hits_asvs_eq_fasta"] = set(hits100) == set(fids)
    sel = recs(json.load(open(f"{cfg['out100']}/asv_reference_selection_d100.json")))
    answered = {a for a, r in sel.items() if r["n_selected"] > 0}
    out["selection_records"] = len(sel)
    out["selection_asvs_eq_hits"] = set(sel) == set(hits100)
    out["answered"] = len(answered)
    man = json.load(open(f"{cfg['syn100']}/synthetic_genome_manifest.json"))
    mrec = recs(man)
    out["manifest_meta"] = man.get("_meta")
    out["manifest_records"] = len(mrec)
    out["manifest_status"] = dict(collections.Counter(
        (v.get("status") if isinstance(v, dict) else str(type(v))) for v in mrec.values()))
    built = {a for a, v in mrec.items() if isinstance(v, dict) and v.get("status") == "built"}
    out["manifest_built_eq_answered"] = built == answered
    out["manifest_member_genomes_failed"] = sum(len(v.get("failed") or []) for v in mrec.values())
    out["manifest_records_with_failed_members"] = sum(1 for v in mrec.values() if v.get("failed"))
    out["manifest_genomes_ne_selection"] = sum(
        1 for a in built if sorted(map(str, mrec[a].get("genomes") or [])) !=
        sorted(str(s["genome_id"]) for s in sel[a]["selected"]))
    out["manifest_n_genomes_ne_selection"] = sum(1 for a in built if mrec[a].get("n_genomes") != sel[a]["n_selected"])
    syn = cfg["syn100"]
    top = os.listdir(syn)
    top_json = {f[:-5] for f in top if f.endswith(".json") and os.path.isfile(os.path.join(syn, f))}
    out["top_level_entries"] = len(top)
    out["top_level_json_files"] = len(top_json)
    out["top_level_non_asv_json"] = sorted(top_json - answered)
    out["top_level_dirs"] = sorted(f for f in top if os.path.isdir(os.path.join(syn, f)))
    out["top_level_asv_json_eq_answered"] = (top_json & set(sel)) == answered
    for sub in ("provenance", "genome_objects", "taxonomies", "ASVset_taxonomies"):
        p = os.path.join(syn, sub)
        ids = {f[:-5] for f in os.listdir(p) if f.endswith(".json")}
        out[f"{sub}_files"] = len(os.listdir(p))
        out[f"{sub}_eq_answered"] = ids == answered
    out["empty_or_zero_byte_genome_json"] = sum(1 for a in answered
                                                if os.path.getsize(os.path.join(syn, a + ".json")) == 0)
    out["tie_cloud_cut_at_100"] = tie_cloud_cut(hits100, 100)
    side = json.load(open(f"{cfg['out100']}/gene_families_tieset.json.tieset.json"))
    out["tieset_n"] = side["n_tieset"]
    out["tieset_present_in_shared"] = side["n_present_in_source"]
    return out, sel


def concordance(hits20, hits100):
    same_list = same_set = same_top1 = same_best = 0
    n = 0
    for a, r in hits20.items():
        n += 1
        A = [(h["genome_id"], h["identity"], h["align_score"]) for h in r["top20"]]
        B = [(h["genome_id"], h["identity"], h["align_score"]) for h in hits100[a]["top20"][:len(A)]]
        same_list += A == B
        same_set += set(A) == set(B)
        same_best += (r["top20"][0]["identity"], r["top20"][0]["align_score"]) == \
            (hits100[a]["top20"][0]["identity"], hits100[a]["top20"][0]["align_score"])
        same_top1 += A[:1] == B[:1]
    return dict(n=n, identical_ordered_top20=same_list, identical_top20_set=same_set,
                identical_rank1_hit=same_top1, identical_best_identity_and_score=same_best,
                hits20_lens=dict(collections.Counter(len(r["top20"]) for r in hits20.values())))


def control_select(hits, cache_path, tie_ok, label, scratch):
    """exact_select's body with P3 on; network disabled; tie-set guard."""
    net = {"calls": 0}

    def no_net(*a, **k):
        net["calls"] += 1
        raise RuntimeError("network disabled in control selection")
    orig_open = urllib.request.urlopen
    urllib.request.urlopen = no_net
    try:
        provider = S.bvbrc_gene_provider(cache_path)
        f = S.gene_count_from_provider(provider)
        g = {"q": 0, "outside": 0}

        def gcf(gid):
            g["q"] += 1
            g["outside"] += str(gid) not in tie_ok
            return f(gid)
        audit, summary = S.build_audit(hits, gene_provider=provider, gene_count_fn=gcf)
    finally:
        urllib.request.urlopen = orig_open
    path = os.path.join(scratch, f"{label}.json")
    with open(path, "w") as fh:
        json.dump({"_meta": {"gene_data": "bvbrc", "generated_from": label, "summary": summary,
                             "control": True}, **audit}, fh)
    del provider
    return path, dict(gene_count_queries=g["q"], queries_outside_tieset=g["outside"],
                      network_calls=net["calls"], abstained=summary["abstained"],
                      mean_handful=summary["mean_handful"])


def run_stats(sel, hits, out):
    r = subprocess.run([PY, os.path.join(PRFBA, "deployment_d100/deployment_stats.py"), "--selection", sel,
                        "--hits", hits, "--out", out, "--force"], capture_output=True, text=True)
    print(r.stdout.strip(), r.stderr.strip()[-2000:], flush=True)
    return json.load(open(out))


def set_changes(selA, selB, hitsA=None, depthA=20):
    """per-ASV selected-set changes from A to B."""
    c = collections.Counter()
    grew_sizes = []
    cut_and_grew = 0
    cutA = set()
    if hitsA is not None:
        for a, r in hitsA.items():
            h = r["top20"]
            if len(h) >= depthA and h[depthA - 1]["identity"] == h[0]["identity"] and \
                    h[depthA - 1]["align_score"] == h[0]["align_score"]:
                cutA.add(a)
    trans = collections.Counter()
    for a in selA:
        ra, rb = selA[a], selB[a]
        ga = {str(s["genome_id"]) for s in ra["selected"]}
        gb = {str(s["genome_id"]) for s in rb["selected"]}
        sa = {s.get("species") for s in ra["selected"]}
        sb = {s.get("species") for s in rb["selected"]}
        na, nb = len(ga), len(gb)
        bucket = lambda n: "0" if n == 0 else "1" if n == 1 else "2" if n == 2 else ">=3"
        trans[(bucket(na), bucket(nb))] += 1
        if na == 0 and nb == 0:
            c["both_abstain"] += 1
            continue
        if na == 0:
            c["answered_only_in_B"] += 1
            continue
        if nb == 0:
            c["answered_only_in_A"] += 1
            continue
        c["both_answered"] += 1
        if ga == gb:
            c["identical_genome_set"] += 1
        else:
            c["changed_genome_set"] += 1
            if nb > na:
                c["changed_grew"] += 1
                grew_sizes.append(nb - na)
                if ga <= gb:
                    c["grew_superset"] += 1
                if a in cutA:
                    cut_and_grew += 1
            elif nb < na:
                c["changed_shrank"] += 1
            else:
                c["changed_same_size"] += 1
            if sa == sb:
                c["changed_genomes_same_species_set"] += 1
        if sa != sb:
            c["changed_species_set"] += 1
        if (na >= 2) != (nb >= 2):
            c["single_multi_flip"] += 1
    out = dict(c)
    out["grew_by"] = dict(sorted(collections.Counter(
        "1" if x == 1 else "2-4" if x <= 4 else "5-9" if x <= 9 else ">=10" for x in grew_sizes).items()))
    out["grew_total_added_genomes"] = sum(grew_sizes)
    out["max_growth"] = max(grew_sizes or [0])
    if hitsA is not None:
        out["A_records_tie_cloud_cut_at_depth"] = len(cutA)
        out["A_cut_records_that_grew"] = cut_and_grew
    out["size_bucket_transitions"] = {f"{x}->{y}": v for (x, y), v in sorted(trans.items())}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--scratch", required=True)
    ap.add_argument("--only", default=None)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.scratch, exist_ok=True)
    for name, cfg in DEP.items():
        if args.only and name != args.only:
            continue
        t0 = time.perf_counter()
        res = {"deployment": name}
        hits100 = json.load(open(f"{cfg['out100']}/asv_top20_alignment_hits.json"))
        res["validation"], sel100 = validate(name, cfg, hits100)
        hits20 = json.load(open(cfg["hits20"]))
        res["hits20_tie_cloud_cut_at_20"] = tie_cloud_cut(hits20, 20)
        res["concordance_d20_vs_d100_first20"] = concordance(hits20, hits100)
        # control tie sets
        side = json.load(open(f"{cfg['out100']}/gene_families_tieset.json.tieset.json"))
        tie100 = set(side["tieset"])
        trunc = {a: dict(r, top20=r["top20"][:20]) for a, r in hits100.items()}
        need20, _ = T.tie_set(hits20)
        needtr, _ = T.tie_set(trunc)
        missing = (need20 | needtr) - tie100
        res["control_tieset"] = dict(d20hits=len(need20), d100top20=len(needtr),
                                     d20hits_not_in_d100_tieset=len(need20 - tie100),
                                     d100top20_not_in_d100_tieset=len(needtr - tie100))
        cache_path = f"{cfg['out100']}/gene_families_tieset.json"
        tie_ok = set(tie100)
        if missing:
            extra, n_seen = T.stream_extract(SHARED_CACHE, missing)
            small = json.load(open(cache_path))
            n_small = len(small)
            small.update({k: v for k, v in extra.items() if k not in small})
            cache_path = os.path.join(args.scratch, f"{name}_control_gene_cache.json")
            with open(cache_path, "w") as fh:
                json.dump(small, fh)
            del small
            tie_ok |= missing
            res["control_tieset"].update(streamed_extra_present=len(extra), small_cache_entries=n_small,
                                         streamed_entries_seen=n_seen)
        ctl = {}
        for label, H in (("ctl_d20hits_p3", hits20), ("ctl_d100top20_p3", trunc)):
            p, info = control_select(H, cache_path, tie_ok, f"{name}_{label}", args.scratch)
            ctl[label] = dict(selection=p, **info)
        res["controls"] = ctl
        del trunc
        # stats for controls
        stats = {}
        for label in ctl:
            hp = cfg["hits20"] if label == "ctl_d20hits_p3" else f"{cfg['out100']}/asv_top20_alignment_hits.json"
            st = run_stats(ctl[label]["selection"], hp, os.path.join(args.out_dir, f"stats_{name}_{label}.json"))
            stats[label] = {k: st[k] for k in ("answered", "abstained", "abstained_by_flag", "n_genomes_1",
                                               "n_genomes_2", "n_genomes_3plus", "mean_n_genomes_answered")}
            stats[label]["cross_phylum"] = st["cross_phylum"]
            stats[label]["max_n"] = max(int(k) for k in st["n_genomes_distribution"])
        res["control_stats"] = stats
        # set changes
        sels = {k: recs(json.load(open(v))) for k, v in cfg["sel20"].items()}
        sels["ctl_d20hits_p3"] = recs(json.load(open(ctl["ctl_d20hits_p3"]["selection"])))
        sels["ctl_d100top20_p3"] = recs(json.load(open(ctl["ctl_d100top20_p3"]["selection"])))
        sels["d100"] = sel100
        ch = {}
        for a in list(cfg["sel20"]) + ["ctl_d20hits_p3", "ctl_d100top20_p3"]:
            ch[f"{a} -> d100"] = set_changes(sels[a], sels["d100"],
                                             hitsA=hits20 if a != "ctl_d100top20_p3" else hits100, depthA=20)
        for a in cfg["sel20"]:
            ch[f"{a} -> ctl_d20hits_p3"] = set_changes(sels[a], sels["ctl_d20hits_p3"])
        ch["ctl_d20hits_p3 -> ctl_d100top20_p3"] = set_changes(sels["ctl_d20hits_p3"], sels["ctl_d100top20_p3"])
        res["set_changes"] = ch
        res["wall_seconds"] = round(time.perf_counter() - t0, 1)
        with open(os.path.join(args.out_dir, f"compare_{name}.json"), "w") as fh:
            json.dump(res, fh, indent=1, default=str)
        print(f"[compare] {name} done in {res['wall_seconds']} s", flush=True)
        del hits100, hits20, sels, sel100


if __name__ == "__main__":
    main()
