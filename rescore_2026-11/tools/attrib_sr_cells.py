#!/usr/bin/env python
"""W8 for the shared-reference cell tables (old = pre_repair copy, new = regenerated).

    attrib_sr_cells.py --root {data,data_family,data_c2} --old OLD.csv.gz --new NEW.csv.gz --out REPORT.json

Rows are joined on (amplicon, fold, density, arm, space). Every changed row gets the first matching reason:
  source_changed            src_genome_id is one of the 17 changed eligible genomes (truth changed)
  prfba_prediction_changed  prFBA arm whose run prediction for this md5 changed (R3.4 attributed every
                            such change to a changed genome among the reps)
  picrust2_ref_changed      PICRUSt2 arm whose reference(s) carry a changed genome (traits re-screened)
  nn1_changed               base_nn1 whose rank-1 neighbour (prFBA predictions' rank1_genome) is changed
  ref_set_changed           base_constant/base_family_donor/base_clade_max whose fold x density reference set
                            (union of the manifests) contains a changed genome (prevalence / donor / clade)
  random_abstain            base_random_abstain whose forced-arm prediction changed or whose fold's shipped
                            coverage (the mask budget) moved
  accessory_vocab_moved     ONLY n_truth_accessory/f1_accessory changed, and the fold x density x space
                            accessory vocabulary (EC prevalence in [0.05, 0.95] over the reference set, as
                            score_folds.FoldContext) differs old vs new (recomputed from the old/new
                            _cache_panel_ec.json, the manifests and function_maps ec_rxn)
  unattributed              none of the above
Strict invariants (counted separately): picrust2 rows on UNTOUCHED refs and prFBA rows with UNCHANGED
predictions must be identical unless the source changed.
"""
import argparse, collections, csv, gzip, io, json, os, sys
import orjson
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference"
B = f"{W}/pre_repair/region_validation/shared_reference"
ch = set(map(str, json.load(open(f"{W}/expected/changed_sets.json"))["changed_eligible"]))
aff = json.load(open(f"{W}/control/r3/r32_affected_refs.json"))
L = lambda p: orjson.loads(open(p, "rb").read())


def read(path):
    """key -> raw line (fields are plain: the writer never quotes; checked below)."""
    rows = {}
    with gzip.open(path, "rt", newline="") as fh:
        hdr = fh.readline().rstrip("\r\n").split(",")
        ix = {h: i for i, h in enumerate(hdr)}
        ka = [ix["amplicon"], ix["fold"], ix["density"], ix["arm"], ix["space"]]
        n = len(hdr)
        for line in fh:
            line = line.rstrip("\r\n")
            if '"' in line:
                raise SystemExit(f"quoted field in {path}: {line[:200]}")
            f = line.split(",")
            if len(f) != n:
                raise SystemExit(f"field count {len(f)} != {n} in {path}")
            k = (f[ka[0]], f[ka[1]], f[ka[2]], f[ka[3]], f[ka[4]])
            if k in rows:
                raise SystemExit(f"duplicate key {k} in {path}")
            rows[k] = line
    return hdr, ix, rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True); ap.add_argument("--old", required=True)
    ap.add_argument("--new", required=True); ap.add_argument("--out", required=True)
    a = ap.parse_args()
    root = a.root
    affected = set(aff[root]["affected"]); untouched = set(aff[root]["unaffected"])
    # reference sets per (fold, density) and per ref name
    refset = {}
    ref_has_changed = {}
    for name in os.listdir(f"{SR}/{root}/refs"):
        mp = f"{SR}/{root}/refs/{name}/manifest.json"
        if not os.path.exists(mp) or name.startswith("recipe"):
            continue
        m = L(mp)
        g = {str(x).removeprefix("g_") for x in m["reference_genomes"]}
        ref_has_changed[name] = bool(g & ch)
        dens = "%03d" % int(round(float(m["density"]) * 100))
        if root == "data_c2":
            for f in range(5):
                refset.setdefault((str(f), dens), set()).update(g)
        else:
            refset.setdefault((str(int(m["fold"])), dens), set()).update(g)
    # accessory vocabulary per (fold, density, space), old vs new (score_folds.FoldContext semantics)
    croot = "data" if root == "data_c2" else root
    uni = set(json.load(open(f"{SR}/{croot}/panel_ec_universe.json"))["ec"])
    ecn = L(f"{SR}/{croot}/_cache_panel_ec.json"); eco = L(f"{B}/{croot}/_cache_panel_ec.json")
    ec_rxn = {e: set(r) for e, r in L("/home/freiburger/Documents/prFBA/region_validation/picrust2_comparison/data/function_maps.json")["ec_rxn"].items()}
    def acc(cache, ids):
        out = {}
        sets = {g: set(cache[g]["ec"]) & uni for g in ids if g in cache}
        n = max(1, len(sets))
        for space in ("ec", "rxn"):
            cnt = collections.Counter()
            for g, es in sets.items():
                cnt.update(es if space == "ec" else set().union(*[ec_rxn.get(e, set()) for e in es]) if es else set())
            out[space] = {k for k, v in cnt.items() if 0.05 <= v / n <= 0.95}
        return out
    acc_moved = {}
    for (fold, dens), ids in refset.items():
        ao, an = acc(eco, ids), acc(ecn, ids)
        for space in ("ec", "rxn"):
            acc_moved[(fold, dens, space)] = ao[space] != an[space]
    ACC = {"n_truth_accessory", "f1_accessory"}
    # prFBA prediction changes per run/arm and rank1 genomes
    pf_changed, rank1, cov_moved = {}, {}, {}
    for run in sorted(os.listdir(f"{SR}/{root}/prfba")):
        po = L(f"{B}/{root}/prfba/{run}/predictions.json")["pred"]
        pn = L(f"{SR}/{root}/prfba/{run}/predictions.json")["pred"]
        for arm in pn:
            pf_changed[(run, arm)] = {m for m in set(po[arm]) | set(pn[arm]) if po[arm].get(m) != pn[arm].get(m)}
            cov = lambda p: sum(1 for v in p.values() if not v["abstained"])
            cov_moved[(run, arm)] = cov(po[arm]) != cov(pn[arm])
            for m, v in pn[arm].items():
                rank1.setdefault(run, {})[m] = {str(v.get("rank1_genome")), str((po[arm].get(m) or {}).get("rank1_genome"))}
        del po, pn
    hdr, ix, old = read(a.old)
    hdr2, ix2, new = read(a.new)
    rep = {"root": root, "old": a.old, "new": a.new, "header_equal": hdr == hdr2,
           "n_old": len(old), "n_new": len(new),
           "one_sided": len(set(old) ^ set(new)), "n_changed": 0, "reasons": collections.Counter(),
           "by_arm": collections.defaultdict(collections.Counter), "strict_violations": collections.Counter(),
           "unattributed_examples": [], "changed_columns": collections.Counter()}
    for k in old.keys() & new.keys():
        if old[k] == new[k]:
            continue
        fo, fn = old[k].split(","), new[k].split(",")
        rep["n_changed"] += 1
        for i, h in enumerate(hdr):
            if fo[i] != fn[i]:
                rep["changed_columns"][h] += 1
        amp, fold, dens_s, arm, space = k
        src = fn[ix["src_genome_id"]]; md5 = fn[ix["md5"]]
        dens = "%03d" % int(round(float(dens_s) * 100))
        if root == "data_c2":
            run = f"incself_f{int(fold):02d}_d100"
        else:
            run = f"f{int(fold):02d}_d{dens}" + ("_forced" if arm.endswith("_forced") else "")
        refs_of_cell = set()
        dom = "arc" if fn[ix["domain"]] == "Archaea" else "bac"
        if arm == "picrust2":
            refs_of_cell = {f"incself_{dom}_d100"} if root == "data_c2" else {f"f{int(fold):02d}_{dom}_d{dens}"}
        elif arm == "picrust2_nsti_routed":
            refs_of_cell = {f"f{int(fold):02d}_bac_d100", f"f{int(fold):02d}_arc_d100"}
        reason = None
        cols = {hdr[i] for i in range(len(hdr)) if fo[i] != fn[i]}
        acc_only = cols <= ACC and acc_moved.get((fold, dens, space), False)
        if src in ch:
            reason = "source_changed"
        elif acc_only:
            reason = "accessory_vocab_moved"
        elif arm.startswith("prfba_"):
            key = arm[len("prfba_"):]
            if md5 in pf_changed.get((run, key), ()):
                reason = "prfba_prediction_changed"
            else:
                rep["strict_violations"]["prfba_unchanged_prediction_but_row_changed"] += 1
        elif arm.startswith("picrust2"):
            if any(ref_has_changed.get(r) for r in refs_of_cell):
                reason = "picrust2_ref_changed"
            elif refs_of_cell & untouched:
                rep["strict_violations"]["picrust2_untouched_ref_row_changed"] += 1
        elif arm == "base_nn1":
            pr = run if not run.endswith("_forced") else run[:-7]
            if rank1.get(pr, {}).get(md5, set()) & ch:
                reason = "nn1_changed"
        elif arm in ("base_constant", "base_family_donor", "base_clade_max"):
            if refset.get((fold, dens), set()) & ch:
                reason = "ref_set_changed"
        elif arm == "base_random_abstain":
            fr = f"f{int(fold):02d}_d100_forced" if root != "data_c2" else run
            shipped = f"f{int(fold):02d}_d100" if root != "data_c2" else run
            if (md5 in pf_changed.get((fr, "exact_tie_forced"), ()) or cov_moved.get((shipped, "exact_tie"))
                    or cov_moved.get((fr, "exact_tie_forced"))):
                reason = "random_abstain"
        reason = reason or "unattributed"
        rep["reasons"][reason] += 1
        rep["by_arm"][arm][reason] += 1
        if reason == "unattributed" and len(rep["unattributed_examples"]) < 15:
            rep["unattributed_examples"].append({"key": k, "old": dict(zip(hdr, fo)), "new": dict(zip(hdr, fn))})
    # strict identity: picrust2 rows of untouched refs (any, not only changed) -> count changed among them
    rep["untouched_refs"] = sorted(untouched)
    rep["accessory_vocab_moved"] = {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in sorted(acc_moved.items())}
    rep["n_unattributed"] = rep["reasons"].get("unattributed", 0)
    rep["pass"] = (rep["header_equal"] and rep["one_sided"] == 0 and rep["n_unattributed"] == 0
                   and not rep["strict_violations"])
    json.dump(rep, open(a.out, "w"), indent=1, default=lambda o: dict(o) if isinstance(o, collections.Counter) else str(o))
    print(json.dumps({k: rep[k] for k in ("root", "n_old", "n_new", "one_sided", "n_changed", "n_unattributed", "pass")}
                     | {"reasons": dict(rep["reasons"]), "strict": dict(rep["strict_violations"])}))
    return 0 if rep["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
