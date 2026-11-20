#!/usr/bin/env python
"""R3.2 checks: refs refreshed without re-running FastTree / hmmalign / probe placement.

- sha256 before/after over every ref file: the changed set per affected ref must be within
  {traits_EC.txt.gz, manifest.json, traits_16S.txt.gz, hsp_screen/*, <name>.hsp_screen.log}; unaffected refs untouched;
- traits_16S.txt.gz decompressed-equal to the old copy;
- manifest: tip_labels_sha256, model_string, reference_genomes, n_tips, msa_width, set_equality_ok equal;
- traits_EC decompressed: rows differ only for the 17 changed eligible genomes (label -> genome via panel);
  column deltas (prevalence screen) reported;
- build logs: every build read 'EC sets from <panel dir>/_cache_panel_ec.json'; none derived or fell through.
"""
import gzip, json, os, re, sys, collections
W = "/home/freiburger/Documents/prFBA_rescore_2026-11"
SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference"
B = f"{W}/pre_repair/region_validation/shared_reference"
C = f"{W}/control/r3"
ch = set(map(str, json.load(open(f"{W}/expected/changed_sets.json"))["changed_eligible"]))
aff = json.load(open(f"{C}/r32_affected_refs.json"))


def sha(path):
    d = {}
    for line in open(path):
        h, p = line.rstrip("\n").split("  ", 1)
        d[p] = h
    return d


import argparse
_ap = argparse.ArgumentParser(); _ap.add_argument("--after", default=f"{C}/r32_refs_after.sha256")
_ap.add_argument("--roots", default="data,data_family,data_c2"); _ap.add_argument("--out", default=f"{C}/r32_checks.json")
_a = _ap.parse_args(); ROOTS = _a.roots.split(",")
before, after = sha(f"{C}/r32_refs_before.sha256"), sha(_a.after)
rep = {"refs": {}, "unaffected_changed_files": [], "problems": []}
allowed = lambda name, rel: rel in ("traits_EC.txt.gz", "manifest.json", "traits_16S.txt.gz", f"{name}.hsp_screen.log") or rel.startswith("hsp_screen/")
changed = sorted(p for p in set(before) | set(after) if before.get(p) != after.get(p))
by_ref = collections.defaultdict(list)
for p in changed:
    parts = p.split("/")            # root/refs/name/rel...
    by_ref[(parts[0], parts[2])].append("/".join(parts[3:]))
for root in ROOTS:
    for name in aff[root]["unaffected"]:
        if (root, name) in by_ref:
            rep["unaffected_changed_files"].append([root, name, by_ref[(root, name)]])
for (root, name), rels in by_ref.items():
    if name not in aff[root]["affected"] and not name.startswith("recipe") and "/" not in name:
        pass
for root in ROOTS:
    panel = json.load(open(f"{SR}/{root if root != 'data_c2' else 'data'}/panel_eligible.json"))["genomes"]
    lab2g = {v["label"]: g for g, v in panel.items()}
    for name in aff[root]["affected"]:
        r = {}
        rels = by_ref.get((root, name), [])
        bad = [x for x in rels if not allowed(name, x)]
        r["changed_files"] = sorted(rels)
        r["disallowed_changes"] = bad
        d, b = f"{SR}/{root}/refs/{name}", f"{B}/{root}/refs/{name}"
        r["traits_16S_decompressed_equal"] = gzip.open(f"{d}/traits_16S.txt.gz").read() == gzip.open(f"{b}/traits_16S.txt.gz").read()
        mo, mn = json.load(open(f"{b}/manifest.json")), json.load(open(f"{d}/manifest.json"))
        keys = ("tip_labels_sha256", "model_string", "reference_genomes", "n_tips", "msa_width", "set_equality_ok",
                "hmm_leng", "msa_mean_gap_frac", "mean_pairwise_tip_distance", "genus_monophyly_frac", "family_monophyly_frac")
        r["manifest_equal"] = {k: mo.get(k) == mn.get(k) for k in keys}
        r["set_equality_ok_new"] = mn.get("set_equality_ok")
        r["trait_cols"] = [mo.get("trait_cols"), mn.get("trait_cols")]
        r["n_dropped"] = [mo.get("n_dropped"), mn.get("n_dropped")]
        r["retry_used"] = [mo.get("retry_used"), mn.get("retry_used")]
        r["fasttree_seconds_new"] = mn.get("fasttree_seconds")
        def rows(p):
            with gzip.open(p, "rt") as fh:
                hdr = fh.readline().rstrip("\n").split("\t")[1:]
                out = {}
                for line in fh:
                    f = line.rstrip("\n").split("\t")
                    out[f[0]] = {e for e, v in zip(hdr, f[1:]) if v == "1"}
            return hdr, out
        ho, ro = rows(f"{b}/traits_EC.txt.gz"); hn, rn = rows(f"{d}/traits_EC.txt.gz")
        cols_same = set(ho) == set(hn)
        r["traits_cols_added"] = sorted(set(hn) - set(ho)); r["traits_cols_removed"] = sorted(set(ho) - set(hn))
        common = set(ho) & set(hn)
        diff_rows = sorted(l for l in set(ro) | set(rn) if (ro.get(l, set()) & common) != (rn.get(l, set()) & common))
        r["rows_same_set"] = set(ro) == set(rn)
        r["n_rows_differ_on_common_cols"] = len(diff_rows)
        r["rows_differ_genomes_not_in_17"] = sorted(lab2g.get(l, l) for l in diff_rows if lab2g.get(l, l).removeprefix("g_") not in ch)
        r["n_changed_genomes_in_ref"] = len({str(g).removeprefix("g_") for g in mn["reference_genomes"]} & ch)
        log = open(f"{W}/logs/r3/build_{root}_{name}.log").read()
        panel_dir = f"{SR}/{root if root != 'data_c2' else 'data'}"
        r["ec_source_ok"] = f"EC sets from {panel_dir}/_cache_panel_ec.json" in log and "deriving EC sets" not in log and "falling through" not in log
        r["fasttree_reused"] = "FastTree: reusing" in log
        r["exit"] = re.findall(r"^EXIT=(\d+)", log, re.M)
        ok = (not bad and r["traits_16S_decompressed_equal"] and all(r["manifest_equal"][k] for k in ("tip_labels_sha256", "model_string", "reference_genomes", "n_tips", "msa_width", "set_equality_ok"))
              and r["set_equality_ok_new"] and r["rows_same_set"] and not r["rows_differ_genomes_not_in_17"] and r["ec_source_ok"] and r["fasttree_reused"] and r["exit"] == ["0"])
        r["ok"] = ok
        if not ok:
            rep["problems"].append(f"{root}/{name}")
        rep["refs"][f"{root}/{name}"] = r
rep["n_affected"] = len(rep["refs"]); rep["n_ok"] = sum(r["ok"] for r in rep["refs"].values())
rep["panel_ec_sets_written"] = [p for p in (f"{SR}/data/panel_ec_sets.json", f"{SR}/data_family/panel_ec_sets.json") if os.path.exists(p)]
rep["roots"] = ROOTS
json.dump(rep, open(_a.out, "w"), indent=1)
print(json.dumps({k: rep[k] for k in ("n_affected", "n_ok", "problems", "unaffected_changed_files", "panel_ec_sets_written")}))
summ = collections.Counter()
for k, r in rep["refs"].items():
    summ["cols_changed"] += bool(r["traits_cols_added"] or r["traits_cols_removed"])
    summ["rows_differ"] += r["n_rows_differ_on_common_cols"]
    summ["manifest_descriptives_differ"] += not all(r["manifest_equal"].values())
print(dict(summ))
sys.exit(0 if not rep["problems"] and not rep["unaffected_changed_files"] and not rep["panel_ec_sets_written"] else 1)
