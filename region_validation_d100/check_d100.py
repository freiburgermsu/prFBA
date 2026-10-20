#!/usr/bin/env python
"""Section-8 checks for the depth-100 rerun (see RUNBOOK_d100.md). Read-only: it never
writes anything except the optional --json report.

    check_d100.py align  <tree>     checks 1-4 on a union tree's alignment
    check_d100.py copy   <tree>     a legacy tree's hits are byte copies of its union tree's
    check_d100.py post   <tree>     check 5 (+ the selection rule each output recorded)
    check_d100.py baseline          check 6: sha256 of the cited depth-20 artefacts
    check_d100.py all               every check whose inputs exist

<tree> is a directory name under ~/Documents/prFBA, e.g. region_validation_kmeans_d100.

Checks (RUNBOOK section "Checks"):
  1  every record has exactly 100 hits; method_selection params k2=100, gpu_cand=500,
     prefilter_k=500
  2  the first 20 align_score vectors equal the depth-20 record for 100% of md5s, and a hit
     may differ only where its score's tie group crosses GPU rank 100 (the depth-20 run
     rescored only the GPU top max(gpu_cand=100, k2=20) = 100, so a tie spanning that cut was
     sampled arbitrarily); verified against the persisted shortlist + GPU scores
  3  (original panel) the 316 depth-probe md5s equal hits_k500[:100] exactly
  4  records whose 100th hit is still >= 0.987 (plan_folds._truncation_violations)
  5  taxacc _meta.n_empty_after_exclude per tau; gene-capture pgfam_missing == 0 in both
     self modes; select_mode recorded by selection / taxacc / genecap matches the tree rule
  6  sha256sum -c of the depth-20 baseline
"""
import argparse, csv, hashlib, json, os, subprocess, sys
from collections import Counter

import numpy as np
import orjson

R = "/home/freiburger/Documents/prFBA"
BASELINE = os.path.join(R, "region_validation_d100", "d20_baseline.sha256")
PROBE_HITS = os.path.join(R, "region_validation", "picrust2_comparison", "data", "depth_probe",
                          "hits_k500", "asv_top20_alignment_hits.json")
UNION = ("region_validation_d100", "region_validation_kmeans_d100", "region_validation_random_d100")
TAU = 0.987
D20_GPU_KEEP = 100          # depth-20 run: keep = max(gpu_cand=100, k2=20)


def load(p):
    with open(p, "rb") as fh:
        return orjson.loads(fh.read())


def d20_tree(tree):
    return tree.replace("_legacy", "").replace("_d100", "")


def union_tree(tree):
    return tree.replace("_legacy", "")


def is_legacy(tree):
    return tree.endswith("_legacy")


def report(name, ok, detail):
    print("[%s] %-52s %s" % ("PASS" if ok else "FAIL", name, detail), flush=True)
    return ok


def check_align(tree):
    raw = os.path.join(R, tree, "data", "hits", "raw")
    ms = json.load(open(os.path.join(raw, "method_selection.json")))
    prm = ms["params"]
    ok = report("1 params k2/gpu_cand/prefilter_k", (prm["k2"], prm["gpu_cand"], prm["prefilter_k"]) == (100, 500, 500),
                "k2=%s gpu_cand=%s prefilter_k=%s wall=%ss n_asvs=%s" % (
                    prm["k2"], prm["gpu_cand"], prm["prefilter_k"], ms.get("wall_seconds_total"), ms.get("n_asvs")))
    hits = load(os.path.join(raw, "asv_top20_alignment_hits.json"))
    lens = Counter(len(v["top20"]) for v in hits.values())
    ok &= report("1 every record has 100 hits", set(lens) == {100}, "record lengths %s" % dict(lens))

    # --- 2: prefix vs the depth-20 record, tie-crossing explained by the GPU rank-100 cut
    fasta = os.path.join(R, tree, "data", "combined_amplicons.unique.fasta")
    ids = [l[1:].split()[0] for l in open(fasta) if l.startswith(">")]
    z = np.load(os.path.join(raw, "edlib_top500_shortlist.npz"))
    n = z["n"]; off = np.concatenate([[0], np.cumsum(n)])
    gs = np.load(os.path.join(raw, "gpu_scores.npy"))
    assert len(ids) == len(n) == len(hits), (len(ids), len(n), len(hits))
    d20 = load(os.path.join(R, d20_tree(tree), "data", "hits", "raw", "asv_top20_alignment_hits.json"))
    n_score_bad = n_exact = n_explained = 0
    unexplained = []
    for row, md5 in enumerate(ids):
        a = d20[md5]["top20"]; b = hits[md5]["top20"][:len(a)]
        sa = [h["align_score"] for h in a]
        if sa != [h["align_score"] for h in b]:
            n_score_bad += 1; unexplained.append((md5, "score vector")); continue
        if [h["md5"] for h in a] == [h["md5"] for h in b]:
            n_exact += 1; continue
        g = gs[off[row]:off[row + 1]]
        diff = [j for j in range(len(a)) if a[j]["md5"] != b[j]["md5"]]
        if all(int((g >= sa[j]).sum()) > D20_GPU_KEEP for j in diff):
            n_explained += 1
        else:
            unexplained.append((md5, "hit differs outside a rank-100-crossing tie"))
    N = len(ids)
    ok &= report("2 top-20 score vectors equal depth 20", n_score_bad == 0,
                 "%d/%d equal" % (N - n_score_bad, N))
    ok &= report("2 hit differences only in rank-100-crossing ties", not unexplained,
                 "identical %d, tie-explained %d, unexplained %d %s" % (
                     n_exact, n_explained, len(unexplained), unexplained[:3]))
    del d20

    # --- 3: the depth probe (original panel only: its 316 md5s are in this FASTA)
    if d20_tree(tree) == "region_validation":
        pr = load(PROBE_HITS)
        missing = [k for k in pr if k not in hits]
        eq = sum(pr[k]["top20"][:100] == hits[k]["top20"] for k in pr if k in hits)
        ok &= report("3 probe md5s == hits_k500[:100]", eq == len(pr) and not missing,
                     "%d/%d identical, %d absent from this record" % (eq, len(pr), len(missing)))

    # --- 4: truncation violations (same definition as plan_folds._truncation_violations)
    sys.path.insert(0, os.path.join(R, "region_validation", "shared_reference", "scripts"))
    viol = sum(1 for v in hits.values() if len(v["top20"]) >= 100 and (v["top20"][-1].get("identity") or 0.0) >= TAU)
    try:
        import plan_folds
        nv, _ = plan_folds._truncation_violations(os.path.join(raw, "asv_top20_alignment_hits.json"),
                                                  100, ms["n_unique_refs"], TAU)
        same = nv == viol
    except Exception as exc:          # informational; the inline count is authoritative
        nv, same = "n/a (%s)" % exc, True
    report("4 records with 100th hit >= %.3f" % TAU, same,
           "%d of %d md5 records (%.1f%%); plan_folds: %s" % (viol, N, 100 * viol / N, nv))
    return ok


def check_copy(tree):
    if not is_legacy(tree):
        return report("copy %s" % tree, False, "not a legacy tree")
    ok = True
    for f in ("asv_top20_alignment_hits.json", "asv_alignment_summary.csv", "method_selection.json", "align.done"):
        a = os.path.join(R, union_tree(tree), "data", "hits", "raw", f)
        b = os.path.join(R, tree, "data", "hits", "raw", f)
        same = os.path.exists(b) and not os.path.islink(b) and \
            subprocess.run(["cmp", "-s", a, b]).returncode == 0
        nlink = os.stat(b).st_nlink if os.path.exists(b) else None
        ok &= report("copy %s" % f, same and nlink == 1, "identical=%s links=%s" % (same, nlink))
    return ok


def check_post(tree):
    d = os.path.join(R, tree, "data")
    want = "legacy" if is_legacy(tree) else "default_union"
    want_sel = "legacy" if is_legacy(tree) else "exact-tie"
    ok = True
    sel_mode = None
    with open(os.path.join(d, "selection.json"), "rb") as fh:
        head = fh.read(1 << 20).decode("utf-8", "ignore")
    import re
    m = re.search(r'"select_mode":\s*"([^"]+)"', head)
    sel_mode = m.group(1) if m else None
    ok &= report("5 selection.json select_mode", sel_mode == want_sel, "%s (want %s)" % (sel_mode, want_sel))
    tx = json.load(open(os.path.join(d, "taxacc_region_summary.json")))["_meta"]
    ok &= report("5 taxacc select_mode", tx.get("select_mode") == want, "%s (want %s)" % (tx.get("select_mode"), want))
    report("5 n_empty_after_exclude", True, "%s of n_total %s" % (tx.get("n_empty_after_exclude"), tx.get("n_total")))
    gc = json.load(open(os.path.join(d, "gene_capture_by_region.json")))["_meta"]
    ok &= report("5 genecap select_mode", gc.get("select_mode") == want, "%s (want %s)" % (gc.get("select_mode"), want))
    miss = Counter()
    with open(os.path.join(d, "gene_capture_cells.csv")) as fh:
        for r in csv.DictReader(fh):
            if r["outcome"] == "pgfam_missing":
                miss[r["self_mode"]] += 1
    ok &= report("5 pgfam_missing == 0", not miss, "pgfam_missing cells %s" % dict(miss))
    ok &= report("5 gap close ran", os.path.exists(os.path.join(d, "close_exclude_rep_gap.done")), "")
    return ok


def check_baseline():
    r = subprocess.run(["sha256sum", "-c", "--quiet", BASELINE], cwd=R, capture_output=True, text=True)
    n = sum(1 for _ in open(BASELINE))
    return report("6 depth-20 baseline sha256 (%d files)" % n, r.returncode == 0,
                  (r.stdout + r.stderr).strip() or "all OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=["align", "copy", "post", "baseline", "all"])
    ap.add_argument("tree", nargs="?")
    a = ap.parse_args()
    ok = True
    if a.what == "align":
        ok = check_align(a.tree)
    elif a.what == "copy":
        ok = check_copy(a.tree)
    elif a.what == "post":
        ok = check_post(a.tree)
    elif a.what == "baseline":
        ok = check_baseline()
    else:
        for t in UNION:
            if os.path.exists(os.path.join(R, t, "data", "hits", "raw", "align.done")):
                print("== align %s" % t); ok &= check_align(t)
        for t in UNION:
            lt = t + "_legacy"
            if os.path.exists(os.path.join(R, lt, "data", "hits", "raw", "align.done")):
                print("== copy %s" % lt); ok &= check_copy(lt)
        for t in UNION + tuple(u + "_legacy" for u in UNION):
            if os.path.exists(os.path.join(R, t, "data", "gene_capture_by_region.json")):
                print("== post %s" % t); ok &= check_post(t)
        ok &= check_baseline()
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
