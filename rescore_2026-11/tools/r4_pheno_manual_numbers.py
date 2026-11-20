#!/usr/bin/env python
"""R4.1: the manuscript's phenotype numbers that no summary key stores (claims_spec 'manual' parts of PH-01..PH-05,
abs-phenotype-readout, disc.refusal_credit), recomputed with one recipe from the A-C13 ATP pass's per_genome.jsonl
and summary, OLD (pre-repair, moved aside) and NEW. The recipe is accepted only if OLD reproduces the printed values.

  n_v4         ATP-informative paired genomes (controlled, prFBA and PICRUSt2 both ok) whose amplicon is V4 (747)
  ecf1_delta   mean per-genome EC F1, prFBA exact-tie minus PICRUSt2, over those genomes (-0.028)
  media_signal media with ATP > 1e-6 in any row of truth / prFBA / PICRUSt2 / nn1 / family donor / constant (22 of 54)
  spearman     Spearman (EC-F1 delta, ATP-Hamming delta) and its p over the paired genomes (-0.14, p = 2e-4)
  mcnemar      media with q_bh <= 0.05 in mcnemar_per_medium, and which method they favour (2 of 22, PICRUSt2)
usage: r4_pheno_manual_numbers.py OUT_JSON
"""
import json
import sys

from scipy.stats import spearmanr

SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference/data"
PR = "/home/freiburger/Documents/prFBA_rescore_2026-11/pre_repair/region_validation/shared_reference/data"
MAIN = ("truth", "prfba_exact_tie", "picrust2", "base_nn1", "base_family_donor", "base_constant")
T = 1e-6


def one(d):
    S = json.load(open(f"{d}/phenotype_atp_c13/phenotype_summary.json"))
    rows = [json.loads(l) for l in open(f"{d}/phenotype_atp_c13/per_genome.jsonl") if l.strip()]
    by = {}
    for r in rows:
        by[(r["src_genome_id"], r["arm"], r["pipeline"])] = r
    genomes = sorted({r["src_genome_id"] for r in rows})
    media = sorted({m for r in rows for m in ((r.get("atp") or {}).get("vector") or {})})
    sig = sorted({m for r in rows if r["arm"] in MAIN for m, v in ((r.get("atp") or {}).get("vector") or {}).items()
                  if v is not None and v > T})
    pairs = []
    for g in genomes:
        t, p, q = (by.get((g, a, "controlled")) for a in ("truth", "prfba_exact_tie", "picrust2"))
        if not (t and p and q) or any(x.get("status") != "ok" for x in (t, p, q)):
            continue
        tv = (t.get("atp") or {}).get("vector") or {}
        if not any(v is not None and v > T for v in tv.values()):
            continue                              # informative-truth-only denominator
        hp, hq = (p.get("vs_truth") or {}).get("atp_hamming"), (q.get("vs_truth") or {}).get("atp_hamming")
        pairs.append({"g": g, "region": p.get("region"), "f1p": p.get("ec_f1"), "f1q": q.get("ec_f1"), "dh": hp - hq})
    f1 = [x for x in pairs if x["f1p"] is not None and x["f1q"] is not None]
    rho, pv = spearmanr([x["f1p"] - x["f1q"] for x in f1], [x["dh"] for x in f1])
    c = S["contrasts"]["prfba_exact_tie__vs__picrust2"]["controlled"]
    mc = c.get("mcnemar_per_medium") or {}
    sigm = {m: v for m, v in mc.items() if isinstance(v, dict) and v.get("q_bh") is not None and v["q_bh"] <= 0.05}
    return {"n_genomes": len(genomes), "n_media_in_vector": len(media), "media_with_signal": len(sig),
            "media_with_signal_list": sig, "n_pairs_informative": len(pairs),
            "n_pairs_summary": c.get("n_pairs"), "n_v4": sum(1 for x in pairs if x["region"] == "V4"),
            "ecf1_delta_mean": sum(x["f1p"] - x["f1q"] for x in f1) / len(f1), "n_ecf1": len(f1),
            "spearman_rho": rho, "spearman_p": pv, "spearman_rho_summary": c.get("ec_f1_delta_vs_phenotype_spearman"),
            "mcnemar_media_tested": len(mc), "mcnemar_q_le_0.05": {m: v for m, v in sigm.items()}}


out = {"old": one(PR), "new": one(SR)}
json.dump(out, open(sys.argv[1], "w"), indent=1, default=str)
for k in ("n_pairs_informative", "n_pairs_summary", "n_v4", "ecf1_delta_mean", "media_with_signal", "n_media_in_vector",
          "spearman_rho", "spearman_rho_summary", "spearman_p", "mcnemar_media_tested"):
    print(f"{k:24s} {out['old'][k]!s:>24} -> {out['new'][k]!s}")
print("mcnemar q<=0.05 old:", json.dumps(out["old"]["mcnemar_q_le_0.05"])[:400])
print("mcnemar q<=0.05 new:", json.dumps(out["new"]["mcnemar_q_le_0.05"])[:400])
