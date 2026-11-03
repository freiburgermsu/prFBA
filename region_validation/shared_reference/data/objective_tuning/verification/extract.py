"""Independent extraction of the Stage O rows into a compact pickle (no author code imported).

Reads: data/objective_tuning/per_amplicon.jsonl (streamed), data/weighted_objective/
per_amplicon.jsonl (prior run, streamed), data/weighted_objective/population.json (meta only).
Writes: rows.pkl, prior.pkl in this directory.
"""
import json
import os
import pickle
import sys

SR = "/home/freiburger/Documents/prFBA/region_validation/shared_reference"
OUT = os.path.dirname(os.path.abspath(__file__))

ARM_FIELDS = ("status", "objective", "biomass", "biomass_ratio", "max_biomass",
              "taf_face_min", "taf_face_mid", "taf_face_max", "taf_face_width",
              "absent_flux_min", "absent_flux_max", "taf_gb", "taf_point_D",
              "flux_exchange", "flux_gene_backed", "flux_absent", "flux_gapfilled",
              "flux_rxn", "precision_flux", "gapfilled_flux_face_min",
              "gapfilled_flux_face_max", "internal_term", "exchange_term", "internal_share",
              "internal_term_split", "n_active_gb", "n_active_absent", "n_active_gapfilled",
              "face_lp_failed", "face_lp_retried", "g_value", "x_value", "factors",
              "futile_offset_bound")
TOP = ("key", "src_genome_id", "region", "md5", "fold", "group_id", "template", "genus",
       "n_reps", "stratum", "reps", "novelty_bin", "rank1_identity", "tier", "half", "status",
       "reason", "gapfill", "model", "prior_match", "est_info", "determinism",
       "released_check", "shipped_check")


def main():
    rows = []
    with open(os.path.join(SR, "data/objective_tuning/per_amplicon.jsonl")) as fh:
        for line in fh:
            r = json.loads(line)
            c = {k: r.get(k) for k in TOP}
            c["arms"] = {a: {k: v.get(k) for k in ARM_FIELDS} for a, v in
                         (r.get("arms") or {}).items()}
            c["n_arm_keys"] = {a: len(v) for a, v in (r.get("arms") or {}).items()}
            rows.append(c)
            del r
    pickle.dump(rows, open(os.path.join(OUT, "rows.pkl"), "wb"))
    prior = []
    with open(os.path.join(SR, "data/weighted_objective/per_amplicon.jsonl")) as fh:
        for line in fh:
            r = json.loads(line)
            c = {k: r.get(k) for k in ("key", "group_id", "stratum", "status", "gapfill",
                                       "model", "src_genome_id", "genus", "md5", "fold")}
            c["arms"] = {a: {k: v.get(k) for k in ("objective", "taf_face_min",
                                                   "taf_face_mid", "taf_face_max", "taf_gb",
                                                   "flux_gene_backed", "flux_absent",
                                                   "status")}
                         for a, v in (r.get("arms") or {}).items()}
            prior.append(c)
            del r
    pickle.dump(prior, open(os.path.join(OUT, "prior.pkl"), "wb"))
    print(len(rows), "rows;", len(prior), "prior rows")


if __name__ == "__main__":
    sys.exit(main())
