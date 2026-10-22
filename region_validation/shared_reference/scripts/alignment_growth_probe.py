#!/usr/bin/env python
"""Which database-direction restrictions stop a gapfilled model from growing?

Takes a model rebuilt with the five-reaction template fix and the gapfiller's reversals
integrated (data/egc_analysis/models_tv_PR1v2_rev, or --models DIR). Every reaction in it
that the complete database alignment restricts (tiers B and C in
data/template_direction_audit/direction_audit.tsv) is restricted to the database
direction, and growth on the complete gapfilling medium is reported with all
restrictions applied and with each one alone. A reaction that stops growth by itself is
one the model needs in the database-forbidden direction.

    python alignment_growth_probe.py [--models DIR] MODEL_KEY [MODEL_KEY ...] \
        > data/template_direction_audit/alignment_growth_probe.json

Prints one JSON list, one record per model.
"""
import argparse
import csv
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import phenotype_fold as PF  # noqa: E402
import cobra  # noqa: E402

SR = os.path.abspath(os.path.join(HERE, ".."))
TEMPLATE_FILE = {"archaea": "ArchaeaTemplateV6.json", "gramneg": "GramNegModelTemplateV7.json",
                 "grampos": "GramPosModelTemplateV7.json"}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--models", default=os.path.join(SR, "data/egc_analysis/models_tv_PR1v2_rev"))
    ap.add_argument("keys", nargs="+")
    a = ap.parse_args()
    m = PF.pa()
    _, gf_cpds = m.complete_media()
    audit = list(csv.DictReader(open(os.path.join(SR, "data/template_direction_audit/direction_audit.tsv")),
                                delimiter="\t"))

    def grow(md):
        return float(m.growth_vector(md, {"complete": gf_cpds}, thresh=PF.GROWTH_THRESH).get("complete", 0.0))

    records = []
    for key in a.keys:
        base = os.path.join(a.models, key)
        st = json.load(open(base + ".stages.json"))
        tfile = TEMPLATE_FILE[st["template"]]
        tgt = {r["reaction"]: r["target"] for r in audit
               if r["template"] == tfile and r["tier"] in ("B", "C")}
        mdl = cobra.io.load_json_model(base + ".model.json")
        hit = [r.id for r in mdl.reactions if r.id[:-1] in tgt]   # rxnNNNNN_c0 -> rxnNNNNN_c

        def restrict(md, ids):
            for rid in ids:
                r = md.reactions.get_by_id(rid)
                if tgt[rid[:-1]] == ">":
                    r.lower_bound = max(r.lower_bound, 0.0)
                else:
                    r.upper_bound = min(r.upper_bound, 0.0)

        g0 = grow(mdl)
        with mdl:
            restrict(mdl, hit)
            g_all = grow(mdl)
        killers = []
        for rid in hit:
            r = mdl.reactions.get_by_id(rid)
            as_built = r.reaction   # bounds as gapfilled, before the restriction
            with mdl:
                restrict(mdl, [rid])
                stops = grow(mdl) <= PF.GROWTH_THRESH
            if stops:
                killers.append({"reaction": rid[:-3], "name": r.name.replace(" [c0]", ""),
                                "database_direction": tgt[rid[:-1]], "model_reaction": as_built})
        records.append({"model": key, "models_dir": os.path.relpath(a.models, SR), "template": tfile,
                        "growth": round(g0, 4),
                        "n_restricted": len(hit), "growth_all_restricted": round(g_all, 4),
                        "single_restriction_stops_growth": killers})
    print(json.dumps(records, indent=1))


if __name__ == "__main__":
    main()
