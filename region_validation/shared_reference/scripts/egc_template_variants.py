#!/usr/bin/env python
"""Template variants for the draft-level scope scan: ModelSEEDTemplates 21a8084 (the
five-reaction fix) with more reactions restricted to their forward direction, edited the way
21a8084 edits a '=' reaction (direction and GapfillDirection '>', lower_bound 0,
maxrevflux 0). These are probes of where the raw drafts' ATP-from-nothing loops come from,
not proposed templates: nothing here was validated through gapfilling.

  ccp       rxn00014 (cytochrome c peroxidase): reversed, it closes the only draft loop the
            five-reaction fix opens (with catalase forward, cbb3 oxidase and F1-ATPase)
  ccp_gapn_mdg
            ccp, plus rxn23850 (GAPN, the first support of 130 of the 136 leaking truth
            drafts) and the maltodextrin glucosidases rxn05746, rxn08933-rxn08936

    python egc_template_variants.py --base ~/Documents/ModelSEEDTemplates-c12run/templates/v7.0 \
        --out data/egc_analysis/templates_scope
"""
import argparse
import hashlib
import json
import os

VARIANTS = {
    "ccp": ["rxn00014"],
    "ccp_gapn_mdg": ["rxn00014", "rxn23850", "rxn05746", "rxn08933", "rxn08934",
                     "rxn08935", "rxn08936"],
}
FILES = ("GramNegModelTemplateV7.json", "GramPosModelTemplateV7.json",
         "ArchaeaTemplateV6.json", "Core-V6.json")


def forward_only(rxn):
    rxn["direction"] = ">"
    rxn["GapfillDirection"] = ">"
    rxn["lower_bound"] = 0
    rxn["maxrevflux"] = 0.0


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--base", required=True, help="templates/v7.0 directory to start from")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    manifest = {"base": os.path.abspath(args.base), "variants": {}}
    for name, rids in VARIANTS.items():
        d = os.path.join(args.out, name)
        os.makedirs(d, exist_ok=True)
        rec = {"reactions": rids, "edited": {}, "sha256": {}}
        for f in FILES:
            t = json.load(open(os.path.join(args.base, f)))
            hit = []
            for r in t["reactions"]:
                if r["id"].rsplit("_", 1)[0] in rids:
                    forward_only(r)
                    hit.append(r["id"])
            path = os.path.join(d, f)
            with open(path, "w") as fh:
                json.dump(t, fh)
            rec["edited"][f] = sorted(hit)
            rec["sha256"][f] = hashlib.sha256(open(path, "rb").read()).hexdigest()
        manifest["variants"][name] = rec
        print(name, {f: len(v) for f, v in rec["edited"].items()})
    with open(os.path.join(args.out, "manifest.json"), "w") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)


if __name__ == "__main__":
    main()
