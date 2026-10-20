#!/usr/bin/env python
"""Build the shared functional namespaces that let prFBA and PICRUSt2 be compared.

prFBA predicts a genome's PGFam gene families; PICRUSt2 predicts EC numbers (and
KOs) per sequence. The only vocabulary both sides can be expressed in without
inventing a mapping is the EC number, and from EC the ModelSEED reaction space
that prFBA's synthetic genomes are ultimately reconstructed in:

  PGFam --family_product--> EC numbers            (parsed from the product string)
  PGFam --family_product--> ModelSEED roles --> complexes --> reactions  (templates)
  EC ------------------------------------------> ModelSEED reactions  (Aliases)

Role handling follows ModelSEED/KBase reconstruction: ``split_role`` splits
multi-function products, ``convert_to_search_role`` normalises each part, and a
reaction is reachable if any triggering role of any of its complexes is present
(ModelSEEDpy's permissive complex rule, so prFBA is never under-credited).

Outputs ``--out`` (JSON):
  pgfam_ec      {PGF_id: [EC, ...]}          complete 4-level ECs only
  pgfam_rxn     {PGF_id: [rxn, ...]}         via ModelSEED template roles
  ec_rxn        {EC: [rxn, ...]}             ModelSEED reaction EC aliases
  role_hits     coverage stats for the mapping itself
"""
import argparse
import glob
import json
import os
import re
import sys

import orjson

MSD = "/home/freiburger/Documents/ModelSEEDDatabase"
MST = "/home/freiburger/Documents/ModelSEEDTemplates/templates/v7.0"
TEMPLATES = ["GramNegModelTemplateV7.json", "GramPosModelTemplateV7.json",
             "ArchaeaTemplateV6.json"]

# EC in a BV-BRC/RAST product string: "(EC 1.2.3.4)" / "(EC: 1.2.3.4)".
EC_RE = re.compile(r"\(EC:?\s*([0-9]+\.[0-9]+\.[0-9]+\.[0-9]+)\)", re.I)
COMPLETE_EC_RE = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+\.[0-9]+$")


def split_role(role):
    """ModelSEEDpy annotation_ontology.split_role (multi-function annotations)."""
    return re.split(r"\s*;\s+|\s+[@/]\s+", role)


def convert_to_search_role(role):
    """ModelSEEDpy annotation_ontology.convert_to_search_role."""
    role = role.lower()
    role = re.sub(r"\s", "", role)
    role = re.sub(r"[\d\-]+\.[\d\-]+\.[\d\-]+\.[\d\-]*", "", role)
    role = re.sub(r"\#.*$", "", role)
    role = re.sub(r"\(ec:*\)", "", role)
    role = re.sub(r"[\(\)\[\],-]", "", role)
    return role


def load_templates():
    """search_role -> {base rxn ids}, from role -> complex -> template reaction."""
    role_rxn = {}
    stats = {"templates": [], "n_roles": 0, "n_rxn_with_complex": 0}
    for name in TEMPLATES:
        t = json.load(open(os.path.join(MST, name)))
        roles = {r["id"]: r["name"] for r in t["roles"]}
        cpx_roles = {}
        for c in t["complexes"]:
            ids = [cr["templaterole_ref"].rsplit("/", 1)[-1] for cr in c["complexroles"]
                   if cr.get("triggering", 1)]
            cpx_roles[c["id"]] = ids
        n_rxn = 0
        for rxn in t["reactions"]:
            refs = rxn.get("templatecomplex_refs") or []
            if not refs:
                continue
            n_rxn += 1
            base = rxn["id"].rsplit("_", 1)[0]
            for ref in refs:
                for rid in cpx_roles.get(ref.rsplit("/", 1)[-1], []):
                    search = convert_to_search_role(roles.get(rid, ""))
                    if search:
                        role_rxn.setdefault(search, set()).add(base)
        stats["templates"].append({"file": name, "n_rxn_with_complex": n_rxn})
        stats["n_roles"] = len(roles)
        stats["n_rxn_with_complex"] = max(stats["n_rxn_with_complex"], n_rxn)
    return role_rxn, stats


def load_ec_rxn():
    """EC -> {rxn}: ModelSEED reaction EC aliases (obsolete reactions dropped)."""
    obsolete = set()
    for path in glob.glob(os.path.join(MSD, "Biochemistry", "reaction_*.tsv")):
        with open(path) as fh:
            header = fh.readline().rstrip("\n").split("\t")
            i_id, i_obs = header.index("id"), header.index("is_obsolete")
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if f[i_obs] == "1":
                    obsolete.add(f[i_id])
    ec_rxn = {}
    with open(os.path.join(MSD, "Biochemistry", "Aliases",
                          "Unique_ModelSEED_Reaction_ECs.txt")) as fh:
        fh.readline()
        for line in fh:
            rxn, ec, _src = line.rstrip("\n").split("\t")[:3]
            if COMPLETE_EC_RE.match(ec) and rxn not in obsolete:
                ec_rxn.setdefault(ec, set()).add(rxn)
    return ec_rxn, len(obsolete)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    here = os.path.dirname(os.path.abspath(__file__))
    ap.add_argument("--products", default=os.path.join(here, "..", "data",
                                                       "pgfam_products.json"))
    ap.add_argument("--out", default=os.path.join(here, "..", "data", "function_maps.json"))
    args = ap.parse_args()

    role_rxn, tstats = load_templates()
    ec_rxn, n_obsolete = load_ec_rxn()
    products = orjson.loads(open(args.products, "rb").read())
    print(f"[function_maps] {len(products)} PGFams, {len(role_rxn)} template search-roles, "
          f"{len(ec_rxn)} ECs with reactions ({n_obsolete} obsolete reactions dropped)")

    pgfam_ec, pgfam_rxn = {}, {}
    n_prod = n_ec = n_role = n_hypo = 0
    for pgf, rec in products.items():
        prod = (rec or {}).get("family_product")
        if not prod:
            continue
        n_prod += 1
        if "hypothetical protein" == prod.strip().lower():
            n_hypo += 1
        ecs = sorted({e for e in EC_RE.findall(prod)})
        if ecs:
            pgfam_ec[pgf] = ecs
            n_ec += 1
        rxns = set()
        for part in split_role(prod):
            search = convert_to_search_role(part)
            if search:
                rxns |= role_rxn.get(search, set())
        if rxns:
            pgfam_rxn[pgf] = sorted(rxns)
            n_role += 1

    out = {
        "pgfam_ec": pgfam_ec,
        "pgfam_rxn": pgfam_rxn,
        "ec_rxn": {k: sorted(v) for k, v in ec_rxn.items()},
        "_meta": {
            "n_pgfams": len(products), "n_with_product": n_prod,
            "n_hypothetical": n_hypo,
            "n_pgfam_with_ec": n_ec, "n_pgfam_with_template_rxn": n_role,
            "n_distinct_ec": len({e for v in pgfam_ec.values() for e in v}),
            "n_distinct_rxn_role": len({r for v in pgfam_rxn.values() for r in v}),
            "n_ec_with_rxn": len(ec_rxn),
            "modelseed_db": MSD, "modelseed_templates": tstats,
            "role_rule": "any triggering role of any complex of the reaction",
        },
    }
    with open(args.out, "wb") as fh:
        fh.write(orjson.dumps(out))
    print(f"[function_maps] -> {args.out}: {json.dumps(out['_meta'])[:600]}")


if __name__ == "__main__":
    main()
