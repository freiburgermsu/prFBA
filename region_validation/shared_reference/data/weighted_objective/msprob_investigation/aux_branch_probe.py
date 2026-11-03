#!/usr/bin/env python
"""aux_branch_probe.py: does MSProbability.prFBA run on another branch's tree?

Extract the branch read-only (no checkout) and run this with that tree first on PYTHONPATH:

  mkdir -p /tmp/br && git -C ~/Documents/ModelSEEDpy archive <ref> modelseedpy | tar -x -C /tmp/br
  PYTHONPATH=/tmp/br ~/Documents/py_venv/bin/python -W ignore aux_branch_probe.py

Used for cshenry/main (b758b44; live tip ac02db7 on 2026-11-02 has a byte-identical
msprobability.py, last changed at 13274cc, 2024-03-17) and freiburgermsu origin/main (971df5d).
Logs: aux_branch_probe_cshenry_main.log, aux_branch_probe_origin_main.log.

cshenry/main carries a prFBA with the `if environment is not None` guard and a single
ElementUptakePkg build, but its own MSModelUtil has a different API (no add_medium,
add_minimal_objective_cons or add_objective; __init__(self, model) only). origin/main carries the
c2e31a0 defects plus an import of the deleted modelseedpy.community.commhelper.
"""
import contextlib
import gzip
import io
import os
import traceback

import modelseedpy
import cobra
import cobra.io

HERE = os.path.dirname(os.path.abspath(__file__))
print("modelseedpy from", modelseedpy.__file__)
try:
    from modelseedpy.core.msprobability import MSProbability  # noqa: E402
except Exception as exc:                                       # noqa: BLE001
    print("import of modelseedpy.core.msprobability FAILED:", type(exc).__name__, exc)
    raise SystemExit(0)
from modelseedpy.core.msmodelutl import MSModelUtil        # noqa: E402

print("MSModelUtil.__init__ signature:", MSModelUtil.__init__.__code__.co_varnames[:MSModelUtil.__init__.__code__.co_argcount])
print("MSModelUtil has add_medium / add_minimal_objective_cons / add_objective:",
      [hasattr(MSModelUtil, a) for a in ("add_medium", "add_minimal_objective_cons", "add_objective")])


def load():
    with gzip.open(os.path.join(HERE, "fixture_g100053.8_gf.json.gz"), "rt") as fh:
        m = cobra.io.load_json_model(fh)
    m.objective = "bio1"
    for r in m.reactions:
        if r.id.startswith("rxn"):
            r.notes["probability"] = str(1.0 if r.genes else 0.0)
    return m


def member(mid, carbon, uptake):
    form = {carbon: "C6H12O6", "cpd00002": "C10H13N5O13P3", "cpd00001": "H2O",
            "cpd00008": "C10H13N5O10P2", "cpd00009": "HO4P", "cpd00067": "H", "cpd11416": "R"}
    m = cobra.Model(mid)
    mets = {k: cobra.Metabolite(k, name=k, compartment=k[-2:], formula=form[k[:-3]])
            for k in (f"{carbon}_e0", f"{carbon}_c0", "cpd00002_c0", "cpd00001_c0",
                      "cpd00008_c0", "cpd00009_c0", "cpd00067_c0", "cpd11416_c0")}
    m.add_metabolites(list(mets.values()))
    ex = cobra.Reaction(f"EX_{carbon}_e0", lower_bound=-float(uptake), upper_bound=1000)
    ex.add_metabolites({mets[f"{carbon}_e0"]: -1})
    tr = cobra.Reaction("rxn10000_c0", lower_bound=0, upper_bound=1000)
    tr.add_metabolites({mets[f"{carbon}_e0"]: -1, mets[f"{carbon}_c0"]: 1})
    tr.notes["probability"] = "0.5"
    bio = cobra.Reaction("bio1", lower_bound=0, upper_bound=1000)
    bio.add_metabolites({mets[f"{carbon}_c0"]: -1, mets["cpd11416_c0"]: 1})
    sk = cobra.Reaction("SK_cpd11416_c0", lower_bound=0, upper_bound=1000)
    sk.add_metabolites({mets["cpd11416_c0"]: -1})
    m.add_reactions([ex, tr, bio, sk])
    m.objective = "bio1"
    return m


ENV = {c.strip(): 100.0 for c in open(os.path.join(HERE, "fixture_kbase_complete_media.txt")) if c.strip()}
for label, fn in (("prFBA(model, environment=ENV)", lambda: MSProbability.prFBA(load(), environment=ENV)),
                  ("prFBA([model], environment=ENV)", lambda: MSProbability.prFBA([load()], environment=ENV)),
                  ("prFBA([memA, memB], environment=...)",
                   lambda: MSProbability.prFBA([member("memA", "cpd00027", 10), member("memB", "cpd00082", 4)],
                                               environment={"cpd00027": 10, "cpd00082": 4}))):
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            fn()
        print("{:<40s} returned".format(label))
    except Exception as exc:                        # noqa: BLE001
        frames = traceback.extract_tb(exc.__traceback__)
        where = " <- ".join("{}:{}".format(f.filename.split("modelseedpy/")[-1] if "modelseedpy/" in f.filename
                                           else os.path.basename(f.filename), f.lineno) for f in frames[::-1][:3])
        print("{:<40s} RAISED {}: {} | at {}".format(label, type(exc).__name__, str(exc)[:120], where))
