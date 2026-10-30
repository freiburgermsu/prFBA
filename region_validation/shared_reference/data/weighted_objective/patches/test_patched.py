"""Run the PATCHED msprobability.py (scratch copy, both patches) on the probe model."""
import sys, json, io, contextlib, importlib.util
import cobra.io
spec = importlib.util.spec_from_file_location("msprob_patched", sys.argv[1])
M = importlib.util.module_from_spec(spec); spec.loader.exec_module(M)
path = sys.argv[2]
meta = json.load(open(path + ".meta.json")); new = set(meta["block"]["new_rxn"])
env = {c.strip(): 100.0 for c in open("/home/freiburger/Documents/ModelSEEDDatabase/Media/KBaseMedia.cpd") if c.strip()}
def fresh(strip_gf=False):
    m = cobra.io.load_json_model(path); m.objective = "bio1"
    for r in m.reactions:
        if strip_gf and r.id[:-1] in new: continue
        r.notes["probability"] = str(1.0 if r.genes else 0.0)
    return m
for label, fn in [("bare Model", lambda: M.MSProbability.prFBA(fresh(), environment=env)),
                  ("[Model]", lambda: M.MSProbability.prFBA([fresh()], environment=env)),
                  ("bare Model, gapfilled rxns un-noted", lambda: M.MSProbability.prFBA(fresh(True), environment=env)),
                  ("bare Model, environment=None", lambda: M.MSProbability.prFBA(fresh()))]:
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        s = fn()
    print("%-40s status=%s obj=%.6f bio1=%.6f printed_chars=%d" % (label, s.status, s.objective_value, s.fluxes["bio1"], len(buf.getvalue())))
