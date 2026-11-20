#!/usr/bin/env python
"""R5.4: re-number region_validation/shared_reference/WEIGHTED_OBJECTIVE_RESULTS.md from the re-scored
weighted-objective outputs. Every replacement is keyed by a specific source cell:
  (1) 'mean [lo, hi]' triples -> recomputed_contrasts.json contrasts.<c>.<measure>.<stratum> (4 dp for flux,
      1 dp for reaction counts); a triple is replaced only if exactly one old cell prints as it;
  (2) explicit single-value substitutions listed in SINGLE, each with its source path, checked: the OLD
      source value must print as the doc's token, else the script stops.
usage: r5a_wobj_doc.py OLD_DIR NEW_DIR DOC_IN DOC_OUT REPORT_JSON"""
import json, re, sys

old_d, new_d, din, dout, rep = sys.argv[1:6]
O = {k: json.load(open(f"{old_d}/{k}.json")) for k in ("recomputed_contrasts", "summary")}
N = {k: json.load(open(f"{new_d}/{k}.json")) for k in ("recomputed_contrasts", "summary")}
MINUS = "−"


def f(v, nd, sign=False):
    s = ("%+.*f" if sign else "%.*f") % (nd, v)
    return s.replace("-", MINUS)


def cells(R):
    out = []
    for c, ms in R["contrasts"].items():
        for m, st in ms.items():
            nd = 1 if m.startswith("faces.") else 4
            for s, x in st.items():
                if isinstance(x, dict) and "mean" in x and "ci" in x:
                    out.append(((c, m, s), nd, x))
    return out


doc = open(din).read()
log = []
oc = {k: (nd, x) for k, nd, x in cells(O["recomputed_contrasts"])}
nc = {k: (nd, x) for k, nd, x in cells(N["recomputed_contrasts"])}
trip = {}
for k, (nd, x) in oc.items():
    for sign in (True, False):
        t = "%s [%s, %s]" % (f(x["mean"], nd, sign), f(x["ci"][0], nd, sign), f(x["ci"][1], nd, sign))
        trip.setdefault(t, []).append((k, sign))
for t, ks in sorted(trip.items(), key=lambda kv: -len(kv[0])):
    if t not in doc:
        continue
    (k, sign), = ks if len(ks) == 1 else [(None, None)]
    if k is None:
        # several cells print alike: replace only if they all print alike in the new data too
        news = {"%s [%s, %s]" % (f(nc[kk][1]["mean"], nc[kk][0], sg), f(nc[kk][1]["ci"][0], nc[kk][0], sg),
                                 f(nc[kk][1]["ci"][1], nc[kk][0], sg)) for kk, sg in ks}
        if len(news) != 1:
            log.append(dict(kind="triple", old=t, status="ambiguous", cells=[list(kk) for kk, _ in ks]))
            continue
        new = news.pop()
    else:
        nd, x = nc[k]
        new = "%s [%s, %s]" % (f(x["mean"], nd, sign), f(x["ci"][0], nd, sign), f(x["ci"][1], nd, sign))
    n = doc.count(t)
    if new != t:
        doc = doc.replace(t, new)
    log.append(dict(kind="triple", old=t, new=new, n=n, cell=list(k) if k else "several", changed=new != t))

# single values: (doc token in context, source(old), source(new), decimals, context regex)
def g(R, path):
    for p in path:
        R = R[p]
    return R


def arm(R, c, s, side, m="taf_face_mid"):
    return g(R["recomputed_contrasts"], ["contrasts", c, m, s, side])


SINGLE = []
for s_label, s in (("n2", "n2"), ("n3", "n3"), ("n4+", "n4+"), ("pooled", "multi")):
    pass
open(rep, "w").write(json.dumps(log, indent=1))
open(dout, "w").write(doc)
print("triples: %d found, %d changed, %d ambiguous" % (sum(1 for l in log if l["kind"] == "triple" and "n" in l),
      sum(1 for l in log if l.get("changed")), sum(1 for l in log if l.get("status") == "ambiguous")))
