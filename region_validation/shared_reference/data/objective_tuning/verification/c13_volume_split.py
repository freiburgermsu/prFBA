"""Exact per-row split of the change in truth-absent gene-backed flux (units of D, returned
vertex) into a volume term s_b*(gb_a-gb_b) and a selectivity term gb_a*(s_a-s_b), where
s = flux_absent/flux_gene_backed; pooled A with the cluster bootstrap."""
import json, os
import vlib as L
rows = L.load_rows(); WA, WB, _ = L.weights(); P = {"A": WA, "B": WB}
conf = [r for r in rows if r["half"] == "confirm"]
def parts(r, a, b):
    D = r["model"]["D_flux_gb_two_tier"]
    ga, gb_ = L.v(r, a, "flux_gene_backed") / D, L.v(r, b, "flux_gene_backed") / D
    sa, sb = L.v(r, a, "flux_absent") / L.v(r, a, "flux_gene_backed"), L.v(r, b, "flux_absent") / L.v(r, b, "flux_gene_backed")
    return sb * (ga - gb_), ga * (sa - sb)
out = {}
for name, (a, b) in {"selected-shipped": (L.SELECTED, L.SHIPPED), "est@X100-shipped": ("est/lin/G0/X100", L.SHIPPED),
                     "est X1-X100": (L.SELECTED, "est/lin/G0/X100"), "kn X1-X100": ("kn/lin/G0/X1", L.SHIPPED)}.items():
    vol = L.paired(conf, lambda r: parts(r, a, b)[0], None, pools=P)["pooled"]
    sel = L.paired(conf, lambda r: parts(r, a, b)[1], None, pools=P)["pooled"]
    out[name] = {"volume_A": vol["A"], "selectivity_A": sel["A"], "volume_B": vol["B"], "selectivity_B": sel["B"]}
    print("%-18s volume A %s | selectivity A %s | volume B %s | selectivity B %s" % (name, L.fmt(vol["A"]), L.fmt(sel["A"]), L.fmt(vol["B"]), L.fmt(sel["B"])))
json.dump(out, open(os.path.join(L.HERE, "out_c13.json"), "w"), indent=1)
