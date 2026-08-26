"""Rigorous correctness checks for query.py (the cosine-similarity reference tool).

Fixtures (``/tmp/prFBA_test`` + ``/tmp/sample_fasta.json``) are built by
``tests/make_fixtures.py`` and auto-generated below if absent.
"""
import os, sys, json, subprocess, hashlib, numpy as np, pandas as pd, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import nt_embed

if not (os.path.exists("/tmp/prFBA_test/embeddings.f16.npy")
        and os.path.exists("/tmp/sample_fasta.json")):
    import make_fixtures
    make_fixtures.main()

STORE = "/tmp/prFBA_test"
PY = "/home/freiburger/Documents/py_venv/bin/python"
QP = "/home/freiburger/Documents/prFBA/query.py"
fails = []
def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok: fails.append(name)

# ---- build queries: exact-in-store (full 16S), 1bp variant, short frag, random novel ----
fa = json.load(open("/tmp/sample_fasta.json"))
full = [v.upper() for v in fa.values() if 1400 <= len(v) <= 1600]
exact = full[0]
variant = exact[:700] + ("A" if exact[700] != "A" else "C") + exact[701:]
short = [v.upper() for v in fa.values() if 200 <= len(v) <= 350][0]
novel = "ACGT" * 375
queries = {"exact": exact, "variant": variant, "short": short, "novel": novel}
with open("/tmp/qv.fasta", "w") as f:
    for k, v in queries.items():
        f.write(f">{k}\n{v}\n")
json.dump(queries, open("/tmp/qv.json", "w"))

# ---- independent brute-force reference ----
idx = pd.read_parquet(f"{STORE}/index.parquet")
E = np.load(f"{STORE}/embeddings.f16.npy").astype(np.float32)
tok, model = nt_embed.load_model()
qv = nt_embed.embed_batch(list(queries.values()), tok, model).float().cpu().numpy()  # (4,1024)
S = qv @ E.T                                                                          # (4, N) cosine
bf_top = {name: np.argsort(-S[i])[:5] for i, name in enumerate(queries)}

# ---- run the real CLI ----
out = subprocess.run([PY, QP, "--fasta", "/tmp/qv.fasta", "--store", STORE,
                      "--topk", "5", "--out", "/tmp/qv_hits.csv"],
                     capture_output=True, text=True)
hits = pd.read_csv("/tmp/qv_hits.csv")

print("# 1. reported cosine == independent brute force, ordering, top-k membership")
for i, name in enumerate(queries):
    h = hits[hits["query"] == name].sort_values("rank")
    # cosine values match brute force at the reported row_ids
    bf_at_rows = S[i][h["row_id"].values]
    cos_ok = np.max(np.abs(h["cosine"].values - bf_at_rows)) < 3e-3
    desc_ok = np.all(np.diff(h["cosine"].values) <= 1e-6)
    # top-1 row matches brute-force argmax (allow fp ties at near-identical sims)
    top1_ok = (h["row_id"].iloc[0] == bf_top[name][0]) or abs(S[i][h["row_id"].iloc[0]] - S[i][bf_top[name][0]]) < 3e-3
    check(f"{name}: cosine==brute-force", cos_ok, f"maxΔ={np.max(np.abs(h['cosine'].values-bf_at_rows)):.1e}")
    check(f"{name}: ranks descending", desc_ok)
    check(f"{name}: top-1 == brute-force argmax", top1_ok)

print("\n# 2. exact-match semantics")
he = hits[hits["query"] == "exact"].sort_values("rank").iloc[0]
emd5 = hashlib.md5(exact.encode()).hexdigest()
check("exact: cosine ~1.0", he["cosine"] > 0.999, f"cos={he['cosine']:.4f}")
check("exact: md5 flag True", bool(he["exact_md5_match"]))
check("exact: matched row has that md5", idx.iloc[int(he["row_id"])]["md5"] == emd5)
hn = hits[hits["query"] == "novel"].sort_values("rank").iloc[0]
check("novel: flagged sparse (top cos < 0.8)", hn["cosine"] < 0.8, f"cos={hn['cosine']:.4f}")

print("\n# 3. CLI behaviors: batch==single, min-sim filter, JSON input parsing")
# single-seq invocation must reproduce the batch row for 'exact'
o2 = subprocess.run([PY, QP, "--seq", exact, "--store", STORE, "--topk", "5", "--out", "/tmp/qv_single.csv"],
                    capture_output=True, text=True)
single = pd.read_csv("/tmp/qv_single.csv").sort_values("rank")
batch_exact = hits[hits["query"] == "exact"].sort_values("rank")
# Invariant: batching does not change the SCORES. The encoder's fp16 mean-pool is
# not bit-identical across batch shapes (METHODS Section 5), so sub-3e-3 cosine
# ties at deep ranks can reorder harmlessly; we therefore require an identical
# top-1 row and an identical cosine spectrum (within tolerance), not a brittle
# exact ordering of near-tied ranks.
top1_ok = int(single["row_id"].iloc[0]) == int(batch_exact["row_id"].iloc[0])
cos_ok = np.max(np.abs(np.sort(single["cosine"].values)[::-1]
                       - np.sort(batch_exact["cosine"].values)[::-1])) < 3e-3
check("batch == single (top-1 + cosine spectrum within 3e-3)", top1_ok and cos_ok,
      f"top1_match={top1_ok}")
# min-sim filter
o3 = subprocess.run([PY, QP, "--fasta", "/tmp/qv.fasta", "--store", STORE, "--topk", "5",
                     "--min-sim", "0.95", "--out", "/tmp/qv_min.csv"], capture_output=True, text=True)
hmin = pd.read_csv("/tmp/qv_min.csv")
check("min-sim filter: all hits >= 0.95", bool((hmin["cosine"] >= 0.95).all()), f"min={hmin['cosine'].min():.3f}")
# JSON input parses to same as FASTA
o4 = subprocess.run([PY, QP, "--fasta", "/tmp/qv.json", "--store", STORE, "--topk", "5", "--out", "/tmp/qv_json.csv"],
                    capture_output=True, text=True)
hjson = pd.read_csv("/tmp/qv_json.csv")
same = np.array_equal(hjson.sort_values(["query","rank"])["row_id"].values,
                      hits.sort_values(["query","rank"])["row_id"].values)
check("JSON input == FASTA input", same)

print("\n# 4. encoder determinism")
a = nt_embed.embed_batch([exact], tok, model).float().cpu().numpy()
b = nt_embed.embed_batch([exact], tok, model).float().cpu().numpy()
cos_ab = float((a @ b.T)[0, 0])
check("same sequence -> identical vector (cos==1)", cos_ab > 0.99999, f"cos={cos_ab:.6f}")

print("\nRESULT:", "ALL PASS" if not fails else f"{len(fails)} FAILED: {fails}")
sys.exit(1 if fails else 0)
