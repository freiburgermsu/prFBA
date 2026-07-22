"""Rigorous correctness checks for analyze_space.py internals.

Fixtures (``/tmp/prFBA_test``) are built by ``tests/make_fixtures.py`` and
auto-generated below if absent, so this runs from a clean checkout.
"""
import os, sys, numpy as np, torch
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, "/home/freiburger/Documents/prFBA")
import analyze_space as A

if not os.path.exists("/tmp/prFBA_test/embeddings.f16.npy"):
    import make_fixtures
    make_fixtures.main()

rng = np.random.default_rng(0)
fails = []
def check(name, ok, detail=""):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name} {detail}")
    if not ok: fails.append(name)

# ---------------------------------------------------------------- kNN vs brute force
print("\n# 1. GPU kNN vs CPU brute force (real 691x1024 store)")
E = np.load("/tmp/prFBA_test/embeddings.f16.npy").astype(np.float32)
E /= np.linalg.norm(E, axis=1, keepdims=True)
Eg = torch.from_numpy(E.astype(np.float16)).cuda()
k = 10; radii = [0.99, 0.97, 0.95, 0.90, 0.80]
nn1, nnk, nn1i, rad = A.gpu_knn(Eg, k, chunk=128, radii=radii)
# brute force in float32
S = E @ E.T
np.fill_diagonal(S, -np.inf)
bf_nn1 = S.max(1)
bf_topk = np.sort(S, axis=1)[:, ::-1][:, :k]
bf_nnk = bf_topk.mean(1)
Sr = E @ E.T; np.fill_diagonal(Sr, -1.0)
bf_rad = np.stack([(Sr >= r).sum(1) for r in radii], 1).astype(float)
check("NN1 matches brute force", np.max(np.abs(nn1 - bf_nn1)) < 3e-3, f"maxdiff={np.max(np.abs(nn1-bf_nn1)):.1e}")
check("mean-kNN matches brute force", np.max(np.abs(nnk - bf_nnk)) < 3e-3, f"maxdiff={np.max(np.abs(nnk-bf_nnk)):.1e}")
check("radius counts ~match (mean within 0.5)", np.abs(rad.mean(0) - bf_rad.mean(0)).max() < 0.5,
      f"gpu={rad.mean(0).round(2)} bf={bf_rad.mean(0).round(2)}")

# ---------------------------------------------------------------- chunk invariance
print("\n# 2. Chunk invariance (auto-chunk stable to fp16 matmul-tiling noise)")
nn1b, nnkb, _, radb = A.gpu_knn(Eg, k, chunk=691, radii=radii)
check("NN1 stable across chunk sizes (<2e-3)", np.max(np.abs(nn1 - nn1b)) < 2e-3,
      f"maxdiff={np.max(np.abs(nn1-nn1b)):.1e}")
check("radius-count means stable across chunk sizes (<0.05)",
      np.abs(rad.mean(0) - radb.mean(0)).max() < 0.05,
      f"maxΔmean={np.abs(rad.mean(0)-radb.mean(0)).max():.3f}")

# ---------------------------------------------------------------- PCA vs numpy float64
print("\n# 3. PCA spectrum vs numpy float64 eig")
Xc = E.astype(np.float64) - E.astype(np.float64).mean(0)
cov = Xc.T @ Xc / (len(E) - 1)
ev = np.linalg.eigvalsh(cov)[::-1].clip(min=0)
pr_ref = ev.sum()**2 / (ev**2).sum()
cum = np.cumsum(ev / ev.sum())
n90_ref = int(np.searchsorted(cum, 0.90) + 1)
# analyze_space computes the same on GPU float32 inside main(); replicate its path
Ef = Eg.float(); Xc2 = Ef - Ef.mean(0, keepdim=True)
covg = (Xc2.T @ Xc2) / (len(E) - 1)
evg = torch.linalg.eigvalsh(covg).flip(0).clamp(min=0).cpu().numpy()
pr_gpu = evg.sum()**2 / (evg**2).sum()
n90_gpu = int(np.searchsorted(np.cumsum(evg/evg.sum()), 0.90) + 1)
check("participation ratio matches numpy", abs(pr_gpu - pr_ref) / pr_ref < 0.02, f"gpu={pr_gpu:.2f} ref={pr_ref:.2f}")
check("n_pcs_90% matches numpy", abs(n90_gpu - n90_ref) <= 1, f"gpu={n90_gpu} ref={n90_ref}")

# ---------------------------------------------------------------- TwoNN ground truth
print("\n# 4. TwoNN recovers KNOWN intrinsic dimension")
def sphere(n, k, D=1024, seed=1):
    g = np.random.default_rng(seed)
    x = np.zeros((n, D), np.float32)
    x[:, :k] = g.standard_normal((n, k))           # uniform on S^{k-1} after norm -> intrinsic dim k-1
    x /= np.linalg.norm(x, axis=1, keepdims=True)
    return torch.from_numpy(x.astype(np.float16)).cuda()
ests = {}
for k_dim in (4, 6, 11, 21):
    est = A.twonn(sphere(12000, k_dim), sample=12000)
    expected = k_dim - 1
    ests[expected] = est
    ok = abs(est - expected) / expected < 0.20      # approximate estimator; ±20%
    check(f"S^{k_dim-1} (intrinsic dim {expected})", ok, f"-> TwoNN={est:.2f}")
check("TwoNN monotonic in true dim", all(x < y for x, y in zip(
      [ests[d] for d in sorted(ests)], [ests[d] for d in sorted(ests)][1:])))

print("\nRESULT:", "ALL PASS" if not fails else f"{len(fails)} FAILED: {fails}")
sys.exit(1 if fails else 0)
