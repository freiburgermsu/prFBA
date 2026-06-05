#!/usr/bin/env python
"""
analyze_space.py — characterize the density / sparsity / geometry of the 16S embedding space.

All vectors are unit-norm, so they live on the surface of the 1023-sphere and cosine
similarity == dot product; cosine *distance* = 1 - cos. We quantify how the cloud fills
that sphere: where it is dense (tight clusters of near-identical / con-generic 16S) vs
sparse (isolated, novel lineages), and the manifold's effective dimensionality.

Computed (GPU, chunked exact kNN — no approximation):
  * NN1 / k-NN cosine per point        -> packing tightness
  * neighbors within cosine radii      -> local density counts
  * sampled pairwise cosine            -> global background concentration
  * PCA spectrum (full covariance)     -> effective / participation-ratio dimensionality
  * TwoNN intrinsic dimension (sample) -> manifold dimensionality
  * k-NN taxonomic purity (genus)      -> does proximity == taxonomy
  * isolation fraction (NN1 < thr)     -> size of the sparse frontier
  * UMAP 2-D (subsample)               -> visualization, colored by genus + local density

Outputs to the store dir: density_stats.json + PNG figures + umap_coords.npy.
"""
from __future__ import annotations
import argparse, json, time
from pathlib import Path
import numpy as np

STORE = Path(__file__).resolve().parent


def log(m): print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def gpu_knn(Eg, k, chunk, radii):
    """Exact cosine kNN on GPU. Returns nn1, nnk(mean of top-k), and counts within each radius."""
    import torch
    N = Eg.shape[0]
    nn1 = torch.empty(N, device=Eg.device)
    nnk = torch.empty(N, device=Eg.device)
    nn1_idx = torch.empty(N, dtype=torch.long, device=Eg.device)
    rad = torch.zeros((N, len(radii)), device=Eg.device)
    radii_t = torch.tensor(radii, device=Eg.device)
    kk = min(k + 1, N)
    for s in range(0, N, chunk):
        q = Eg[s:s + chunk]
        sims = (q @ Eg.T).float()                      # (b, N) cosine
        b = sims.shape[0]
        sims[torch.arange(b, device=Eg.device), torch.arange(s, s + b, device=Eg.device)] = -1.0  # mask self
        top = torch.topk(sims, k=kk - 1, dim=1)
        nn1[s:s + b] = top.values[:, 0]
        nn1_idx[s:s + b] = top.indices[:, 0]
        nnk[s:s + b] = top.values[:, :k].mean(dim=1)
        for j, r in enumerate(radii):                  # per-radius to bound memory
            rad[s:s + b, j] = (sims >= r).sum(dim=1).float()
    return (nn1.cpu().numpy(), nnk.cpu().numpy(), nn1_idx.cpu().numpy(),
            rad.cpu().numpy())


def twonn(Eg, sample, seed=0, trim=0.05, chunk=2048):
    """
    TwoNN intrinsic-dimension estimator (Facco et al. 2017) on chord (Euclidean-on-sphere)
    distance, which is what the estimator's locally-Euclidean derivation assumes. Uses the
    maximum-likelihood form d = n / sum(ln mu), mu = r2/r1, after discarding the top `trim`
    fraction of mu (the heavy upper tail from near-duplicate / outlier pairs). Validated to
    recover known intrinsic dimensions of synthetic S^d manifolds to ~10% (mean rel. error
    over d in {3,5,10,20}); trim=0.05 was selected as the least-biased setting.
    """
    import torch
    g = torch.Generator(device="cpu").manual_seed(seed)
    sel = torch.randperm(Eg.shape[0], generator=g)[:sample].to(Eg.device)
    tops = []
    for s in range(0, len(sel), chunk):                # chunk queries to bound memory
        si = sel[s:s + chunk]
        sims = (Eg[si] @ Eg.T).float().clamp(-1, 1)
        sims[torch.arange(len(si), device=Eg.device), si] = 2.0   # exclude self
        tops.append(torch.topk(sims, k=3, dim=1).values)
    top = torch.cat(tops, 0)                            # [self(2.0), nn1_cos, nn2_cos]
    r1 = torch.sqrt((2 - 2 * top[:, 1]).clamp(min=0))   # chord distance to 1st NN
    r2 = torch.sqrt((2 - 2 * top[:, 2]).clamp(min=0))   # to 2nd NN
    mu = (r2 / r1.clamp(min=1e-9)).cpu().numpy()
    mu = mu[np.isfinite(mu) & (mu > 1.0 + 1e-9)]
    mu.sort()
    cut = int(len(mu) * (1 - trim))                    # drop heavy upper tail (near-dup/outlier pairs)
    mu = mu[:cut]
    return float(len(mu) / np.log(mu).sum())           # MLE of the Pareto(d) ratio law


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--store", type=Path, default=STORE)
    ap.add_argument("--k", type=int, default=10, help="k for k-NN density / purity")
    ap.add_argument("--chunk", type=int, default=2048)
    ap.add_argument("--radii", type=float, nargs="+", default=[0.99, 0.97, 0.95, 0.90, 0.80])
    ap.add_argument("--iso-thr", type=float, default=0.80, help="NN1 below this = isolated/sparse")
    ap.add_argument("--pair-samples", type=int, default=5_000_000)
    ap.add_argument("--twonn-sample", type=int, default=20000)
    ap.add_argument("--umap-sample", type=int, default=60000)
    ap.add_argument("--min-len", type=int, default=None, help="restrict analyzed set to seq_len >= this")
    ap.add_argument("--max-len", type=int, default=None, help="restrict analyzed set to seq_len <= this")
    ap.add_argument("--exclude-atypical", action="store_true", help="drop atypical-length rows from the analysis")
    ap.add_argument("--tag", default="", help="suffix for output files so filtered runs don't clobber the full run")
    ap.add_argument("--device", default="cuda")
    args = ap.parse_args()

    import torch, pandas as pd
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    t0 = time.time()
    idx = pd.read_parquet(args.store / "index.parquet")
    E = np.load(args.store / "embeddings.f16.npy", mmap_mode="r")
    N0 = len(idx)
    # optional length / atypical filter on the analyzed reference set
    sl_all = idx["seq_len"].values
    keep = np.ones(N0, bool)
    if args.min_len: keep &= sl_all >= args.min_len
    if args.max_len: keep &= sl_all <= args.max_len
    if args.exclude_atypical: keep &= ~idx["atypical_len"].values.astype(bool)
    E = np.ascontiguousarray(np.asarray(E)[keep] if not keep.all() else np.asarray(E))
    idx = idx[keep].reset_index(drop=True)
    N, D = E.shape
    tag = f"_{args.tag}" if args.tag else ""
    def P(name):                                                       # tag-aware output path
        stem, dot, ext = name.rpartition(".")
        return args.store / (f"{stem}{tag}.{ext}" if dot else f"{name}{tag}")
    log(f"loaded {N:,} x {D} embeddings" + (f"  (filtered from {N0:,})" if N != N0 else ""))
    device = args.device if torch.cuda.is_available() else "cpu"
    Eg = torch.from_numpy(E).to(device)                               # f16 on GPU
    stats = {"n_points": int(N), "n_points_unfiltered": int(N0), "dim": int(D), "k": args.k,
             "filter": {"min_len": args.min_len, "max_len": args.max_len,
                        "exclude_atypical": args.exclude_atypical}}
    # auto-shrink the query chunk so a (chunk x N) float32 score block stays <= ~2.0 GB
    eff_chunk = max(256, min(args.chunk, int(2.0e9 / (N * 4))))
    if eff_chunk != args.chunk:
        log(f"auto chunk: {eff_chunk} (N={N:,})")

    # ---- kNN density ----------------------------------------------------------
    log("exact GPU kNN (NN1, mean-kNN, radius counts) ...")
    nn1, nnk, nn1_idx, rad = gpu_knn(Eg, args.k, eff_chunk, args.radii)
    torch.cuda.empty_cache()
    stats["nn1_cosine"] = pct(nn1)
    stats["mean_kNN_cosine"] = pct(nnk)
    stats["isolated_fraction(NN1<%.2f)" % args.iso_thr] = float((nn1 < args.iso_thr).mean())
    stats["neighbors_within_radius_mean"] = {f"cos>={r}": float(rad[:, i].mean()) for i, r in enumerate(args.radii)}
    stats["singletons_no_neighbor_above_0.95"] = float((rad[:, args.radii.index(0.95)] == 0).mean()) \
        if 0.95 in args.radii else None

    # ---- global pairwise sample ----------------------------------------------
    log("sampling global pairwise cosine ...")
    g = torch.Generator(device="cpu").manual_seed(0)
    a = torch.randint(0, N, (args.pair_samples,), generator=g)
    b = torch.randint(0, N, (args.pair_samples,), generator=g)
    keep = a != b; a, b = a[keep], b[keep]
    parts = []
    for s in range(0, len(a), 1_000_000):              # chunk gather to bound GPU memory
        ai = a[s:s + 1_000_000].to(device); bi = b[s:s + 1_000_000].to(device)
        parts.append((Eg[ai] * Eg[bi]).sum(1).float().cpu().numpy())
    pa = np.concatenate(parts)
    stats["global_pairwise_cosine"] = pct(pa)
    torch.cuda.empty_cache()

    # ---- PCA spectrum (effective dimensionality) ------------------------------
    log("PCA spectrum on full covariance ...")
    Ef = Eg.float()
    mean = Ef.mean(0, keepdim=True)
    Xc = Ef - mean
    cov = (Xc.T @ Xc) / (N - 1)                                        # 1024 x 1024
    del Ef, Xc; torch.cuda.empty_cache()
    evals = torch.linalg.eigvalsh(cov).flip(0).clamp(min=0).cpu().numpy()
    ev_ratio = evals / evals.sum()
    cum = np.cumsum(ev_ratio)
    stats["pca"] = {
        "var_explained_top1": float(ev_ratio[0]),
        "var_explained_top10": float(cum[9]),
        "n_pcs_50pct": int(np.searchsorted(cum, 0.50) + 1),
        "n_pcs_90pct": int(np.searchsorted(cum, 0.90) + 1),
        "n_pcs_95pct": int(np.searchsorted(cum, 0.95) + 1),
        "n_pcs_99pct": int(np.searchsorted(cum, 0.99) + 1),
        "participation_ratio": float((evals.sum() ** 2) / (np.square(evals).sum())),
    }

    # ---- TwoNN intrinsic dimension -------------------------------------------
    log("TwoNN intrinsic dimension ...")
    stats["intrinsic_dim_twoNN"] = twonn(Eg, min(args.twonn_sample, N), chunk=eff_chunk)
    torch.cuda.empty_cache()

    # ---- taxonomic purity of neighborhoods -----------------------------------
    log("k-NN taxonomic purity (genus) ...")
    def _genus(org):
        if not isinstance(org, str) or not org.strip():
            return ""
        toks = org.replace("[", "").replace("]", "").split()
        if not toks:
            return ""
        g = toks[0]
        if g.lower() in ("candidatus", "uncultured", "unclassified") and len(toks) > 1:
            g = toks[1]
        return g
    genus = np.array([_genus(o) for o in idx["organism"].values])    # genus label per row
    has_g = genus != ""
    # recompute k neighbors' indices for purity on a sample to bound memory
    sel = np.random.default_rng(0).choice(np.where(has_g)[0], size=min(20000, has_g.sum()), replace=False)
    Qs = Eg[torch.tensor(sel, device=device)]
    purity = []
    for s in range(0, len(sel), eff_chunk):
        qs = Qs[s:s + eff_chunk]
        sims = (qs @ Eg.T).float()
        for r, gi in enumerate(sel[s:s + eff_chunk]):
            sims[r, gi] = -1
        nbr = torch.topk(sims, k=args.k, dim=1).indices.cpu().numpy()
        for r, gi in enumerate(sel[s:s + eff_chunk]):
            ng = genus[nbr[r]]
            purity.append(np.mean(ng == genus[gi]))
    purity = np.array(purity)
    stats["knn_genus_purity"] = pct(purity)

    # ---- multiplicity / length context ---------------------------------------
    sl = idx["seq_len"].values.astype(float)
    stats["seq_len"] = pct(sl)
    stats["multiplicity"] = pct(idx["multiplicity"].values.astype(float))
    stats["pct_atypical_len"] = float(idx["atypical_len"].mean())

    # ---- length-stratified metrics -------------------------------------------
    log("length-stratified metrics ...")
    r95 = rad[:, args.radii.index(0.95)] if 0.95 in args.radii else rad[:, 0]
    BANDS = [("<800", 0, 800), ("800-1200", 800, 1200), ("1200-1400", 1200, 1400),
             ("1400-1600", 1400, 1600), (">=1600", 1600, 10**18)]
    sl_sel = sl[sel]                                                   # lengths of the purity-sample points
    strata = {}
    for name, lo, hi in BANDS:
        m = (sl >= lo) & (sl < hi)
        if not m.any():
            continue
        d = {"n": int(m.sum()), "frac": float(m.mean()),
             "nn1_cos_mean": float(nn1[m].mean()), "nn1_cos_median": float(np.median(nn1[m])),
             "mean_kNN_cos": float(nnk[m].mean()),
             "isolated_frac": float((nn1[m] < args.iso_thr).mean()),
             "neighbors_cos95_mean": float(r95[m].mean())}
        ms = (sl_sel >= lo) & (sl_sel < hi)
        if ms.sum() > 30:
            d["genus_purity_mean"] = float(purity[ms].mean())
        strata[name] = d
    stats["length_strata"] = strata

    # cross-band: do partial fragments (<1400 bp) land near the full-length core (1400-1600 bp)?
    full_m = (sl >= 1400) & (sl < 1600)
    part_m = sl < 1400
    best = None
    if full_m.any() and part_m.any():
        Ef_full = Eg[torch.tensor(np.where(full_m)[0], device=device)]
        pidx = np.where(part_m)[0]
        chunks = []
        for s in range(0, len(pidx), eff_chunk):
            q = Eg[torch.tensor(pidx[s:s + eff_chunk], device=device)]
            chunks.append((q @ Ef_full.T).float().max(1).values.cpu().numpy())
        best = np.concatenate(chunks)
        stats["partial_to_fulllength"] = {
            "n_partials": int(part_m.sum()), "n_fulllength": int(full_m.sum()),
            "best_cos_to_fulllength": pct(best),
            "frac_partials_match_fulllength_ge_0.95": float((best >= 0.95).mean()),
            "frac_partials_orphan_lt_0.80": float((best < 0.80).mean())}
        del Ef_full; torch.cuda.empty_cache()
    stats["wall_seconds"] = round(time.time() - t0, 1)

    json.dump(stats, open(P("density_stats.json"), "w"), indent=2,
              default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else o)
    log("wrote " + P("density_stats.json").name)

    # ---- figures --------------------------------------------------------------
    log("figures ...")
    _hist(nn1, "cosine to nearest neighbor (NN1)", "packing tightness — right=dense, left=sparse/novel",
          P("fig_nn1_hist.png"), vlines=[args.iso_thr])
    _hist(r95, "# neighbors with cos>=0.95", "local density (log y)",
          P("fig_local_density.png"), logy=True)
    _hist(pa, "cosine of random pairs", "global background similarity", P("fig_pairwise_hist.png"))
    _pca_plot(ev_ratio, cum, stats["pca"], P("fig_pca_spectrum.png"))
    _strata_plot(sl, nn1, BANDS, args.iso_thr, P("fig_length_strata.png"))
    if best is not None:
        _partial_plot(best, P("fig_partial_vs_fulllength.png"))

    # ---- UMAP -----------------------------------------------------------------
    try:
        import umap
        ns = min(args.umap_sample, N)
        rs = np.random.default_rng(0).choice(N, size=ns, replace=False)
        sub = np.ascontiguousarray(E[np.sort(rs)]).astype(np.float32)
        rsi = np.sort(rs)
        log(f"UMAP on {ns:,} sampled points ...")
        reducer = umap.UMAP(n_neighbors=30, min_dist=0.1, metric="cosine", random_state=42, verbose=False)
        xy = reducer.fit_transform(sub)
        np.save(P("umap_coords.npy"), xy)
        np.save(P("umap_row_ids.npy"), rsi)
        dens = rad[rsi, args.radii.index(0.95)] if 0.95 in args.radii else rad[rsi, 0]
        _umap_plots(xy, genus[rsi], np.log1p(dens), P("fig_umap_genus.png"), P("fig_umap_density.png"))
    except Exception as e:
        log(f"UMAP skipped: {e}")
    log(f"analysis DONE in {(time.time()-t0)/60:.1f} min")
    print(json.dumps(stats, indent=2, default=str))


def pct(a):
    a = np.asarray(a, float)
    return {"mean": float(a.mean()), "std": float(a.std()),
            "p1": float(np.percentile(a, 1)), "p5": float(np.percentile(a, 5)),
            "p25": float(np.percentile(a, 25)), "median": float(np.percentile(a, 50)),
            "p75": float(np.percentile(a, 75)), "p95": float(np.percentile(a, 95)),
            "p99": float(np.percentile(a, 99)), "min": float(a.min()), "max": float(a.max())}


def _hist(x, xlabel, title, path, logy=False, vlines=()):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.hist(x, bins=120, color="#2b6cb0", alpha=0.85)
    if logy: ax.set_yscale("log")
    for v in vlines: ax.axvline(v, color="crimson", ls="--", lw=1)
    ax.set_xlabel(xlabel); ax.set_ylabel("count"); ax.set_title(title)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def _pca_plot(ev_ratio, cum, p, path):
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    ax[0].plot(np.arange(1, 51), ev_ratio[:50], "o-", ms=3, color="#2b6cb0")
    ax[0].set_xlabel("principal component"); ax[0].set_ylabel("variance ratio"); ax[0].set_title("PCA scree (top 50)")
    ax[1].plot(np.arange(1, len(cum) + 1), cum, color="#2b6cb0")
    for q, c in [(0.5, "#888"), (0.9, "crimson"), (0.95, "darkorange")]:
        ax[1].axhline(q, ls="--", lw=0.8, color=c)
    ax[1].set_xlabel("# components"); ax[1].set_ylabel("cumulative variance")
    ax[1].set_title(f"PR≈{p['participation_ratio']:.0f}; 90% at {p['n_pcs_90pct']} PCs")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def _umap_plots(xy, genus, logdens, path_genus, path_density):
    import matplotlib.pyplot as plt
    # by top genera
    import collections
    counts = collections.Counter(g for g in genus if g)
    topg = [g for g, _ in counts.most_common(12)]
    fig, ax = plt.subplots(figsize=(8, 7))
    other = ~np.isin(genus, topg)
    ax.scatter(xy[other, 0], xy[other, 1], s=1, c="#dddddd", alpha=0.4, label="other")
    cmap = plt.cm.tab20(np.linspace(0, 1, len(topg)))
    for g, col in zip(topg, cmap):
        m = genus == g
        ax.scatter(xy[m, 0], xy[m, 1], s=2, color=col, label=f"{g} ({counts[g]})", alpha=0.7)
    ax.legend(markerscale=4, fontsize=7, loc="best", ncol=2)
    ax.set_title("16S embedding (UMAP, cosine) — colored by top genera"); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(path_genus, dpi=150); plt.close(fig)
    # by local density
    fig, ax = plt.subplots(figsize=(8, 7))
    sc = ax.scatter(xy[:, 0], xy[:, 1], s=2, c=logdens, cmap="viridis", alpha=0.7)
    fig.colorbar(sc, label="log(1 + #neighbors cos>=0.95)")
    ax.set_title("16S embedding (UMAP) — local density"); ax.set_xticks([]); ax.set_yticks([])
    fig.tight_layout(); fig.savefig(path_density, dpi=150); plt.close(fig)


def _strata_plot(sl, nn1, bands, iso_thr, path):
    """NN1 distribution and isolated-fraction per length band."""
    import matplotlib.pyplot as plt
    present = [(name, lo, hi) for name, lo, hi in bands if ((sl >= lo) & (sl < hi)).any()]
    data = [nn1[(sl >= lo) & (sl < hi)] for _, lo, hi in present]
    labels = [f"{name}\n(n={len(d)})" for (name, _, _), d in zip(present, data)]
    iso = [float((d < iso_thr).mean()) for d in data]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    ax[0].boxplot(data, tick_labels=labels, showfliers=False)
    ax[0].set_ylabel("NN1 cosine"); ax[0].set_title("packing tightness by length band")
    ax[0].tick_params(axis="x", labelsize=8)
    ax[1].bar(range(len(present)), iso, color="#2b6cb0")
    ax[1].set_xticks(range(len(present))); ax[1].set_xticklabels([p[0] for p in present], rotation=20)
    ax[1].set_ylabel(f"isolated fraction (NN1<{iso_thr})"); ax[1].set_title("sparsity by length band")
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


def _partial_plot(best, path):
    """Best cosine of each partial fragment to the full-length core."""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6.5, 4))
    ax.hist(best, bins=120, color="#2b6cb0", alpha=0.85)
    ax.axvline(0.95, color="crimson", ls="--", lw=1, label="0.95 (species)")
    ax.axvline(0.80, color="darkorange", ls="--", lw=1, label="0.80 (orphan)")
    ax.set_xlabel("best cosine of a partial (<1400 bp) to any full-length (1400–1600 bp) 16S")
    ax.set_ylabel("count"); ax.set_title("do partial fragments map onto the full-length core?")
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=140); plt.close(fig)


if __name__ == "__main__":
    import pandas as pd
    main()
