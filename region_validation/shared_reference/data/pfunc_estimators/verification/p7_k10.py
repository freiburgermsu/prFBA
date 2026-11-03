"""Why the 200-hit K=10 q>=0.5 macro F1 differs from calibrate_pfunc by 1.3e-6: summation
order at q = 0.5 (dot product vs cumulative sum)."""
import numpy as np
import p2_pass as P
P.setup()
G = P.G
E = G["E"]
f1_dot, f1_cum, n_near = [], [], 0
for i in range(len(G["cells"])):
    f = int(G["fold"][i])
    g2, i2 = G["l200"][i]
    kk = min(10, len(g2))
    X = G["Ef"][f][g2[:100]].astype(np.float64)
    w = np.exp(-(i2[0] - i2[:100]) / 0.01)
    qd = (w[:kk] @ X[:kk]) / w[:kk].sum()
    cum = np.cumsum(X * w[:, None], axis=0)
    qc = cum[kk - 1] / np.cumsum(w)[kk - 1]
    y = E[G["src"][i]].astype(bool)
    T = y.sum()
    n_near += int((np.abs(qd - 0.5) < 1e-12).sum())
    for q, lst in ((qd, f1_dot), (qc, f1_cum)):
        q = np.where(q >= 1 - 1e-9, 1.0, q)
        p = q >= 0.5
        lst.append(2 * (p & y).sum() / (p.sum() + T))
print("macro F1 dot-product:", np.mean(f1_dot), " cumulative-sum:", np.mean(f1_cum),
      " EC points with |q-0.5|<1e-12:", n_near)
