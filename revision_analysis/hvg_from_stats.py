"""streaming の平均・分散から scanpy の HVG(flavor='seurat') を再現する。

scanpy._highly_variable_genes_single_batch(flavor='seurat') と同じ手順:
  CP10K の mean/var(ddof=1) → dispersion=var/mean → log → mean=log1p(mean)
  → mean を等幅20ビンに切り、ビン内で dispersion を z 化 → 上位 n_top_genes
"""
import numpy as np, pandas as pd

def hvg_from_stats(mean_cp, var_cp, n_top_genes, n_bins=20):
    mean = np.asarray(mean_cp, dtype=float).copy()
    var  = np.asarray(var_cp, dtype=float).copy()
    mean[mean == 0] = 1e-12
    disp = var / mean
    disp[disp == 0] = np.nan
    disp = np.log(disp)
    mean = np.log1p(mean)
    df = pd.DataFrame({"means": mean, "dispersions": disp})
    df["mean_bin"] = pd.cut(df.means, bins=n_bins)
    g = df.groupby("mean_bin", observed=True)["dispersions"]
    dm, ds = g.mean(), g.std(ddof=1)
    one = ds.isnull()
    ds[one] = dm[one]; dm[one] = 0
    df["dispersions_norm"] = (df.dispersions.values - dm.loc[df.mean_bin].values) / ds.loc[df.mean_bin].values
    order = df.dispersions_norm.fillna(-np.inf).sort_values(ascending=False).index
    hv = np.zeros(len(df), dtype=bool); hv[order[:n_top_genes]] = True
    return hv, df
