"""Stage 4b: Focused deep-dive on anchor lncRNA panel.

For each anchor lncRNA:
  - per-lineage mean expression vs pseudotime bin
  - cell-type composition where lncRNA is detected (count > 0)
  - top 20 co-expressed genes per lineage
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc

from pipeline import config, lncrna_utils


def per_lineage_trajectory(a, gene, lineages, n_bins=10):
    rows = []
    if gene not in a.var_names:
        return pd.DataFrame()
    g_idx = a.var_names.get_loc(gene)
    X = a.X[:, g_idx].toarray().ravel() if hasattr(a.X, "toarray") else np.asarray(a.X[:, g_idx]).ravel()
    for lin in lineages:
        mask = (a.obs["lineage"] == lin) & ~a.obs["pseudotime"].isna()
        n = int(mask.sum())
        if n < 100:
            continue
        pt = a.obs.loc[mask, "pseudotime"].values
        expr = X[mask.values]
        bins = pd.qcut(pt, n_bins, labels=False, duplicates="drop")
        means = pd.Series(expr).groupby(bins).mean()
        for b, m in means.items():
            rows.append({"gene": gene, "lineage": lin, "bin": int(b), "mean_expr": float(m), "n_cells": n})
    return pd.DataFrame(rows)


def composition(a, gene):
    if gene not in a.var_names:
        return pd.DataFrame()
    g_idx = a.var_names.get_loc(gene)
    X = a.X[:, g_idx].toarray().ravel() if hasattr(a.X, "toarray") else np.asarray(a.X[:, g_idx]).ravel()
    detected = X > 0
    df = pd.DataFrame({"lineage": a.obs["lineage"].values, "detected": detected})
    return df.groupby("lineage")["detected"].agg(["sum", "count"]).rename(columns={"sum": "n_detected", "count": "n_total"}).reset_index().assign(gene=gene)


def coexpression_top(a, gene, lineage, top_n=20):
    if gene not in a.var_names:
        return pd.DataFrame()
    sub = a[a.obs["lineage"] == lineage]
    if sub.n_obs < 200:
        return pd.DataFrame()
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    g_idx = sub.var_names.get_loc(gene)
    g_vec = X[:, g_idx]
    if np.std(g_vec) == 0:
        return pd.DataFrame()
    from scipy.stats import spearmanr
    rhos = np.zeros(sub.n_vars)
    for i in range(sub.n_vars):
        if i == g_idx or np.std(X[:, i]) == 0:
            continue
        rhos[i], _ = spearmanr(g_vec, X[:, i])
    df = pd.DataFrame({"gene": gene, "lineage": lineage, "partner": sub.var_names, "rho": rhos})
    return df.assign(abs_rho=df["rho"].abs()).nlargest(top_n, "abs_rho").drop(columns=["abs_rho"])


def main():
    a = sc.read_h5ad(config.RESULTS_DIR / "liver_anchor_pseudotime.h5ad")
    name_to_id = lncrna_utils.get_lncrna_name_to_id()

    traj_rows = []
    comp_rows = []
    coex_rows = []
    for name in config.ANCHOR_LNCRNAS:
        target = name if name in a.var_names else name_to_id.get(name)
        if target is None or target not in a.var_names:
            print(f"  {name}: NOT IN ATLAS, skipping")
            continue
        traj_rows.append(per_lineage_trajectory(a, target, config.PSEUDOTIME["lineages"]))
        comp_rows.append(composition(a, target))
        for lin in config.PSEUDOTIME["lineages"]:
            coex_rows.append(coexpression_top(a, target, lin))
        print(f"  {name} ({target}): done")

    if traj_rows:
        pd.concat(traj_rows, ignore_index=True).to_csv(config.RESULTS_DIR / "anchor_trajectory.csv", index=False)
    if comp_rows:
        pd.concat(comp_rows, ignore_index=True).to_csv(config.RESULTS_DIR / "anchor_composition.csv", index=False)
    if coex_rows:
        pd.concat(coex_rows, ignore_index=True).to_csv(config.RESULTS_DIR / "anchor_coexpression.csv", index=False)
    print("Saved anchor_trajectory.csv, anchor_composition.csv, anchor_coexpression.csv")
    return 0


if __name__ == "__main__":
    sys.exit(main())
