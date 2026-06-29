"""Stage 4c fallback: compute pseudotime peak timing for progression-coupled lncRNAs
WITHOUT BRIC-seq half-life data. All cells get stability_class = 'unknown'.

This lets Fig 3 (waves) and Fig 5 (stability) render with partial information.
When the Tani 2012 half_lives.csv becomes available, run 04c_stability.py for the
full version.
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc

from pipeline import config


def find_peak_pseudotime(a, gene, lineage, n_bins=10):
    if gene not in a.var_names:
        return None
    g_idx = a.var_names.get_loc(gene)
    mask = (a.obs["lineage"] == lineage) & ~a.obs["pseudotime"].isna()
    if mask.sum() < 100:
        return None
    pt = a.obs.loc[mask, "pseudotime"].values
    X = a.X[:, g_idx].toarray().ravel() if hasattr(a.X, "toarray") else np.asarray(a.X[:, g_idx]).ravel()
    expr = X[mask.values]
    bins = pd.qcut(pt, n_bins, labels=False, duplicates="drop")
    means = pd.Series(expr).groupby(bins).mean()
    if means.empty:
        return None
    return int(means.idxmax())


def main():
    a = sc.read_h5ad(config.RESULTS_DIR / "liver_anchor_pseudotime.h5ad")
    progression = pd.read_csv(config.RESULTS_DIR / "progression_lncrna_table.csv")
    progression = progression[progression["is_progression_coupled"]]

    rows = []
    for _, r in progression.iterrows():
        gene = r["gene"]
        peak = find_peak_pseudotime(a, gene, r["lineage"], n_bins=config.SCREENING["n_bins"])
        rows.append({
            "gene_id": gene,
            "gene_name": gene,  # var_names are already symbols in this build
            "lineage": r["lineage"],
            "stability_class": "unknown",
            "peak_bin": peak,
            "spearman_rho": r["spearman_rho"],
            "log2fc_max_min": r["log2fc_max_min"],
        })

    out = pd.DataFrame(rows)
    out.to_csv(config.RESULTS_DIR / "stability_x_pseudotime.csv", index=False)
    print("Saved (all stability_class='unknown' placeholder):")
    print(out.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
