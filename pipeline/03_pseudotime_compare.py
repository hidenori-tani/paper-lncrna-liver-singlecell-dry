"""Stage 3 prelim: compare Monocle3 vs PAGA+DPT pseudotime on Hepatocyte subset.

Output: results/pseudotime_compare.csv with columns: tool, n_cells, smoothness, donor_consistency
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc

from pipeline import config


def subset_hepatocyte(a):
    return a[a.obs["lineage"] == "Hepatocyte"].copy()


def run_paga_dpt(a):
    sc.tl.paga(a, groups="leiden")
    healthy_cells = a.obs.index[a.obs.get("condition", "healthy") == "healthy"]
    if len(healthy_cells) == 0:
        healthy_cells = a.obs.index[:100]
    a.uns["iroot"] = a.obs_names.get_loc(np.random.choice(healthy_cells))
    sc.tl.dpt(a)
    return a.obs["dpt_pseudotime"].values


def run_monocle3_via_r(a):
    """Calls Monocle3 in R via rpy2. Fallback: skip and return NaN if rpy2/R unavailable."""
    try:
        import rpy2.robjects as ro
        from rpy2.robjects import pandas2ri
        pandas2ri.activate()
        from anndata2ri import py2rpy
        sce = py2rpy(a)
        ro.r("library(monocle3)")
        ro.globalenv["sce"] = sce
        ro.r("cds <- as.cell_data_set(sce)")
        ro.r("cds <- preprocess_cds(cds, num_dim = 50)")
        ro.r("cds <- reduce_dimension(cds, reduction_method='UMAP')")
        ro.r("cds <- cluster_cells(cds)")
        ro.r("cds <- learn_graph(cds, use_partition=FALSE)")
        ro.r("cds <- order_cells(cds, root_cells=colnames(cds)[1])")
        pt = ro.r("pseudotime(cds)")
        return np.array(pt)
    except Exception as e:
        print(f"  Monocle3 unavailable ({e}); returning NaN")
        return np.full(a.n_obs, np.nan)


def smoothness_score(pt):
    """Lower = smoother. RMS of adjacent pseudotime differences after sorting."""
    s = np.sort(pt[~np.isnan(pt)])
    if len(s) < 2:
        return np.nan
    return float(np.sqrt(np.mean(np.diff(s) ** 2)))


def donor_consistency(a, pt, donor_key):
    """Spearman rho between donor medians and global rank — higher = more consistent."""
    from scipy.stats import spearmanr
    df = pd.DataFrame({"pt": pt, "donor": a.obs[donor_key].values})
    medians = df.groupby("donor")["pt"].median().dropna()
    if len(medians) < 3:
        return np.nan
    ranks = medians.rank()
    rho, _ = spearmanr(medians, ranks)
    return float(rho)


def main():
    a = sc.read_h5ad(config.RESULTS_DIR / "liver_anchor_atlas.h5ad")
    h = subset_hepatocyte(a)
    print(f"Hepatocyte subset: {h.shape}")

    donor_key = "donor_id" if "donor_id" in h.obs.columns else "sample"

    pt_paga = run_paga_dpt(h.copy())
    pt_mon = run_monocle3_via_r(h.copy())

    rows = []
    for name, pt in [("paga_dpt", pt_paga), ("monocle3", pt_mon)]:
        rows.append({
            "tool": name,
            "n_cells": int((~np.isnan(pt)).sum()),
            "smoothness_lower_better": smoothness_score(pt),
            "donor_consistency_higher_better": donor_consistency(h, pt, donor_key),
        })
    out = pd.DataFrame(rows)
    out.to_csv(config.RESULTS_DIR / "pseudotime_compare.csv", index=False)
    print(out.to_string(index=False))
    print("\nDecision criterion: choose the tool with higher donor_consistency AND smoothness within 2x of other.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
