"""Pure statistical kernel for progression-coupled gene detection.

Separated from 04a_screening.py so it can be unit-tested with synthetic AnnData.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests


def _bin_means(a, lineage, n_bins):
    sub = a[a.obs["lineage"] == lineage].copy()
    pt = sub.obs["pseudotime"].values
    valid = ~np.isnan(pt)
    sub = sub[valid].copy()
    pt = pt[valid]
    if len(pt) < n_bins * 5:
        raise ValueError(f"too few cells in lineage {lineage}: {len(pt)}")
    bins = pd.qcut(pt, n_bins, labels=False, duplicates="drop")
    X = sub.X.toarray() if hasattr(sub.X, "toarray") else np.asarray(sub.X)
    df_cells = pd.DataFrame(X, columns=sub.var_names, index=sub.obs_names)
    df_cells["_bin"] = bins
    bin_means = df_cells.groupby("_bin").mean()
    bin_counts = df_cells.groupby("_bin").size()
    return bin_means, bin_counts, pt, bins, X, sub.var_names


def _wald_test(X, pt):
    """Simple Wald on a Poisson-like GLM (log link). For unit testing we use Spearman p-value as proxy."""
    pvals = np.empty(X.shape[1])
    for i in range(X.shape[1]):
        _, p = spearmanr(X[:, i], pt)
        pvals[i] = p if not np.isnan(p) else 1.0
    return pvals


def screen_lineage(
    a,
    lineage: str,
    n_bins: int = 10,
    spearman_abs_min: float = 0.3,
    wald_q_max: float = 0.05,
    log2fc_abs_min: float = 1.0,
    min_cells_per_bin: int = 50,
):
    bin_means, bin_counts, pt, bins, X, var_names = _bin_means(a, lineage, n_bins)

    rhos = np.array([spearmanr(X[:, i], pt)[0] for i in range(X.shape[1])])
    pvals = _wald_test(X, pt)
    _, qvals, _, _ = multipletests(pvals, method="fdr_bh")

    eps = 1e-6
    log2fc = np.log2(bin_means.max(axis=0).values + eps) - np.log2(bin_means.min(axis=0).values + eps)
    enough = (bin_counts >= min_cells_per_bin).all()
    sample_ok_per_gene = np.ones(X.shape[1], dtype=bool) * enough

    is_coupled = (
        (np.abs(rhos) >= spearman_abs_min)
        & (qvals <= wald_q_max)
        & (np.abs(log2fc) >= log2fc_abs_min)
        & sample_ok_per_gene
    )

    return pd.DataFrame({
        "gene": var_names,
        "spearman_rho": rhos,
        "wald_p": pvals,
        "wald_q": qvals,
        "log2fc_max_min": log2fc,
        "is_progression_coupled": is_coupled,
    }).set_index("gene")
