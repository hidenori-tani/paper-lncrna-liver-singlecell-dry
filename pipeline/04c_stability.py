"""Stage 4c: Join BRIC-seq stability class with pseudotime peak timing."""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc

from pipeline import config, lncrna_utils, stability_core


def find_peak_pseudotime(a, gene, lineage, n_bins=10):
    """Return pseudotime bin index where this gene peaks in this lineage."""
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
    return int(means.idxmax())


def main():
    a = sc.read_h5ad(config.RESULTS_DIR / "liver_anchor_pseudotime.h5ad")
    progression = pd.read_csv(config.RESULTS_DIR / "progression_lncrna_table.csv")
    progression = progression[progression["is_progression_coupled"]]

    half = pd.read_csv(config.DATA_DIR / "tani2012" / "half_lives.csv")
    half_classified = stability_core.classify(half)
    name_to_class = dict(zip(half_classified["gene_symbol"].str.upper(), half_classified["stability_class"]))

    id_to_name = lncrna_utils.get_lncrna_name_mapping()
    rows = []
    for _, r in progression.iterrows():
        gene_id = r["gene"]
        gene_name = id_to_name.get(gene_id, gene_id)
        stab = name_to_class.get(gene_name.upper(), "unknown")
        peak = find_peak_pseudotime(a, gene_id if gene_id in a.var_names else gene_name, r["lineage"], n_bins=config.SCREENING["n_bins"])
        rows.append({
            "gene_id": gene_id,
            "gene_name": gene_name,
            "lineage": r["lineage"],
            "stability_class": stab,
            "peak_bin": peak,
            "spearman_rho": r["spearman_rho"],
            "log2fc_max_min": r["log2fc_max_min"],
        })

    out = pd.DataFrame(rows)
    out.to_csv(config.RESULTS_DIR / "stability_x_pseudotime.csv", index=False)
    print(out.groupby(["stability_class", "lineage"])["peak_bin"].agg(["count", "median"]).to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
