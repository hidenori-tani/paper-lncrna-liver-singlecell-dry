"""Stage 5: Project validation cells onto anchor atlas, replicate progression trajectories.

Input:  results/liver_anchor_pseudotime.h5ad (anchor)
        data/GSE136103/ramachandran_cirrhosis.h5ad (validation)
Output: results/validation_projection.csv
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr

from pipeline import config, lncrna_utils


def main():
    anchor = sc.read_h5ad(config.RESULTS_DIR / "liver_anchor_pseudotime.h5ad")
    val = sc.read_h5ad(config.VALIDATION_PRIMARY_DIR / "ramachandran_cirrhosis.h5ad")
    print(f"anchor: {anchor.shape}, val: {val.shape}")

    try:
        import symphonypy as sp
        sp.tl.map_embedding(adata_ref=anchor, adata_query=val, use_rep="X_pca_harmony")
        print("symphonypy: projection done")
    except Exception as e:
        print(f"symphonypy unavailable ({e}); fallback to KNN nearest-cell in PCA space")

    from sklearn.neighbors import NearestNeighbors
    if "X_pca_harmony" in val.obsm:
        nn = NearestNeighbors(n_neighbors=15).fit(anchor.obsm["X_pca_harmony"])
        _, idx = nn.kneighbors(val.obsm["X_pca_harmony"])
        val.obs["transferred_lineage"] = anchor.obs["lineage"].values[idx[:, 0]]
        val.obs["transferred_pseudotime"] = np.nanmedian(anchor.obs["pseudotime"].values[idx], axis=1)

    progression = pd.read_csv(config.RESULTS_DIR / "progression_lncrna_table.csv")
    progression = progression[progression["is_progression_coupled"]]
    id_to_name = lncrna_utils.get_lncrna_name_mapping()

    rows = []
    for _, r in progression.iterrows():
        gid = r["gene"]
        gname = id_to_name.get(gid, gid)
        target = gid if gid in val.var_names else (gname if gname in val.var_names else None)
        if target is None:
            continue
        mask = (val.obs.get("transferred_lineage") == r["lineage"]) & ~val.obs["transferred_pseudotime"].isna()
        if mask.sum() < 100:
            continue
        gix = val.var_names.get_loc(target)
        expr = val.X[mask.values, gix].toarray().ravel() if hasattr(val.X, "toarray") else np.asarray(val.X[mask.values, gix]).ravel()
        pt = val.obs.loc[mask, "transferred_pseudotime"].values
        rho_v, _ = spearmanr(expr, pt)
        rows.append({
            "gene": gid,
            "gene_name": gname,
            "lineage": r["lineage"],
            "anchor_rho": r["spearman_rho"],
            "val_rho": rho_v,
            "replicated": (np.sign(rho_v) == np.sign(r["spearman_rho"])) and (abs(rho_v) > 0.15),
        })

    out = pd.DataFrame(rows)
    out.to_csv(config.RESULTS_DIR / "validation_projection.csv", index=False)
    n_total = len(out)
    n_repl = int(out["replicated"].sum())
    print(f"\nReplication: {n_repl} / {n_total} = {100 * n_repl / max(n_total, 1):.1f}%")
    return 0


if __name__ == "__main__":
    sys.exit(main())
