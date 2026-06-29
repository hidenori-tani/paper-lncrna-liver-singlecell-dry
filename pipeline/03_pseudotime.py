"""Stage 3: Per-lineage pseudotime on full atlas using chosen tool from config.

Input:  results/liver_anchor_atlas.h5ad
Output: results/liver_anchor_pseudotime.h5ad
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc

from pipeline import config

np.random.seed(config.RANDOM_SEED)


def find_root(a_lineage):
    """Healthy-side median cell of the lineage as pseudotime root."""
    cond_key = next((k for k in ("condition", "disease", "diagnosis", "health_status") if k in a_lineage.obs.columns), None)
    if cond_key is None:
        return 0
    healthy_mask = a_lineage.obs[cond_key].astype(str).str.lower().str.contains("health|control|normal|lean|hca")
    if healthy_mask.sum() == 0:
        return 0
    cands = np.where(healthy_mask)[0]
    return int(np.random.choice(cands))


def pseudotime_paga_dpt(a_lineage):
    """Per-lineage DPT pseudotime. PAGA is skipped because per-lineage subsets
    often have insufficient leiden clusters for PAGA's connectivity matrix."""
    sc.pp.neighbors(a_lineage, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.diffmap(a_lineage)
    a_lineage.uns["iroot"] = find_root(a_lineage)
    sc.tl.dpt(a_lineage)
    return a_lineage.obs["dpt_pseudotime"].values


PSEUDOTIME_FUNCS = {
    "paga_dpt": pseudotime_paga_dpt,
    # "monocle3": run_monocle3_via_r (import from 03_pseudotime_compare.py if chosen)
}


def main():
    a = sc.read_h5ad(config.RESULTS_DIR / "liver_anchor_atlas.h5ad")
    tool = config.PSEUDOTIME.get("tool", "paga_dpt")
    func = PSEUDOTIME_FUNCS[tool]

    pt = pd.Series(np.full(a.n_obs, np.nan), index=a.obs_names)
    for lin in config.PSEUDOTIME["lineages"]:
        mask = a.obs["lineage"] == lin
        n = int(mask.sum())
        if n < 200:
            print(f"  skipping {lin}: only {n} cells")
            continue
        sub = a[mask].copy()
        pt_sub = func(sub)
        pt.loc[sub.obs_names] = pt_sub
        print(f"  {lin}: {n} cells, pseudotime range [{np.nanmin(pt_sub):.3f}, {np.nanmax(pt_sub):.3f}]")

    a.obs["pseudotime"] = pt.values
    out = config.RESULTS_DIR / "liver_anchor_pseudotime.h5ad"
    a.write_h5ad(out)
    print(f"Saved {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
