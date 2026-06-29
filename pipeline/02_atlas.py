"""Stage 2: Atlas construction with lncRNA-forced HVG, Harmony, Leiden, annotation.

Input:  results/liver_anchor_clean.h5ad
Output: results/liver_anchor_atlas.h5ad
"""
from __future__ import annotations

import sys

import numpy as np
import scanpy as sc

from pipeline import config, lncrna_utils

sc.settings.verbosity = 2
np.random.seed(config.RANDOM_SEED)


def normalize_log(a):
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def select_hvg_with_lncrna_forced(a):
    sc.pp.highly_variable_genes(a, n_top_genes=config.ATLAS["n_hvg"], subset=False)
    if config.ATLAS["force_include_lncrna"]:
        lncrna_ids = set(lncrna_utils.get_lncrna_ids())
        lncrna_names = set(lncrna_utils.get_lncrna_name_mapping().values())
        var_match = (
            a.var_names.isin(lncrna_ids)
            | a.var_names.isin(lncrna_names)
        )
        a.var["highly_variable"] = a.var["highly_variable"] | var_match
        print(f"  forced lncRNA inclusion: +{int(var_match.sum())} genes")
    print(f"  total HVG: {int(a.var['highly_variable'].sum())}")
    return a


def harmony_correct(a):
    # Avoid sc.pp.scale full-matrix densification (152k cells * 24k genes ~14GB).
    # Use zero_center=False on sparse data — recommended for large sc atlases.
    sc.tl.pca(a, n_comps=config.ATLAS["n_pcs"], use_highly_variable=True, zero_center=False)
    batch_key = config.ATLAS["harmony_batch_key"]
    if batch_key not in a.obs.columns:
        # Try fallback names
        for alt in ("sample", "sample_id", "patient", "patient_id"):
            if alt in a.obs.columns:
                batch_key = alt
                print(f"  using fallback batch key: {batch_key}")
                break
    import harmonypy as hm
    ho = hm.run_harmony(a.obsm["X_pca"], a.obs, batch_key, random_state=config.RANDOM_SEED)
    # Detect the right orientation depending on harmonypy version
    z = ho.Z_corr
    n_cells = a.n_obs
    if z.ndim == 2 and z.shape[0] == n_cells:
        corrected = z
    elif z.ndim == 2 and z.shape[1] == n_cells:
        corrected = z.T
    else:
        raise RuntimeError(f"Unexpected ho.Z_corr shape {z.shape} for {n_cells} cells")
    a.obsm["X_pca_harmony"] = corrected
    print(f"  X_pca_harmony shape: {corrected.shape}")
    return a, batch_key


def cluster_and_umap(a):
    sc.pp.neighbors(a, use_rep="X_pca_harmony", n_neighbors=config.ATLAS["n_neighbors"])
    sc.tl.leiden(a, resolution=config.ATLAS["leiden_resolution"], random_state=config.RANDOM_SEED)
    sc.tl.umap(a, random_state=config.RANDOM_SEED)
    return a


# Marker genes per lineage — sourced from Guilliams 2022 + standard liver scRNA-seq references
LINEAGE_MARKERS = {
    "Hepatocyte": ["ALB", "TF", "APOB", "CYP3A4", "HNF4A"],
    "HSC": ["LRAT", "ACTA2", "COL1A1", "COL3A1", "PDGFRB"],
    "Kupffer": ["CD68", "MARCO", "VSIG4", "CD163", "TIMD4"],
    "LSEC": ["PECAM1", "STAB2", "CLEC4G", "OIT3", "FCN3"],
    "Cholangiocyte": ["KRT19", "KRT7", "EPCAM", "SOX9", "CFTR"],
    "Tcell": ["CD3D", "CD3E", "CD8A", "CD4"],
    "Bcell": ["CD79A", "MS4A1", "IGHM"],
    "NK": ["KLRD1", "NKG7", "GNLY"],
}


def annotate_lineages(a):
    """Score each cluster against marker sets, assign top-scoring lineage."""
    for name, markers in LINEAGE_MARKERS.items():
        valid = [m for m in markers if m in a.var_names]
        if not valid:
            continue
        sc.tl.score_genes(a, valid, score_name=f"score_{name}")
    score_cols = [f"score_{n}" for n in LINEAGE_MARKERS if f"score_{n}" in a.obs.columns]
    cluster_scores = a.obs.groupby("leiden")[score_cols].mean()
    cluster_to_lineage = cluster_scores.idxmax(axis=1).str.replace("score_", "")
    a.obs["lineage"] = a.obs["leiden"].map(cluster_to_lineage).astype("category")
    return a


def main():
    f_in = config.RESULTS_DIR / "liver_anchor_clean.h5ad"
    a = sc.read_h5ad(f_in)
    print(f"Loaded {f_in} shape={a.shape}")

    a = normalize_log(a)
    a = select_hvg_with_lncrna_forced(a)
    a, batch_key = harmony_correct(a)
    a = cluster_and_umap(a)
    a = annotate_lineages(a)

    out = config.RESULTS_DIR / "liver_anchor_atlas.h5ad"
    a.write_h5ad(out)
    print(f"Saved {out}")
    print(a.obs["lineage"].value_counts())
    return 0


if __name__ == "__main__":
    sys.exit(main())
