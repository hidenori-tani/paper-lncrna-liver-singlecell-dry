"""Re-annotate the anchor h5ad with proper Cell Ranger aggr suffix mapping.

Background: matrix barcodes use Cell Ranger aggr suffix (`-N` where N is the
sample index in `sampleComp_humanAll.txt` order). annot_humanAll.csv has
a `sample` column and a `cell` column with per-library `-1` suffix.

Correct join key: (sample_name, barcode_without_suffix).

Input:  data/GSE192742/human_allcells.h5ad         (from 00_build_anchor_h5ad.py)
        data/GSE192742/annot_humanAll.csv
        data/GSE192742/rawData_human_extracted/rawData_human/sampleComp_humanAll.txt
Output: data/GSE192742/human_allcells.h5ad         (overwritten with .obs populated)
"""
from __future__ import annotations

import re
import sys

import anndata as ad
import pandas as pd

from pipeline import config


def parse_sample_order(path) -> dict[int, str]:
    """Parse `- SAMPLENAME: URL` lines in sampleComp_humanAll.txt.

    Returns dict mapping 1-based aggr index to sample name.
    """
    samples = []
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^\s*-\s+([A-Za-z0-9]+):", line)
            if m:
                samples.append(m.group(1))
    return {i + 1: s for i, s in enumerate(samples)}


def main():
    anchor_dir = config.ANCHOR_DIR
    h5ad_path = anchor_dir / "human_allcells.h5ad"
    annot_path = anchor_dir / "annot_humanAll.csv"
    samplecomp_path = anchor_dir / "rawData_human_extracted" / "rawData_human" / "sampleComp_humanAll.txt"

    print(f"Loading {h5ad_path}...")
    a = ad.read_h5ad(h5ad_path)
    print(f"  shape: {a.shape}")

    print(f"Parsing sample order from {samplecomp_path}...")
    sample_order = parse_sample_order(samplecomp_path)
    print(f"  {len(sample_order)} samples in aggr order: 1->{sample_order[1]}, ..., {len(sample_order)}->{sample_order[len(sample_order)]}")

    # Extract aggr suffix and barcode root from matrix barcode
    bc_series = pd.Series(a.obs_names, index=a.obs_names)
    split = bc_series.str.rsplit("-", n=1, expand=True)
    bc_root = split[0]
    aggr_idx = pd.to_numeric(split[1], errors="coerce")
    print(f"  aggr suffix distribution: min={int(aggr_idx.min())}, max={int(aggr_idx.max())}, unique={aggr_idx.nunique()}")

    a.obs["aggr_idx"] = aggr_idx.values
    a.obs["sample"] = aggr_idx.map(sample_order).values
    a.obs["barcode_root"] = bc_root.values

    n_unmapped = a.obs["sample"].isna().sum()
    print(f"  cells without sample mapping: {n_unmapped} ({100*n_unmapped/a.n_obs:.2f}%)")

    print(f"Loading {annot_path}...")
    annot = pd.read_csv(annot_path)
    print(f"  annot rows: {len(annot)}, columns: {list(annot.columns)}")

    # Strip per-library suffix from annot cell barcode for consistent matching
    annot["barcode_root"] = annot["cell"].str.rsplit("-", n=1).str[0]

    # Join key: (sample, barcode_root)
    annot_indexed = annot.set_index(["sample", "barcode_root"])
    a_obs_keys = pd.MultiIndex.from_arrays([a.obs["sample"].values, a.obs["barcode_root"].values],
                                            names=["sample", "barcode_root"])
    print("Joining annot on (sample, barcode_root)...")
    annot_aligned = annot_indexed.reindex(a_obs_keys)
    annot_aligned.index = a.obs_names

    n_matched = annot_aligned["cell"].notna().sum()
    print(f"  matched: {n_matched} ({100*n_matched/a.n_obs:.2f}%) — annot retains only QC-passing cells")
    print(f"  expected: ~167599 (annot row count) out of {a.n_obs}")

    # Merge — keep matrix cells (left), attach annot fields
    for col in ["cluster", "annot", "patient", "digest", "typeSample", "diet", "UMAP_1", "UMAP_2"]:
        if col in annot_aligned.columns:
            a.obs[col] = annot_aligned[col].values

    # Standardize obs column names for downstream pipeline
    if "patient" in a.obs.columns:
        a.obs["donor_id"] = a.obs["patient"].astype("string")
    if "diet" in a.obs.columns:
        a.obs["condition"] = a.obs["diet"].astype("string")  # Lean / Obese
    if "annot" in a.obs.columns:
        a.obs["cell_type"] = a.obs["annot"].astype("string")

    # Filter to QC-passing cells (those Guilliams' team retained in annot_humanAll.csv)
    keep_mask = a.obs["cell_type"].notna()
    n_before = a.n_obs
    a = a[keep_mask].copy()
    print(f"Filtered to QC-passing cells (Guilliams team annotation): {n_before} -> {a.n_obs}")

    print("Final obs columns:", list(a.obs.columns))
    print("Donor distribution (top 20):")
    print(a.obs["donor_id"].value_counts().head(20))
    print("Condition distribution:")
    print(a.obs["condition"].value_counts(dropna=False))

    print(f"Saving {h5ad_path}...")
    a.write_h5ad(h5ad_path)
    print("Done.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
