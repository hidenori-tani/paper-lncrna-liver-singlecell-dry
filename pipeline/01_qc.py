"""Stage 1: QC for Guilliams 2022 anchor dataset.

Input:  data/GSE192742/human_allcells.h5ad  (produced by pipeline/00_build_anchor_h5ad.py)
Output: results/liver_anchor_clean.h5ad
        results/qc_report.csv
"""
from __future__ import annotations

import sys

import numpy as np
import pandas as pd
import scanpy as sc

from pipeline import config

sc.settings.verbosity = 2
np.random.seed(config.RANDOM_SEED)


def load_raw() -> "sc.AnnData":
    f = config.ANCHOR_DIR / "human_allcells.h5ad"
    print(f"Loading {f}...")
    a = sc.read_h5ad(f)
    print(f"  shape: {a.shape}")
    return a


def calc_qc_metrics(a):
    a.var["mt"] = a.var_names.str.upper().str.startswith("MT-")
    a.var["ribo"] = a.var_names.str.upper().str.startswith(("RPS", "RPL"))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt", "ribo"], inplace=True, log1p=False)
    return a


def filter_cells_genes(a):
    n_before = a.n_obs
    sc.pp.filter_cells(a, min_genes=config.QC["min_genes"])
    a = a[a.obs["n_genes_by_counts"] < config.QC["max_genes"]].copy()
    a = a[a.obs["pct_counts_mt"] < config.QC["mt_pct_max"]].copy()
    sc.pp.filter_genes(a, min_cells=3)
    print(f"  filter: {n_before} -> {a.n_obs} cells")
    return a


def doublet_removal(a):
    """Use scanpy's built-in scrublet if available; fallback to no-op + flag."""
    try:
        sc.pp.scrublet(a, batch_key=None)
        a = a[~a.obs["predicted_doublet"]].copy()
        print(f"  doublet removal: {a.n_obs} cells remaining")
    except Exception as e:
        print(f"  WARN doublet removal skipped: {e}")
        a.obs["predicted_doublet"] = False
    return a


def write_qc_report(a, path):
    donor_col = None
    for c in ("donor_id", "sample", "patient", "patient_id"):
        if c in a.obs.columns:
            donor_col = c
            break
    summary = pd.DataFrame({
        "metric": [
            "n_cells", "n_genes", "median_n_genes_per_cell",
            "median_total_counts_per_cell", "median_pct_mt", "n_donors",
        ],
        "value": [
            a.n_obs,
            a.n_vars,
            float(np.median(a.obs["n_genes_by_counts"])),
            float(np.median(a.obs["total_counts"])),
            float(np.median(a.obs["pct_counts_mt"])),
            a.obs[donor_col].nunique() if donor_col else 0,
        ],
    })
    summary.to_csv(path, index=False)
    print(summary.to_string(index=False))


def main():
    config.RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    a = load_raw()
    a = calc_qc_metrics(a)
    a = filter_cells_genes(a)
    a = doublet_removal(a)
    out = config.RESULTS_DIR / "liver_anchor_clean.h5ad"
    a.write_h5ad(out)
    write_qc_report(a, config.RESULTS_DIR / "qc_report.csv")
    print(f"Saved {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
