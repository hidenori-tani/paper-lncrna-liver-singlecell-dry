"""Single source of truth for paths, seeds, thresholds.

All pipeline scripts MUST import constants from here. Never hardcode.
"""
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
RESULTS_DIR = PROJECT_ROOT / "results"
FIGURES_DIR = PROJECT_ROOT / "figures" / "output"
MANUSCRIPT_DIR = PROJECT_ROOT / "manuscript"

RANDOM_SEED = 42

# Anchor dataset
ANCHOR_GSE = "GSE192742"
ANCHOR_DIR = DATA_DIR / ANCHOR_GSE

# Validation
VALIDATION_PRIMARY_GSE = "GSE136103"  # Ramachandran 2019 Nature
VALIDATION_PRIMARY_DIR = DATA_DIR / VALIDATION_PRIMARY_GSE
VALIDATION_NASH_DIR = DATA_DIR / "NASH_severity"  # selected in Stage 1

# Annotation
GENCODE_VERSION = "v45"
GENCODE_DIR = DATA_DIR / "gencode"

# QC thresholds (Stage 1)
QC = {
    "min_genes": 200,
    "max_genes": 8000,
    "mt_pct_max": 20,
    "doublet_score_max": 0.25,  # scDblFinder threshold
}

# Atlas construction (Stage 2)
ATLAS = {
    "n_hvg": 4000,
    "force_include_lncrna": True,
    "harmony_batch_key": "donor_id",
    "leiden_resolution": 1.0,
    "n_pcs": 50,
    "n_neighbors": 15,
}

# Pseudotime (Stage 3)
PSEUDOTIME = {
    "candidates": ["monocle3", "paga_dpt"],  # decided after subset comparison
    "lineages": ["Hepatocyte", "HSC", "Kupffer", "LSEC", "Cholangiocyte"],
    "root_strategy": "healthy_hepatocyte_median",
}

# Screening (Stage 4a)
SCREENING = {
    "spearman_abs_min": 0.3,
    "wald_q_max": 0.05,
    "log2fc_abs_min": 1.0,
    "min_cells_per_bin": 50,
    "n_bins": 10,
}

# Sensitivity sweep grid for Supp Fig
SCREENING_SENSITIVITY = {
    "spearman_abs_min": [0.2, 0.3, 0.4],
    "log2fc_abs_min": [0.5, 1.0, 1.5],
}

# Stability classes (Tani 2012 Genome Res, h)
STABILITY_CLASSES = {
    "short_max_h": 4.0,
    "long_min_h": 12.0,
    # short:    t1/2 < 4h
    # medium:   4h <= t1/2 < 12h
    # long:     t1/2 >= 12h
}

# Anchor lncRNA panel for focused deep-dive (Stage 4b)
ANCHOR_LNCRNAS = [
    # Author's wet-lab anchors
    "KCNQ1OT1",
    "RMST",
    "IDI2-AS1",
    # Bradford Hill review Tier 1
    "NEAT1",
    "MALAT1",
    "MEG3",
    # Canonical liver-relevant
    "HOTAIR",
    "H19",
    "HULC",
    "MIAT",
    "LINC00261",
    "DANCR",
    "LINC-GINS2",
]
