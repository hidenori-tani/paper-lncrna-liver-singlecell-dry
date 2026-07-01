#!/usr/bin/env python3
"""Independent-cohort replication of the MEG3-in-LSEC finding in MacParland 2018.

Primary cohort: Guilliams 2022 human Liver Cell Atlas (GSE192742, 16 donors) — MEG3 is
donor-reproducibly enriched in liver sinusoidal endothelial cells (LSECs).

This script tests the SAME claim in an INDEPENDENT human liver scRNA-seq atlas
(MacParland et al. 2018, Nat Commun; GEO GSE115469; 5 donors, ~8,444 cells), using the
authors' OWN published cell-type annotation. Statistics are computed at the DONOR level
(n = 5 donors), mirroring the primary analysis (per-donor detection rate + paired Wilcoxon).

Same dependency footprint as scripts/donor_level_analysis.py:
  only gzip, numpy, pandas, scipy, statsmodels, matplotlib. No scanpy / Seurat / R.

Inputs (data/GSE115469/ — fetch from GEO, not redistributed here):
  GSE115469_Data.csv.gz             log-normalized expression, genes x cells (20,007 x 8,444)
  GSE115469_CellClusterType.txt.gz  per-cell annotation: CellName, Sample(=donor), CellType

  Fetch with:
    mkdir -p data/GSE115469 && cd data/GSE115469
    curl -sL -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/GSE115469_Data.csv.gz
    curl -sL -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/GSE115469_CellClusterType.txt.gz

Outputs (results/replication_macparland/):
  macparland_per_donor_detection.csv   donor x cell_type x gene: detection rate + cell counts
  macparland_lsec_enrichment.csv       per gene: LSEC vs pooled-other, per-donor + Wilcoxon
  macparland_celltype_ranking.csv      per gene: median per-donor detection by cell type
  fig_macparland_replication.pdf/png   per-donor detection strip plots (MEG3, KCNQ1OT1)
"""
from __future__ import annotations
import gzip, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent  # repo root
DATA = PROJ / "data" / "GSE115469"
CSV = DATA / "GSE115469_Data.csv.gz"
ANNOT = DATA / "GSE115469_CellClusterType.txt.gz"
OUTDIR = PROJ / "results" / "replication_macparland"
OUTDIR.mkdir(parents=True, exist_ok=True)

# lncRNAs of interest (primary-cohort headline + honest controls)
TARGETS = ["MEG3", "KCNQ1OT1", "ZFAS1", "NEAT1", "MALAT1"]
# canonical LSEC markers, used ONLY as an internal QC that "LSEC"-labelled cells are real LSECs
MARKERS = ["CLEC4G", "STAB2", "OIT3", "FCN3", "FCN2", "CLEC1B", "ALB"]
WANT = set(TARGETS) | set(MARKERS)

# LSEC = liver sinusoidal endothelial cells. MacParland resolves two LSEC zonation clusters.
# Portal_endothelial_Cells is a DISTINCT (non-sinusoidal) vascular endothelial type -> NOT pooled in.
LSEC_TYPES = {"Central_venous_LSECs", "Periportal_LSECs"}


def log(*a): print(*a, flush=True)


def load_matrix_rows(csv_path: Path, want: set):
    """Stream the dense genes x cells CSV; keep only rows whose gene symbol is in `want`."""
    data = {}
    with gzip.open(csv_path, "rt") as fh:
        header = fh.readline().rstrip("\n").split(",")
        cells = [c.strip().strip('"') for c in header[1:]]
        log(f"  matrix header: {len(cells)} cells")
        for ln in fh:
            comma = ln.find(",")
            gene = ln[:comma].strip().strip('"')
            if gene in want:
                vals = np.asarray(ln[comma + 1:].rstrip("\n").split(","), dtype=np.float32)
                data[gene] = vals
    found = sorted(data.keys())
    missing = sorted(want - set(found))
    log(f"  genes found ({len(found)}): {found}")
    if missing:
        log(f"  genes MISSING: {missing}")
    return pd.DataFrame(data, index=cells)


def main():
    log("[1] load annotation")
    ann = pd.read_csv(ANNOT, sep="\t")
    ann.columns = [c.strip() for c in ann.columns]
    ann = ann.rename(columns={"CellName": "cell", "Sample": "donor", "CellType": "cell_type"})
    ann = ann[["cell", "donor", "cell_type"]].set_index("cell")
    log(f"  donors: {sorted(ann['donor'].unique())} | cell types: {ann['cell_type'].nunique()}")

    log("[2] stream matrix (keep only target + marker rows)")
    mat = load_matrix_rows(CSV, WANT)

    log("[3] join expression to annotation on cell name")
    df = mat.join(ann, how="inner")
    df["is_lsec"] = df["cell_type"].isin(LSEC_TYPES)
    log(f"  cells with expression + annotation: {len(df)}")
    log(f"  LSEC cells per donor:\n{df[df['is_lsec']].groupby('donor').size().to_string()}")
    present_targets = [g for g in TARGETS if g in df.columns]

    log("[4] QC: LSEC marker detection (LSEC vs non-LSEC, pooled)")
    for m in [g for g in MARKERS if g in df.columns]:
        dl = (df.loc[df["is_lsec"], m] > 0).mean()
        do = (df.loc[~df["is_lsec"], m] > 0).mean()
        log(f"    {m:8s} LSEC {dl:5.1%}  non-LSEC {do:5.1%}")

    log("[5] per-donor detection rates by cell type")
    rows = []
    for (donor, ct), g in df.groupby(["donor", "cell_type"]):
        for gene in present_targets:
            n = len(g); ndet = int((g[gene] > 0).sum())
            rows.append(dict(donor=donor, cell_type=ct, gene=gene,
                             n_cells=n, n_detected=ndet, detection_rate=ndet / n))
    perdonor = pd.DataFrame(rows)
    perdonor.to_csv(OUTDIR / "macparland_per_donor_detection.csv", index=False)

    log("[6] cell-type ranking by median per-donor detection")
    rankrows = []
    for gene in present_targets:
        sub = perdonor[(perdonor["gene"] == gene) & (perdonor["n_cells"] >= 10)]
        med = sub.groupby("cell_type")["detection_rate"].median().sort_values(ascending=False)
        for rank, (ct, val) in enumerate(med.items(), 1):
            rankrows.append(dict(gene=gene, rank=rank, cell_type=ct, median_detection=val))
    pd.DataFrame(rankrows).to_csv(OUTDIR / "macparland_celltype_ranking.csv", index=False)

    log("[7] donor-level LSEC enrichment (paired one-sided Wilcoxon across 5 donors)")
    enr = []
    donors = sorted(df["donor"].unique())
    for gene in present_targets:
        lsec_rate, other_rate, kept = [], [], []
        for donor in donors:
            gd = df[df["donor"] == donor]
            lc = gd[gd["is_lsec"]]; oc = gd[~gd["is_lsec"]]
            if len(lc) < 10:
                continue
            lsec_rate.append((lc[gene] > 0).mean()); other_rate.append((oc[gene] > 0).mean()); kept.append(donor)
        lsec_rate = np.array(lsec_rate); other_rate = np.array(other_rate); n = len(lsec_rate)
        n_higher = int(np.sum(lsec_rate > other_rate))
        if n >= 5 and not np.allclose(lsec_rate, other_rate):
            try:
                stat, p = wilcoxon(lsec_rate, other_rate, alternative="greater")
            except ValueError:
                stat, p = np.nan, np.nan
        else:
            stat, p = np.nan, np.nan
        enr.append(dict(gene=gene, n_donors=n, n_donors_LSEC_higher=n_higher,
                        median_LSEC=float(np.median(lsec_rate)) if n else np.nan,
                        median_nonLSEC=float(np.median(other_rate)) if n else np.nan,
                        wilcoxon_stat=stat, p_greater=p,
                        per_donor_LSEC=";".join(f"{d}:{r:.3f}" for d, r in zip(kept, lsec_rate)),
                        per_donor_nonLSEC=";".join(f"{d}:{r:.3f}" for d, r in zip(kept, other_rate))))
    enr = pd.DataFrame(enr)
    if enr["p_greater"].notna().any():
        enr["q_bh"] = np.nan
        m = enr["p_greater"].notna()
        enr.loc[m, "q_bh"] = multipletests(enr.loc[m, "p_greater"], method="fdr_bh")[1]
    enr.to_csv(OUTDIR / "macparland_lsec_enrichment.csv", index=False)
    log(enr[["gene", "n_donors", "n_donors_LSEC_higher", "median_LSEC", "median_nonLSEC", "p_greater"]].to_string(index=False))

    log("[8] figure")
    main_genes = [g for g in ["MEG3", "KCNQ1OT1"] if g in present_targets]
    fig, axes = plt.subplots(1, len(main_genes), figsize=(7.0 * len(main_genes), 5.0), squeeze=False)
    rng = np.random.default_rng(0)
    for ax, gene in zip(axes[0], main_genes):
        sub = perdonor[(perdonor["gene"] == gene) & (perdonor["n_cells"] >= 10)]
        order = sub.groupby("cell_type")["detection_rate"].median().sort_values(ascending=False).index.tolist()
        for i, ct in enumerate(order):
            vals = sub[sub["cell_type"] == ct]["detection_rate"].values
            ax.scatter(rng.normal(i, 0.06, size=len(vals)), vals, s=22, alpha=0.75,
                       color="#d62728" if ct in LSEC_TYPES else "#1f77b4", edgecolors="none")
            ax.hlines(np.median(vals), i - 0.25, i + 0.25, color="black", lw=1.5)
        ax.set_xticks(range(len(order))); ax.set_xticklabels(order, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("per-donor detection rate"); ax.set_ylim(-0.02, 1.02)
        ax.set_title(f"{gene}  (MacParland 2018; each point = 1 donor)", fontsize=10)
    fig.suptitle("Independent replication (GSE115469, 5 donors): LSEC clusters in red", fontsize=9, y=1.02)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_macparland_replication.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "fig_macparland_replication.png", dpi=200, bbox_inches="tight")
    log("DONE.")


if __name__ == "__main__":
    sys.exit(main())
