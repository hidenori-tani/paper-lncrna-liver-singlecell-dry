#!/usr/bin/env python3
"""Regenerate Fig 6 with A/B panel labels from the per-donor detection table (no mtx needed)."""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pathlib import Path

OUT = Path(__file__).resolve().parent.parent / "results" / "donor_level"
d = pd.read_csv(OUT / "per_donor_detection.csv")
LSEC = "Endothelial cells"
genes = ["KCNQ1OT1", "MEG3"]
panel = ["A", "B"]
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
for ax, gene, lab in zip(axes, genes, panel):
    sub = d[(d.gene == gene) & (d.n_cells >= 20)]
    order = sub.groupby("cell_type")["detection_rate"].median().sort_values(ascending=False).index.tolist()
    for i, ct in enumerate(order):
        vals = sub[sub.cell_type == ct]["detection_rate"].values
        x = np.full(len(vals), i, dtype=float) + np.linspace(-0.12, 0.12, len(vals))
        color = "#d62728" if ct == LSEC else "#1f77b4"
        ax.scatter(x, vals, s=20, alpha=0.75, color=color, edgecolors="none")
        ax.hlines(np.median(vals), i - 0.28, i + 0.28, color="black", lw=1.6)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=7)
    ax.set_ylabel("per-donor detection rate")
    ax.set_title(f"{gene}", fontsize=11, style="italic")
    ax.set_ylim(-0.02, 1.02)
    ax.text(-0.08, 1.06, lab, transform=ax.transAxes, fontsize=14, fontweight="bold", va="top")
fig.tight_layout()
fig.savefig(OUT / "Figure_6_donor_validation.pdf")
fig.savefig(OUT / "Figure_6_donor_validation.png", dpi=200)
print("wrote Figure_6_donor_validation.{pdf,png} with A/B panels")
