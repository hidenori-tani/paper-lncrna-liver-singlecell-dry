#!/usr/bin/env python3
"""Fig 2 (per-lineage screening volcano) and Fig 5 (peak-bin x BRIC-seq stability),
using neutral 'screen-flagged' terminology, regenerated from the deposited result
tables. Self-contained; reads results/, writes figures/output/.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJ = Path(__file__).resolve().parent.parent  # repo root
RES = PROJ / "results"
OUT = PROJ / "figures" / "output"
OUT.mkdir(parents=True, exist_ok=True)

RHO_THR, Q_THR = 0.3, 0.05
STAB_COLORS = {"short": "#d62728", "medium": "#ff7f0e", "long": "#2ca02c", "unknown": "#7f7f7f"}

# ---- Fig 2: per-lineage screening volcano ----
df = pd.read_csv(RES / "progression_lncrna_table.csv")
lineages, letters = ["Hepatocyte", "Kupffer", "LSEC", "Cholangiocyte"], ["A", "B", "C", "D"]
fig, axes = plt.subplots(2, 2, figsize=(11, 9))
for ax, lin, L in zip(axes.ravel(), lineages, letters):
    sub = df[df["lineage"] == lin].copy()
    sub["neglogq"] = -np.log10(sub["wald_q"].clip(lower=1e-300))
    fl = sub["is_progression_coupled"] == True
    ax.scatter(sub.loc[~fl, "spearman_rho"], sub.loc[~fl, "neglogq"], s=8,
               color="#bdbdbd", alpha=0.5, edgecolors="none")
    ax.scatter(sub.loc[fl, "spearman_rho"], sub.loc[fl, "neglogq"], s=42,
               color="#d62728", edgecolors="black", linewidths=0.4, zorder=5)
    for _, r in sub[fl].iterrows():
        ax.annotate(r["gene"], (r["spearman_rho"], -np.log10(max(r["wald_q"], 1e-300))),
                    fontsize=8, fontweight="bold", xytext=(4, 3), textcoords="offset points")
    ax.axvline(RHO_THR, ls="--", lw=0.7, color="0.4"); ax.axvline(-RHO_THR, ls="--", lw=0.7, color="0.4")
    ax.axhline(-np.log10(Q_THR), ls="--", lw=0.7, color="0.4")
    ax.set_title(f"({L}) {lin}   screen-flagged: n = {int(fl.sum())}", fontsize=10)
    ax.set_xlabel("Spearman ρ (expression vs pseudotime)", fontsize=8)
    ax.set_ylabel("−log₁₀(Benjamini–Hochberg q)", fontsize=8)
    ax.set_xlim(-0.65, 0.65); ax.tick_params(labelsize=7)
fig.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor="#d62728",
                    markeredgecolor="black", markersize=8, label="screen-flagged lncRNA"),
                    Line2D([0], [0], marker="o", color="w", markerfacecolor="#bdbdbd",
                    markersize=6, label="not flagged")],
           loc="lower center", ncol=2, fontsize=9, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("Per-lineage screening for lncRNAs correlated with within-lineage pseudotime", fontsize=11)
fig.tight_layout(rect=[0, 0.03, 1, 0.97])
fig.savefig(OUT / "fig2_screening.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig2_screening.png", dpi=300, bbox_inches="tight")
plt.close(fig)

# ---- Fig 5: peak bin x stability ----
s = pd.read_csv(RES / "stability_x_pseudotime.csv").dropna(subset=["peak_bin"])
lin_order = ["LSEC", "Cholangiocyte"]
ypos = {lin: i for i, lin in enumerate(lin_order[::-1])}
fig, ax = plt.subplots(figsize=(8.5, 3.6))
for _, r in s.iterrows():
    ax.scatter(r["peak_bin"], ypos[r["lineage"]], s=180,
               color=STAB_COLORS.get(r["stability_class"], "#7f7f7f"),
               edgecolors="black", linewidths=0.6, zorder=5)
    ax.annotate(r["gene_name"], (r["peak_bin"], ypos[r["lineage"]]), fontsize=9,
                fontweight="bold", ha="center", va="bottom", xytext=(0, 12), textcoords="offset points")
ax.set_yticks(list(ypos.values())); ax.set_yticklabels(list(ypos.keys()))
ax.set_ylim(-0.6, len(lin_order) - 0.4); ax.set_xlim(-0.6, 9.6); ax.set_xticks(range(10))
ax.set_xlabel("Pseudotime peak bin (0 = earliest, 9 = latest)")
ax.set_title("Peak pseudotime bin and BRIC-seq stability of the screen-flagged lncRNAs", fontsize=10)
labels = {"short": "short (t₁/₂ < 4 h)", "medium": "medium (4–12 h)",
          "long": "long (≥ 12 h)", "unknown": "unknown (absent from BRIC-seq)"}
ax.legend(handles=[Line2D([0], [0], marker="o", color="w", markerfacecolor=STAB_COLORS[k],
                   markeredgecolor="black", markersize=9, label=labels[k])
                   for k in ["short", "medium", "long", "unknown"]],
          title="BRIC-seq stability class", fontsize=8, title_fontsize=8,
          loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
ax.grid(axis="x", ls=":", alpha=0.4)
fig.tight_layout()
fig.savefig(OUT / "fig5_stability.pdf", bbox_inches="tight")
fig.savefig(OUT / "fig5_stability.png", dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote fig2_screening + fig5_stability to figures/output/")
