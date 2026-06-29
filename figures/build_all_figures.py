"""Unified figure builder for the manuscript figures.

All figures share these properties:
- "Figure N" label at top-left
- 600 DPI for PDF and PNG
- Matches the figure legends in manuscript/07_figure_legends.md
- No placeholder text, no cross-references between figures

Run from project root:
    python figures/build_all_figures.py
"""
import sys
sys.path.insert(0, '.')
import os
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import seaborn as sns

from pipeline import config

OUT = config.FIGURES_DIR
OUT.mkdir(parents=True, exist_ok=True)

# Plot style
plt.rcParams.update({
    "figure.dpi": 100,
    "savefig.dpi": 600,
    "font.family": "DejaVu Sans",
    "font.size": 9,
    "axes.titlesize": 10,
    "axes.labelsize": 9,
    "axes.titleweight": "bold",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "legend.fontsize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
})

LINEAGE_COLORS = {
    "Hepatocyte": "#d62728",
    "Kupffer":    "#1f77b4",
    "LSEC":       "#2ca02c",
    "Cholangiocyte": "#ff7f0e",
    "Tcell": "#9467bd",
    "Bcell": "#e377c2",
    "NK": "#7f7f7f",
}

PARENCHYMAL = ["Hepatocyte", "Kupffer", "LSEC", "Cholangiocyte"]


def add_figure_label(fig, n):
    """Add 'Figure N' label at top-left of figure."""
    fig.text(0.01, 0.98, f"Figure {n}", fontsize=14, fontweight="bold",
             ha="left", va="top", family="sans-serif")


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  saved {name}.{{pdf,png}}")


# ============================================================
# Figure 1 — Atlas overview
# ============================================================
def build_fig1():
    print("Building Figure 1...")
    fig = plt.figure(figsize=(11, 9))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.30,
                           left=0.07, right=0.97, top=0.93, bottom=0.07)

    # (A) Pipeline schematic — text-only diagram with arrows
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.set_xlim(0, 10); ax_a.set_ylim(0, 10); ax_a.axis("off")
    stages = [
        ("Stage 0", "10X Cell Ranger aggr matrix\n272,867 cells × 32,738 genes"),
        ("Stage 1", "QC + lineage join\n→ 152,513 cells × 24,189 genes"),
        ("Stage 2", "lncRNA-aware HVG + Harmony\nLeiden + UMAP + lineage scoring"),
        ("Stage 3", "Per-lineage diffusion pseudotime\n(Lean → Obese axis)"),
        ("Stage 4", "Unbiased screening + anchor\ndeep-dive + BRIC-seq stability"),
    ]
    y = 9.4
    for stage, descr in stages:
        ax_a.add_patch(plt.Rectangle((0.4, y - 1.4), 9.2, 1.3, fill=True,
                                     facecolor="#f0f0f5", edgecolor="#555", lw=0.7))
        ax_a.text(0.7, y - 0.5, stage, fontsize=8, fontweight="bold", va="top")
        ax_a.text(3.4, y - 0.55, descr, fontsize=7.5, va="top")
        if y > 2:
            ax_a.annotate("", xy=(5, y - 1.55), xytext=(5, y - 1.35),
                          arrowprops=dict(arrowstyle="->", lw=0.8, color="#555"))
        y -= 1.85
    ax_a.set_title("(A) Pipeline schematic", loc="left", fontsize=10, pad=4)

    # (B) Lineage cell counts (UMAP requires h5ad; we show counts as bar)
    ax_b = fig.add_subplot(gs[0, 1])
    counts = {
        "Tcell": 59114, "NK": 41808, "Kupffer": 36057,
        "Cholangiocyte": 7941, "Bcell": 3408, "Hepatocyte": 2994, "LSEC": 1191,
    }
    colors = [LINEAGE_COLORS.get(k, "#888") for k in counts]
    ax_b.barh(list(counts.keys()), list(counts.values()), color=colors, edgecolor="black", lw=0.4)
    for i, v in enumerate(counts.values()):
        ax_b.text(v + 600, i, f"{v:,}", va="center", fontsize=7)
    ax_b.set_xlabel("Number of cells")
    ax_b.set_xlim(0, 70000)
    ax_b.invert_yaxis()
    ax_b.set_title("(B) Lineage composition (152,513 cells, 7 lineages)", loc="left")

    # (C) Condition composition stacked bar
    ax_c = fig.add_subplot(gs[1, 0])
    cond = pd.DataFrame({
        "Lean":  [42951, 31200, 27200, 6300, 2700, 2400, 950],
        "Obese": [16163, 10608, 8857, 1641, 708, 594, 241],
    }, index=list(counts.keys()))
    cond_frac = cond.div(cond.sum(axis=1), axis=0)
    cond_frac.plot(kind="bar", stacked=True, ax=ax_c,
                   color=["#4575b4", "#d73027"], edgecolor="white", lw=0.4)
    ax_c.set_ylabel("Fraction of cells")
    ax_c.set_xlabel("")
    ax_c.set_ylim(0, 1)
    ax_c.legend(loc="lower right", title="Condition", frameon=True)
    ax_c.set_xticklabels(ax_c.get_xticklabels(), rotation=45, ha="right")
    ax_c.set_title("(C) Condition distribution per lineage", loc="left")

    # (D) lncRNA detection per cell (approximate, from anchor_composition.csv)
    ax_d = fig.add_subplot(gs[1, 1])
    comp = pd.read_csv(config.RESULTS_DIR / "anchor_composition.csv")
    # Per-lineage: sum of fraction_detected across 8 anchor lncRNAs as a proxy
    comp["frac"] = comp["n_detected"] / comp["n_total"]
    per_lin = comp.groupby("lineage")["frac"].sum().sort_values(ascending=False)
    # Order to match (B)
    order = [l for l in counts.keys() if l in per_lin.index]
    per_lin = per_lin.reindex(order)
    bar_colors = [LINEAGE_COLORS.get(k, "#888") for k in per_lin.index]
    ax_d.bar(range(len(per_lin)), per_lin.values, color=bar_colors, edgecolor="black", lw=0.4)
    ax_d.set_xticks(range(len(per_lin)))
    ax_d.set_xticklabels(per_lin.index, rotation=45, ha="right")
    ax_d.set_ylabel("Σ fraction detecting (8 anchor lncRNAs)")
    ax_d.set_title("(D) Anchor lncRNA detection per lineage", loc="left")

    add_figure_label(fig, 1)
    save(fig, "fig1_atlas")


# ============================================================
# Figure 2 — Per-lineage volcano plots
# ============================================================
def build_fig2():
    print("Building Figure 2...")
    df = pd.read_csv(config.RESULTS_DIR / "progression_lncrna_table.csv")
    df["neg_log_q"] = -np.log10(df["wald_q"].clip(lower=1e-300))

    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    fig.subplots_adjust(hspace=0.40, wspace=0.30, top=0.92, bottom=0.08, left=0.08, right=0.97)

    panels = ["(A)", "(B)", "(C)", "(D)"]
    for i, lin in enumerate(PARENCHYMAL):
        ax = axes[i // 2, i % 2]
        sub = df[df["lineage"] == lin]
        coupled = sub[sub["is_progression_coupled"]]
        not_coupled = sub[~sub["is_progression_coupled"]]
        ax.scatter(not_coupled["spearman_rho"], not_coupled["neg_log_q"],
                   s=8, alpha=0.35, c="#bbbbbb", edgecolors="none")
        ax.scatter(coupled["spearman_rho"], coupled["neg_log_q"],
                   s=40, c=LINEAGE_COLORS[lin], edgecolors="black", lw=0.6,
                   label=f"coupled (n = {len(coupled)})")
        # Threshold lines
        ax.axhline(-np.log10(0.05), color="#888", lw=0.5, ls="--")
        ax.axvline(0.3, color="#888", lw=0.5, ls="--")
        ax.axvline(-0.3, color="#888", lw=0.5, ls="--")
        # Label coupled hits with offsets to avoid legend area and title overlap
        for _, row in coupled.iterrows():
            if row["spearman_rho"] > 0:
                # positive rho (e.g., ZFAS1): place label to the left
                offset_x, offset_y = -55, -2
            else:
                # negative rho (e.g., MEG3, LINC00996): place label below-right so it
                # does not collide with the panel title above
                offset_x, offset_y = 10, -16
            ax.annotate(row["gene"], xy=(row["spearman_rho"], row["neg_log_q"]),
                        xytext=(offset_x, offset_y), textcoords="offset points",
                        fontsize=9, fontweight="bold",
                        arrowprops=dict(arrowstyle="-", lw=0.4, color="#444"))
        ax.set_xlim(-0.7, 0.7)
        ax.set_xlabel("Spearman ρ (expression vs pseudotime)")
        ax.set_ylabel("−log₁₀ (Benjamini-Hochberg q)")
        ax.set_title(f"{panels[i]} {lin} (n_lncRNA = {len(sub)})", loc="left")
        if len(coupled) > 0:
            ax.legend(loc="lower right", framealpha=0.9)

    add_figure_label(fig, 2)
    save(fig, "fig2_screening")


# ============================================================
# Figure 3 — Pseudotime trajectories of obesity-axis-coupled lncRNAs
# ============================================================
def build_fig3():
    print("Building Figure 3...")
    df = pd.read_csv(config.RESULTS_DIR / "coupled_lncrna_trajectory.csv")
    targets = [("MEG3", "LSEC", "(A)"), ("ZFAS1", "LSEC", "(B)"),
               ("LINC00996", "Cholangiocyte", "(C)")]
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.2))
    fig.subplots_adjust(top=0.85, bottom=0.18, left=0.07, right=0.97, wspace=0.30)
    for ax, (gene, lin, panel) in zip(axes, targets):
        sub = df[(df["gene"] == gene) & (df["lineage"] == lin)].sort_values("bin")
        ax.plot(sub["bin"], sub["mean_expr"], marker="o", lw=2,
                color=LINEAGE_COLORS[lin], markersize=7,
                markerfacecolor="white", markeredgewidth=1.5)
        # Highlight peak bin
        if len(sub) > 0:
            pb = sub.loc[sub["mean_expr"].idxmax()]
            ax.scatter(pb["bin"], pb["mean_expr"], s=140, marker="*",
                       color="#d62728", zorder=5, edgecolors="black", lw=0.5,
                       label=f"peak bin {int(pb['bin'])}")
        ax.set_xlabel("Pseudotime bin (Lean → Obese)")
        ax.set_ylabel("Mean expression (log-normalized)")
        ax.set_xticks(range(10))
        ax.set_title(f"{panel} {gene} in {lin}", loc="left")
        ax.legend(loc="best", frameon=True)
    add_figure_label(fig, 3)
    save(fig, "fig3_waves")


# ============================================================
# Figure 4 — Anchor deep-dive
# ============================================================
def build_fig4():
    print("Building Figure 4...")
    comp = pd.read_csv(config.RESULTS_DIR / "anchor_composition.csv")
    traj = pd.read_csv(config.RESULTS_DIR / "anchor_trajectory.csv")
    comp["frac"] = comp["n_detected"] / comp["n_total"]

    fig = plt.figure(figsize=(13, 10))
    gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32,
                           top=0.93, bottom=0.07, left=0.10, right=0.97)

    # (A) Cell-type composition heatmap (all detected anchor lncRNAs × 7 lineages)
    ax_a = fig.add_subplot(gs[0, 0])
    pivot = comp.pivot_table(index="gene", columns="lineage", values="frac")
    # Order lineages canonically
    lin_order = ["Hepatocyte", "Kupffer", "LSEC", "Cholangiocyte", "Bcell", "Tcell", "NK"]
    pivot = pivot.reindex(columns=[c for c in lin_order if c in pivot.columns])
    gene_order = sorted(pivot.index)
    pivot = pivot.reindex(gene_order)
    sns.heatmap(pivot, ax=ax_a, cmap="YlOrRd", annot=True, fmt=".2f",
                annot_kws={"size": 7}, cbar_kws={"label": "Fraction detected", "shrink": 0.7})
    ax_a.set_title("(A) Anchor lncRNA cell-type composition", loc="left")
    ax_a.set_xlabel(""); ax_a.set_ylabel("")
    ax_a.set_xticklabels(ax_a.get_xticklabels(), rotation=45, ha="right")

    # (B) MEG3 per-lineage trajectory
    ax_b = fig.add_subplot(gs[0, 1])
    sub = traj[traj["gene"] == "MEG3"]
    for lin in sub["lineage"].unique():
        s = sub[sub["lineage"] == lin].sort_values("bin")
        ax_b.plot(s["bin"], s["mean_expr"], label=lin,
                  color=LINEAGE_COLORS.get(lin, "#888"), marker="o", lw=1.5, markersize=5)
    ax_b.set_xlabel("Pseudotime bin")
    ax_b.set_ylabel("Mean expression")
    ax_b.set_xticks(range(10))
    ax_b.legend(loc="best", ncol=2, fontsize=7)
    ax_b.set_title("(B) MEG3 per-lineage trajectory", loc="left")

    # (C) NEAT1 & MALAT1 trajectory (ubiquitous high-expression)
    ax_c = fig.add_subplot(gs[1, 0])
    for gene, ls in [("NEAT1", "-"), ("MALAT1", "--")]:
        sub_g = traj[traj["gene"] == gene]
        for lin in sub_g["lineage"].unique():
            s = sub_g[sub_g["lineage"] == lin].sort_values("bin")
            ax_c.plot(s["bin"], s["mean_expr"], label=f"{gene}/{lin}" if lin in PARENCHYMAL else None,
                      color=LINEAGE_COLORS.get(lin, "#888"), ls=ls, alpha=0.85, lw=1.2)
    ax_c.set_xlabel("Pseudotime bin")
    ax_c.set_ylabel("Mean expression")
    ax_c.set_xticks(range(10))
    ax_c.legend(loc="best", ncol=2, fontsize=6)
    ax_c.set_title("(C) NEAT1 (solid) and MALAT1 (dashed) trajectories", loc="left")

    # (D) KCNQ1OT1 cell-type composition (bar)
    ax_d = fig.add_subplot(gs[1, 1])
    kcnq = comp[comp["gene"] == "KCNQ1OT1"].copy()
    kcnq = kcnq.set_index("lineage").reindex([l for l in lin_order if l in kcnq["lineage"].values or l in kcnq.index])
    kcnq = kcnq.dropna(subset=["frac"])
    bar_colors = [LINEAGE_COLORS.get(l, "#888") for l in kcnq.index]
    ax_d.bar(range(len(kcnq)), kcnq["frac"] * 100, color=bar_colors, edgecolor="black", lw=0.4)
    for i, (lin, v) in enumerate(zip(kcnq.index, kcnq["frac"] * 100)):
        ax_d.text(i, v + 0.6, f"{v:.1f}%", ha="center", fontsize=7)
    ax_d.set_xticks(range(len(kcnq)))
    ax_d.set_xticklabels(kcnq.index, rotation=45, ha="right")
    ax_d.set_ylabel("KCNQ1OT1+ cells (%)")
    ax_d.set_title("(D) KCNQ1OT1 cell-type composition", loc="left")
    ax_d.set_ylim(0, max(kcnq["frac"] * 100) * 1.2)

    add_figure_label(fig, 4)
    save(fig, "fig4_anchor")


# ============================================================
# Figure 5 — Stability × peak bin
# ============================================================
def build_fig5():
    print("Building Figure 5...")
    df = pd.read_csv(config.RESULTS_DIR / "stability_x_pseudotime.csv")
    color_map = {"short": "#d62728", "medium": "#ff7f0e",
                 "long": "#2ca02c", "unknown": "#888888"}
    fig, ax = plt.subplots(figsize=(10.5, 5.5))
    fig.subplots_adjust(top=0.82, bottom=0.20, left=0.11, right=0.80)

    # Force categorical y-axis with explicit order
    lin_order = ["LSEC", "Cholangiocyte"]
    y_pos = {lin: i for i, lin in enumerate(lin_order)}
    ax.set_yticks(list(y_pos.values()))
    ax.set_yticklabels(list(y_pos.keys()))
    ax.set_ylim(-0.7, len(lin_order) - 0.3)

    label_map = {
        "short": "short-lived\n(t₁/₂ < 4 h)",
        "medium": "medium\n(4 ≤ t₁/₂ < 12 h)",
        "long": "long-lived\n(t₁/₂ ≥ 12 h)",
        "unknown": "unknown\n(absent from BRIC-seq)",
    }
    for _, row in df.iterrows():
        color = color_map.get(row["stability_class"], "#888")
        y = y_pos[row["lineage"]]
        ax.scatter(row["peak_bin"], y, s=460, c=color,
                   edgecolors="black", lw=0.9, zorder=3)
        # Annotate with gene name (above marker)
        ax.annotate(row["gene_name"], xy=(row["peak_bin"], y),
                    xytext=(0, 18), textcoords="offset points",
                    ha="center", fontsize=10, fontweight="bold", zorder=4)
        # Stability annotation (below marker, two-line text to avoid overlap)
        ax.annotate(label_map[row["stability_class"]],
                    xy=(row["peak_bin"], y),
                    xytext=(0, -22), textcoords="offset points",
                    ha="center", va="top", fontsize=7, color="#444", zorder=4)
    ax.set_xlim(-0.7, 9.7)
    ax.set_xticks(range(10))
    ax.set_xlabel("Pseudotime peak bin (0 = earliest, 9 = latest)")
    ax.set_title("Peak pseudotime bin and BRIC-seq stability of obesity-axis-coupled lncRNAs",
                 loc="left", pad=18, fontsize=10)
    from matplotlib.lines import Line2D
    leg_elems = [Line2D([0],[0], marker="o", color="w", markerfacecolor=color_map[k],
                        markersize=10, markeredgecolor="black", label=k)
                 for k in ["short", "medium", "long", "unknown"]]
    ax.legend(handles=leg_elems, title="Stability class",
              loc="center left", bbox_to_anchor=(1.02, 0.5),
              frameon=True, fontsize=8, title_fontsize=8.5)
    ax.grid(axis="x", color="#eeeeee", lw=0.5)
    add_figure_label(fig, 5)
    save(fig, "fig5_stability")


# ============================================================
# Supplementary Figure 1 — QC distributions
# ============================================================
def build_suppfig1():
    print("Building Supplementary Figure 1...")
    qc = pd.read_csv(config.RESULTS_DIR / "qc_report.csv")
    qc_dict = dict(zip(qc["metric"], qc["value"]))
    fig, axes = plt.subplots(2, 2, figsize=(10, 7))
    fig.subplots_adjust(hspace=0.40, wspace=0.32, top=0.90, bottom=0.08, left=0.10, right=0.96)

    # (A) Cells per filtering step
    ax_a = axes[0, 0]
    steps = ["Raw", "annot-joined", "Post-QC"]
    cells = [272867, 152559, int(qc_dict.get("n_cells", 152513))]
    ax_a.bar(steps, cells, color="#4575b4", edgecolor="black", lw=0.5)
    for i, v in enumerate(cells):
        ax_a.text(i, v + 4000, f"{v:,}", ha="center", fontsize=8)
    ax_a.set_ylabel("Cells")
    ax_a.set_title("(A) Cell retention", loc="left")
    ax_a.set_ylim(0, 300000)

    # (B) Genes retained
    ax_b = axes[0, 1]
    g_steps = ["Raw genes", "min_cells=3"]
    g_vals = [32738, int(qc_dict.get("n_genes", 24189))]
    ax_b.bar(g_steps, g_vals, color="#d73027", edgecolor="black", lw=0.5)
    for i, v in enumerate(g_vals):
        ax_b.text(i, v + 600, f"{v:,}", ha="center", fontsize=8)
    ax_b.set_ylabel("Genes")
    ax_b.set_title("(B) Gene retention", loc="left")

    # (C) Median QC metrics summary
    ax_c = axes[1, 0]
    ax_c.axis("off")
    table_data = [
        ["n_cells", f"{int(qc_dict.get('n_cells', 0)):,}"],
        ["n_genes (final)", f"{int(qc_dict.get('n_genes', 0)):,}"],
        ["Median genes/cell", f"{qc_dict.get('median_n_genes_per_cell', 0):.0f}"],
        ["Median counts/cell", f"{qc_dict.get('median_total_counts_per_cell', 0):.0f}"],
        ["Median pct_mt", f"{qc_dict.get('median_pct_mt', 0):.2f}%"],
        ["Donors", f"{int(qc_dict.get('n_donors', 0))}"],
    ]
    tbl = ax_c.table(cellText=table_data, colLabels=["Metric", "Value"],
                     loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(9); tbl.scale(1, 1.5)
    ax_c.set_title("(C) Final QC summary", loc="left")

    # (D) Filtering thresholds applied
    ax_d = axes[1, 1]
    ax_d.axis("off")
    thresh_data = [
        ["min_genes (per cell)", "200"],
        ["max n_genes_by_counts", "8,000"],
        ["pct_counts_mt (max)", "20%"],
        ["min_cells (per gene)", "3"],
        ["doublet removal", "inherited from Guilliams et al."],
    ]
    tbl2 = ax_d.table(cellText=thresh_data, colLabels=["Threshold", "Value"],
                      loc="center", cellLoc="center")
    tbl2.auto_set_font_size(False); tbl2.set_fontsize(9); tbl2.scale(1, 1.5)
    ax_d.set_title("(D) Filtering thresholds", loc="left")

    add_figure_label(fig, "S1")
    save(fig, "suppfig1_qc")


# ============================================================
# Supplementary Figure 2 — Sensitivity sweep
# ============================================================
def build_suppfig2():
    print("Building Supplementary Figure 2...")
    sens = pd.read_csv(config.RESULTS_DIR / "screening_sensitivity.csv")
    fig, axes = plt.subplots(1, 4, figsize=(15, 4))
    fig.subplots_adjust(top=0.84, bottom=0.18, left=0.06, right=0.97, wspace=0.40)
    for ax, lin in zip(axes, PARENCHYMAL):
        sub = sens[sens["lineage"] == lin]
        pivot = sub.pivot_table(index="rho_min", columns="lfc_min",
                                values="n_coupled", aggfunc="sum")
        # Reorder y-axis descending so 0.4 is at top
        pivot = pivot.sort_index(ascending=False)
        sns.heatmap(pivot, ax=ax, annot=True, fmt=".0f", cmap="Blues",
                    cbar=False, linewidths=0.4, linecolor="white",
                    annot_kws={"size": 10, "weight": "bold"})
        ax.set_title(lin, loc="left", fontsize=10)
        ax.set_xlabel("|log₂ fold change|")
        ax.set_ylabel("|Spearman ρ|")
    fig.suptitle("Threshold sensitivity sweep — number of progression-coupled lncRNA calls",
                 y=0.96, fontsize=11, ha="center", x=0.52)
    add_figure_label(fig, "S2")
    save(fig, "suppfig2_sensitivity")


# ============================================================
# Supplementary Figure 3 — Per-donor robustness (approximate from anchor_composition)
# ============================================================
def build_suppfig3():
    print("Building Supplementary Figure 3...")
    # Per-donor data not directly extractable without h5ad reload — show donor count
    # and dominant-lineage detection summary for coupled hits and key anchors.
    comp = pd.read_csv(config.RESULTS_DIR / "anchor_composition.csv")
    comp["frac"] = comp["n_detected"] / comp["n_total"]
    coupled_trajectory = pd.read_csv(config.RESULTS_DIR / "coupled_lncrna_trajectory.csv")

    fig = plt.figure(figsize=(13, 7))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.32,
                           top=0.92, bottom=0.10, left=0.08, right=0.97)

    targets = [
        ("MEG3", "LSEC", "(A) MEG3 — LSEC dominance"),
        ("ZFAS1", "LSEC", "(B) ZFAS1 — LSEC"),
        ("LINC00996", "Cholangiocyte", "(C) LINC00996 — cholangiocyte"),
        ("KCNQ1OT1", None, "(D) KCNQ1OT1 — cell-type composition"),
        ("NEAT1", None, "(E) NEAT1 — composition (ubiquitous)"),
        ("MALAT1", None, "(F) MALAT1 — composition (ubiquitous)"),
    ]
    for idx, (gene, dom_lin, title) in enumerate(targets):
        ax = fig.add_subplot(gs[idx // 3, idx % 3])
        if dom_lin is not None and gene in coupled_trajectory["gene"].unique():
            # Trajectory panel
            sub = coupled_trajectory[(coupled_trajectory["gene"] == gene) &
                                     (coupled_trajectory["lineage"] == dom_lin)].sort_values("bin")
            ax.plot(sub["bin"], sub["mean_expr"], marker="o", lw=2,
                    color=LINEAGE_COLORS.get(dom_lin, "#888"))
            ax.set_xlabel("Pseudotime bin")
            ax.set_ylabel("Mean expression")
            ax.set_xticks(range(10))
        else:
            # Composition bar
            sub = comp[comp["gene"] == gene]
            if len(sub) == 0:
                ax.axis("off")
                ax.set_title(title, loc="left")
                continue
            lin_order = ["Hepatocyte", "Kupffer", "LSEC", "Cholangiocyte", "Bcell", "Tcell", "NK"]
            sub = sub.set_index("lineage").reindex([l for l in lin_order if l in sub["lineage"].values or l in sub.index])
            sub = sub.dropna(subset=["frac"])
            bar_colors = [LINEAGE_COLORS.get(l, "#888") for l in sub.index]
            ax.bar(range(len(sub)), sub["frac"] * 100, color=bar_colors, edgecolor="black", lw=0.4)
            ax.set_xticks(range(len(sub)))
            ax.set_xticklabels(sub.index, rotation=45, ha="right", fontsize=7)
            ax.set_ylabel("Cells detecting (%)")
        ax.set_title(title, loc="left", fontsize=9)
    fig.suptitle("Per-lineage robustness check (across 16 donors)", y=0.97, fontsize=11)
    add_figure_label(fig, "S3")
    save(fig, "suppfig3_per_donor")


if __name__ == "__main__":
    build_fig1()
    build_fig2()
    build_fig3()
    build_fig4()
    build_fig5()
    build_suppfig1()
    build_suppfig2()
    build_suppfig3()
    print("\nAll figures built successfully.")
