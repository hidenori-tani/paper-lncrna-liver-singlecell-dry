"""Fig 2: Unbiased progression-coupled lncRNA screening summary."""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

try:
    from matplotlib_venn import venn2, venn3
    HAS_VENN = True
except Exception:
    HAS_VENN = False

from pipeline import config

df = pd.read_csv(config.RESULTS_DIR / "progression_lncrna_table.csv")
df_c = df[df["is_progression_coupled"]]

fig, axes = plt.subplots(2, 2, figsize=(13, 11))

# (A) per-lineage UMAP — refer to Fig 3
axes[0, 0].axis("off")
axes[0, 0].text(0.5, 0.5, "(A) per-lineage UMAP\n(see Fig 3)", ha="center", va="center")

# (B) Heatmap: top 50 progression-coupled lncRNA × bin
top_n = 50
top_genes = df_c.nlargest(top_n, "spearman_rho")["gene"].tolist()
try:
    traj = pd.read_csv(config.RESULTS_DIR / "anchor_trajectory.csv")
    pivot = traj[traj["gene"].isin(top_genes)].pivot_table(index="gene", columns="bin", values="mean_expr", aggfunc="mean")
    sns.heatmap(pivot, ax=axes[0, 1], cmap="viridis")
    axes[0, 1].set_title(f"(B) Top {top_n} progression-coupled lncRNA")
except Exception:
    axes[0, 1].text(0.5, 0.5, "(B) trajectory CSV missing", ha="center")

# (C) Counts per lineage
counts = df_c.groupby("lineage").size().sort_values(ascending=False)
counts.plot(kind="bar", ax=axes[1, 0], color="steelblue")
axes[1, 0].set_title("(C) Progression-coupled lncRNA count per lineage")
axes[1, 0].set_ylabel("# lncRNA")

# (D) Venn diagram for top 3 lineages
top_lin = counts.head(3).index.tolist()
sets = {lin: set(df_c[df_c["lineage"] == lin]["gene"]) for lin in top_lin}
if HAS_VENN and len(sets) == 3:
    venn3(sets.values(), set_labels=tuple(top_lin), ax=axes[1, 1])
elif HAS_VENN and len(sets) == 2:
    venn2(sets.values(), set_labels=tuple(top_lin), ax=axes[1, 1])
else:
    axes[1, 1].text(0.5, 0.5, "(D) install matplotlib-venn\nfor lineage overlap", ha="center", va="center")
axes[1, 1].set_title("(D) Lineage overlap")

fig.tight_layout()
fig.savefig(config.FIGURES_DIR / "fig2_screening.pdf", bbox_inches="tight")
fig.savefig(config.FIGURES_DIR / "fig2_screening.png", dpi=300, bbox_inches="tight")
print("Saved fig2_screening.{pdf,png}")
