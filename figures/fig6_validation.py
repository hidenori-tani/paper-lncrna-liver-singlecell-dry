"""Fig 6: Cross-dataset validation — projection + replication scatter + table."""
import matplotlib.pyplot as plt
import pandas as pd

from pipeline import config

val = pd.read_csv(config.RESULTS_DIR / "validation_projection.csv")

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# (A) Projection placeholder
axes[0].axis("off")
axes[0].text(0.5, 0.5, "(A) UMAP overlay\n(separately rendered)", ha="center", va="center")

# (B) Replication scatter
axes[1].scatter(val["anchor_rho"], val["val_rho"], alpha=0.5, s=15)
lim = max(abs(val["anchor_rho"]).max(), abs(val["val_rho"]).max()) * 1.1
axes[1].plot([-lim, lim], [-lim, lim], "k--", alpha=0.4)
axes[1].axhline(0, color="gray", linewidth=0.5)
axes[1].axvline(0, color="gray", linewidth=0.5)
axes[1].set_xlabel("Anchor Spearman rho")
axes[1].set_ylabel("Validation Spearman rho")
axes[1].set_title("(B) Replication scatter")

# (C) Replication summary table
summary = val.groupby("lineage")["replicated"].agg(["sum", "count"])
summary["rate"] = summary["sum"] / summary["count"]
summary_str = summary.round(2).to_string()
axes[2].axis("off")
axes[2].text(0.05, 0.5, "(C) Replication rate per lineage\n\n" + summary_str,
             family="monospace", fontsize=10, verticalalignment="center")

fig.tight_layout()
fig.savefig(config.FIGURES_DIR / "fig6_validation.pdf", bbox_inches="tight")
fig.savefig(config.FIGURES_DIR / "fig6_validation.png", dpi=300, bbox_inches="tight")
print("Saved fig6_validation.{pdf,png}")
