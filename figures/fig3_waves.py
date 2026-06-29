"""Fig 3: Wave structure — peak timing dotplot, k-means clusters of trajectory."""
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.cluster import KMeans

from pipeline import config

stab = pd.read_csv(config.RESULTS_DIR / "stability_x_pseudotime.csv")
traj = pd.read_csv(config.RESULTS_DIR / "anchor_trajectory.csv")

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# (A) Peak timing dotplot per lineage
for i, lin in enumerate(config.PSEUDOTIME["lineages"][:5]):
    sub = stab[stab["lineage"] == lin].sort_values("peak_bin")
    axes[0].scatter(sub["peak_bin"], [i] * len(sub), s=10, alpha=0.5, label=lin if i < 5 else None)
axes[0].set_yticks(range(len(config.PSEUDOTIME["lineages"])))
axes[0].set_yticklabels(config.PSEUDOTIME["lineages"])
axes[0].set_xlabel("Pseudotime peak bin")
axes[0].set_title("(A) Peak timing by lineage")

# (B) K-means on lncRNA trajectories
if not traj.empty:
    mat = traj.pivot_table(index="gene", columns="bin", values="mean_expr", aggfunc="mean").fillna(0)
    mat_norm = mat.div(mat.max(axis=1).replace(0, 1), axis=0)
    n_clusters = 4
    km = KMeans(n_clusters=n_clusters, random_state=config.RANDOM_SEED, n_init=10).fit(mat_norm.values)
    labels = km.labels_
    for c in range(n_clusters):
        mean_curve = mat_norm.values[labels == c].mean(axis=0)
        axes[1].plot(mean_curve, label=f"Cluster {c} (n={int((labels == c).sum())})")
    axes[1].legend()
    axes[1].set_xlabel("Pseudotime bin")
    axes[1].set_ylabel("Normalized mean expression")
    axes[1].set_title("(B) Wave clusters")

# (C) Cell-type sequence: median peak bin per lineage
seq = stab.groupby("lineage")["peak_bin"].median().sort_values()
axes[2].barh(seq.index, seq.values, color="coral")
axes[2].set_xlabel("Median peak bin")
axes[2].set_title("(C) Cell-type ordering")

fig.tight_layout()
fig.savefig(config.FIGURES_DIR / "fig3_waves.pdf", bbox_inches="tight")
fig.savefig(config.FIGURES_DIR / "fig3_waves.png", dpi=300, bbox_inches="tight")
print("Saved fig3_waves.{pdf,png}")
