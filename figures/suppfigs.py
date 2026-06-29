"""Supplementary figures: QC + sensitivity + per-donor variability."""
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc

from pipeline import config

# Supp Fig 1: QC report
qc = pd.read_csv(config.RESULTS_DIR / "qc_report.csv")
fig, ax = plt.subplots(figsize=(8, 4))
ax.axis("off")
ax.text(0, 0.5, qc.to_string(index=False), family="monospace", fontsize=10)
fig.savefig(config.FIGURES_DIR / "suppfig1_qc.pdf", bbox_inches="tight")
print("Saved suppfig1_qc.pdf")

# Supp Fig 2: Sensitivity sweep heatmap
try:
    import seaborn as sns
    sens = pd.read_csv(config.RESULTS_DIR / "screening_sensitivity.csv")
    pivot = sens.pivot_table(index="rho_min", columns="lfc_min", values="n_coupled", aggfunc="sum")
    fig, ax = plt.subplots(figsize=(7, 5))
    sns.heatmap(pivot, annot=True, fmt=".0f", ax=ax, cmap="Blues")
    ax.set_title("Sensitivity: # progression-coupled lncRNA")
    fig.savefig(config.FIGURES_DIR / "suppfig2_sensitivity.pdf", bbox_inches="tight")
    print("Saved suppfig2_sensitivity.pdf")
except FileNotFoundError:
    print("screening_sensitivity.csv missing, skipping suppfig2")

# Supp Fig 3: Per-donor variability of anchor lncRNAs
a = sc.read_h5ad(config.RESULTS_DIR / "liver_anchor_pseudotime.h5ad")
donor_key = next((k for k in ("donor_id", "sample", "patient") if k in a.obs.columns), None)
if donor_key:
    fig, axes = plt.subplots(2, 3, figsize=(15, 9))
    for ax, gene in zip(axes.ravel(), ["KCNQ1OT1", "RMST", "NEAT1", "MALAT1", "MEG3", "IDI2-AS1"]):
        if gene not in a.var_names:
            ax.axis("off")
            ax.set_title(f"{gene} (not in atlas)")
            continue
        gi = a.var_names.get_loc(gene)
        X = a.X[:, gi].toarray().ravel() if hasattr(a.X, "toarray") else a.X[:, gi]
        df = pd.DataFrame({"expr": X, "donor": a.obs[donor_key].values})
        df.boxplot(column="expr", by="donor", ax=ax, rot=90)
        ax.set_title(gene)
    fig.tight_layout()
    fig.savefig(config.FIGURES_DIR / "suppfig3_per_donor.pdf", bbox_inches="tight")
    print("Saved suppfig3_per_donor.pdf")
