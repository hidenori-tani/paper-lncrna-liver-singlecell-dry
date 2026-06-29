"""Fig 1: Atlas overview — pipeline schematic + UMAP + composition + lncRNA detection."""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

from pipeline import config, lncrna_utils

a = sc.read_h5ad(config.RESULTS_DIR / "liver_anchor_atlas.h5ad")

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# (A) Schematic placeholder
axes[0, 0].axis("off")
axes[0, 0].text(0.5, 0.5, "(A) Pipeline schematic\n(separately drawn)",
                ha="center", va="center", fontsize=14)

# (B) UMAP by lineage
sc.pl.umap(a, color="lineage", ax=axes[0, 1], show=False, frameon=False)
axes[0, 1].set_title("(B) Lineage UMAP")

# (C) Composition by condition
cond_key = next((k for k in ("condition", "disease", "diagnosis", "health_status") if k in a.obs.columns), None)
if cond_key:
    comp = a.obs.groupby([cond_key, "lineage"]).size().unstack(fill_value=0)
    comp.div(comp.sum(axis=1), axis=0).plot(kind="bar", stacked=True, ax=axes[1, 0])
    axes[1, 0].set_title(f"(C) Composition by {cond_key}")
    axes[1, 0].legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
else:
    axes[1, 0].axis("off")
    axes[1, 0].text(0.5, 0.5, "(C) no condition metadata", ha="center")

# (D) lncRNA detection per cell
a_lnc = lncrna_utils.filter_anndata_to_lncrna(a)
det = (a_lnc.X > 0).sum(axis=1)
det_arr = np.asarray(det).ravel()
df = pd.DataFrame({"n_lncrna_detected": det_arr, "lineage": a.obs["lineage"].values})
df.boxplot(column="n_lncrna_detected", by="lineage", ax=axes[1, 1])
axes[1, 1].set_title("(D) lncRNA detection per cell")
axes[1, 1].set_ylabel("# lncRNA detected")

fig.tight_layout()
config.FIGURES_DIR.mkdir(parents=True, exist_ok=True)
fig.savefig(config.FIGURES_DIR / "fig1_atlas.pdf", bbox_inches="tight")
fig.savefig(config.FIGURES_DIR / "fig1_atlas.png", dpi=300, bbox_inches="tight")
print("Saved fig1_atlas.{pdf,png}")
