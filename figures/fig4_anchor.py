"""Fig 4: Anchor lncRNA deep-dive (Kcnq1ot1 / NEAT1 / MALAT1 / MEG3 / IDI2-AS1)."""
import matplotlib.pyplot as plt
import pandas as pd

from pipeline import config

traj = pd.read_csv(config.RESULTS_DIR / "anchor_trajectory.csv")
comp = pd.read_csv(config.RESULTS_DIR / "anchor_composition.csv")
coex = pd.read_csv(config.RESULTS_DIR / "anchor_coexpression.csv")

fig, axes = plt.subplots(2, 2, figsize=(13, 11))

# (A) KCNQ1OT1 per-lineage trajectory
target = "KCNQ1OT1"
sub = traj[traj["gene"].str.upper().str.contains(target, na=False)]
for lin in sub["lineage"].unique():
    s = sub[sub["lineage"] == lin].sort_values("bin")
    axes[0, 0].plot(s["bin"], s["mean_expr"], label=lin, marker="o")
axes[0, 0].legend()
axes[0, 0].set_xlabel("Pseudotime bin")
axes[0, 0].set_ylabel("Mean expression")
axes[0, 0].set_title(f"(A) {target} trajectory")

# (B) NEAT1 / MALAT1 / MEG3 detection composition
canonical = ["NEAT1", "MALAT1", "MEG3"]
canon_comp = comp[comp["gene"].str.upper().isin(canonical)].copy()
if not canon_comp.empty:
    canon_comp["fraction_detected"] = canon_comp["n_detected"] / canon_comp["n_total"]
    pivot = canon_comp.pivot_table(index="lineage", columns="gene", values="fraction_detected", aggfunc="mean")
    pivot.plot(kind="bar", ax=axes[0, 1])
    axes[0, 1].set_title("(B) Canonical lncRNA detection")
    axes[0, 1].set_ylabel("Fraction of cells detecting")

# (C) IDI2-AS1 hepatic mapping
target3 = "IDI2-AS1"
sub3 = traj[traj["gene"].str.upper().str.contains("IDI2", na=False)]
if not sub3.empty:
    for lin in sub3["lineage"].unique():
        s = sub3[sub3["lineage"] == lin].sort_values("bin")
        axes[1, 0].plot(s["bin"], s["mean_expr"], label=lin, marker="s")
    axes[1, 0].legend()
    axes[1, 0].set_title(f"(C) {target3} hepatic mapping")

# (D) Coexpression module summary — top 5 KCNQ1OT1 partners in Hepatocyte
sub4 = coex[(coex["gene"].str.upper().str.contains("KCNQ1OT1", na=False)) & (coex["lineage"] == "Hepatocyte")]
sub4 = sub4.nlargest(5, "rho")
axes[1, 1].barh(sub4["partner"], sub4["rho"], color="darkgreen")
axes[1, 1].set_title("(D) KCNQ1OT1 top coexpressed (Hepatocyte)")
axes[1, 1].set_xlabel("Spearman rho")

fig.tight_layout()
fig.savefig(config.FIGURES_DIR / "fig4_anchor.pdf", bbox_inches="tight")
fig.savefig(config.FIGURES_DIR / "fig4_anchor.png", dpi=300, bbox_inches="tight")
print("Saved fig4_anchor.{pdf,png}")
