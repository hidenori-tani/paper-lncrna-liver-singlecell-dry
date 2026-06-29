"""Fig 5: BRIC-seq stability class × pseudotime peak timing."""
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from pipeline import config

df = pd.read_csv(config.RESULTS_DIR / "stability_x_pseudotime.csv")
df = df.dropna(subset=["peak_bin"])

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

# (A) Boxplot
sns.boxplot(data=df, x="stability_class", y="peak_bin",
            order=["short", "medium", "long", "unknown"], ax=axes[0])
axes[0].set_title("(A) Peak timing by stability class")
axes[0].set_xlabel("BRIC-seq stability class")
axes[0].set_ylabel("Pseudotime peak bin")

# (B) Wave assignment
df["wave"] = pd.cut(df["peak_bin"], bins=[-0.1, 3, 6, 10], labels=["early", "mid", "late"])
cnt = df.groupby(["stability_class", "wave"]).size().unstack(fill_value=0)
cnt = cnt.reindex(["short", "medium", "long", "unknown"])
cnt.plot(kind="bar", stacked=True, ax=axes[1], colormap="coolwarm")
axes[1].set_title("(B) Wave assignment by stability class")
axes[1].set_ylabel("# lncRNA")

# (C) NEAT1 (short) vs MEG3 (long) Hepatocyte trajectory
traj = pd.read_csv(config.RESULTS_DIR / "anchor_trajectory.csv")
for target, color in [("NEAT1", "tab:red"), ("MEG3", "tab:blue")]:
    sub = traj[traj["gene"].str.upper().str.contains(target, na=False)]
    if sub.empty:
        continue
    hep = sub[sub["lineage"] == "Hepatocyte"].sort_values("bin")
    if not hep.empty:
        axes[2].plot(hep["bin"], hep["mean_expr"], label=target, marker="o", color=color)
axes[2].legend()
axes[2].set_title("(C) NEAT1 (short) vs MEG3 (long), Hepatocyte")
axes[2].set_xlabel("Pseudotime bin")
axes[2].set_ylabel("Mean expression")

fig.tight_layout()
fig.savefig(config.FIGURES_DIR / "fig5_stability.pdf", bbox_inches="tight")
fig.savefig(config.FIGURES_DIR / "fig5_stability.png", dpi=300, bbox_inches="tight")
print("Saved fig5_stability.{pdf,png}")
