"""補足図 S1・S2 を改訂の数値で作り直す（Springer 寸法・figstyle 経由）。

なぜ作り直すか（2026-08-19）
---------------------------
v1 の `Supplementary_figures.pdf` は、改訂本文が撤回した用語と系統名を図の中に
刷り込んでいる：
  - S2 の図中の題が "number of progression-coupled lncRNA calls"
    → `progression-coupled` は改訂本文から 0 件になった語。
  - S3 の panel C が "LINC00996 — cholangiocyte"
    → LINC00996 は撤回、その "cholangiocyte" 系統は胆管細胞を 1 個も含まない。
そのまま納品すると、撤回した主張が掲載物として残る。

Fig. S1  品質管理の要約（投稿版パイプラインの QC ＋ 改訂の再実行との差）
Fig. S2  スクリーニング閾値の感度掃引（再実行の数値・公表アノテ側も併置）
Fig. S3  作らない（撤回した3件が主題だったため。本文の引用も外す）
"""
import pathlib, sys, os, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "lib"))
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from figstyle import use_journal_style, new_figure, savefig, OKABE_ITO as C

use_journal_style("springer")
RA = "revision_analysis"
OUT = "submission/funig_r1/supplementary"
os.makedirs(OUT, exist_ok=True)

# ────────────────────── Fig. S1：品質管理の要約 ──────────────────────
qc = pd.read_csv("submission/funig_v1/supplementary/qc_report.csv").set_index("metric")["value"]
N_RAW, N_JOIN, N_QC = 272_867, 152_559, int(qc["n_cells"])
G_RAW, G_QC = 32_738, int(qc["n_genes"])
G_REEX = len(pd.read_csv(f"{RA}/gene_level_stats.csv").query("n_detected >= 3"))

fig, axes = new_figure("springer", "double", height_mm=74, nrows=1, ncols=2)
fig.subplots_adjust(wspace=0.30, bottom=0.20, top=0.86, left=0.10, right=0.985)

ax = axes[0]
vals = [N_RAW, N_JOIN, N_QC]
ax.bar(range(3), vals, color=[C["blue"], C["sky_blue"], C["vermillion"]], width=0.60)
for i, v in enumerate(vals):
    ax.text(i, v * 1.045, f"{v:,}", ha="center", va="bottom", fontsize=8)
ax.set_xticks(range(3))
ax.set_xticklabels(["Raw", "Annotation-joined\n(this revision)", "Post-QC\n(submitted)"])
ax.set_ylabel("Cells")
ax.set_ylim(0, N_RAW * 1.26)
ax.set_xlim(-0.72, 2.72)
ax.set_title("A  Cell retention", loc="left")

ax = axes[1]
vals = [G_RAW, G_QC, G_REEX]
ax.bar(range(3), vals, color=[C["blue"], C["vermillion"], C["sky_blue"]], width=0.60)
for i, v in enumerate(vals):
    ax.text(i, v * 1.045, f"{v:,}", ha="center", va="bottom", fontsize=8)
ax.set_xticks(range(3))
ax.set_xticklabels(["Raw", "min_cells = 3\n(submitted)", "min_cells = 3\n(this revision)"])
ax.set_ylabel("Genes")
ax.set_ylim(0, G_RAW * 1.26)
ax.set_xlim(-0.72, 2.72)
ax.set_title("B  Gene retention", loc="left")

savefig(fig, f"{OUT}/Figure_S1_quality_control")
plt.close(fig)
print(f"Fig. S1: cells {N_RAW:,}→{N_JOIN:,}→{N_QC:,} / genes {G_RAW:,}→{G_QC:,} (submitted) / {G_REEX:,} (revision)")

# ────────────────────── Fig. S2：閾値の感度掃引 ──────────────────────
d = pd.read_csv(f"{RA}/rescreen_global_table.csv")
RHOS = [0.2, 0.3, 0.4]
LFCS = [0.5, 1.0, 1.5]
PANELS = [("lineage_leiden", "Hepatocyte"), ("lineage_leiden", "Kupffer"),
          ("lineage_leiden", "LSEC"), ("lineage_leiden", "Cholangiocyte"),
          ("lineage_pub", "Kupffer"), ("lineage_pub", "LSEC")]
TITLE = {"lineage_leiden": "Marker-score", "lineage_pub": "Published"}

fig, axes = new_figure("springer", "double", height_mm=108, nrows=2, ncols=3)
fig.subplots_adjust(hspace=0.75, wspace=0.45, bottom=0.13, top=0.93, left=0.12, right=0.985)
axes = np.atleast_2d(axes)
vmax = 0
grids = {}
for ls, ln in PANELS:
    sub = d[(d.label_set == ls) & (d.lineage_name == ln)]
    g = np.zeros((len(LFCS), len(RHOS)), int)
    for i, l in enumerate(LFCS):
        for j, r in enumerate(RHOS):
            g[i, j] = int(((sub.spearman_rho.abs() >= r) & (sub.wald_q < 0.05)
                           & (sub.log2fc_max_min.abs() >= l)).sum())
    grids[(ls, ln)] = g
    vmax = max(vmax, g.max())

for k, (ls, ln) in enumerate(PANELS):
    ax = axes[k // 3, k % 3]
    g = grids[(ls, ln)]
    ax.imshow(g, cmap="Blues", vmin=0, vmax=max(vmax, 1), aspect="auto")
    for i in range(len(LFCS)):
        for j in range(len(RHOS)):
            ax.text(j, i, str(g[i, j]), ha="center", va="center", fontsize=8,
                    color="white" if g[i, j] > max(vmax, 1) * 0.55 else "black")
    ax.set_xticks(range(len(RHOS))); ax.set_xticklabels([f"{r:.1f}" for r in RHOS])
    ax.set_yticks(range(len(LFCS))); ax.set_yticklabels([f"{l:.1f}" for l in LFCS])
    ax.set_title(f"{'ABCDEF'[k]}  {TITLE[ls]}: {ln}", loc="left")
    if k // 3 == 1:
        ax.set_xlabel("|Spearman rho| threshold")
    if k % 3 == 0:
        ax.set_ylabel("|log2 range|\nthreshold", linespacing=1.6)
savefig(fig, f"{OUT}/Figure_S2_threshold_sensitivity")
plt.close(fig)
for ls, ln in PANELS:
    print(f"Fig. S2  {ls}/{ln}: 規定の閾値(0.3,1.0)で {grids[(ls,ln)][1,1]} 件")
