"""Figure 1 を改訂の数値・Springer 規格で作り直す。

なぜ作り直すか（2026-08-19 実測）
--------------------------------
引き継いだ v1 の Figure_1_atlas.pdf は3つの点で出せない：
 1. 数値が本文と食い違う。図は投稿版の系統内訳（LSEC 1,191／Cholangiocyte 7,941
    ／Bcell 3,408 …）を表示するが、改訂本文は再実行の値（1,190／7,432／4,626 …）
    を報告し、「All numbers in this revision come from the re-execution」と書いている。
 2. 凡例と中身が違う。凡例は「(B) UMAP of 152,513 cells」と書くが、実物の panel B は
    棒グラフで、UMAP のパネルは無い。
 3. 誌の規定に反する。幅 277.8 mm・最小 4.4 pt・Type3 フォント。F&IG の規定は逐語で
    "size figures to fit in the column width ... 174 mm"、"lettering ... usually about
    2–3 mm (8–12 pt)"、"legible at final size"。174 mm に縮めると 2.8 pt になる。
"""
import pathlib, sys, os, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "lib"))
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from figstyle import use_journal_style, new_figure, savefig, OKABE_ITO as C

use_journal_style("springer")
RA = "revision_analysis"
OUT = "submission/funig_r1/figures"

lab = pd.read_csv(f"{RA}/global_labels.csv")[["col", "lineage_leiden", "diet"]]
frac = pd.read_csv(f"{RA}/anchor_detection_by_lineage_reexec.csv", index_col=0)

ORDER = ["Tcell", "NK", "Kupffer", "Cholangiocyte", "Bcell", "Hepatocyte", "LSEC"]
NICE = {"Tcell": "T cell", "NK": "NK", "Kupffer": "Kupffer",
        "Cholangiocyte": '"Cholangiocyte"', "Bcell": "B cell",
        "Hepatocyte": '"Hepatocyte"', "LSEC": "LSEC"}
n = lab.lineage_leiden.value_counts()
N_TOTAL = len(lab)

fig, axes = new_figure("springer", "double", height_mm=152, nrows=2, ncols=2)
fig.subplots_adjust(wspace=0.52, hspace=0.44, left=0.135, right=0.985,
                    bottom=0.115, top=0.945)

# ───────── A：解析の流れ ─────────
ax = axes[0, 0]
ax.axis("off")
ax.set_xlim(0, 1); ax.set_ylim(0, 1)
STAGES = ["Stage 0   10X Cell Ranger aggr matrix\n272,867 cells x 32,738 genes",
          f"Stage 1   Annotation join\n{N_TOTAL:,} cells x 24,190 genes",
          "Stage 2   lncRNA-aware HVG + Harmony\nLeiden + marker-score lineage",
          "Stage 3   Within-lineage diffusion\npseudotime (Lean-side root)",
          "Stage 4   Screening + anchor deep-dive\n+ BRIC-seq stability"]
h, gap = 0.165, 0.028
for i, body in enumerate(STAGES):
    y = 0.96 - i * (h + gap) - h
    ax.add_patch(FancyBboxPatch((0.01, y), 0.98, h, boxstyle="round,pad=0.004",
                                linewidth=0.6, edgecolor="#666666", facecolor="#F0F0F0"))
    ax.text(0.05, y + h / 2, body, fontsize=8, va="center", linespacing=1.30)
    if i < len(STAGES) - 1:
        ax.annotate("", xy=(0.5, y - 0.004), xytext=(0.5, y - gap + 0.004),
                    arrowprops=dict(arrowstyle="-|>", linewidth=0.6, color="#666666"))
ax.set_title("A  Analysis pipeline", loc="left")

# ───────── B：系統の内訳 ─────────
ax = axes[0, 1]
vals = [n[k] for k in ORDER]
ypos = np.arange(len(ORDER))[::-1]
ax.barh(ypos, vals, color=C["blue"], height=0.68)
for y_, v in zip(ypos, vals):
    ax.text(v + N_TOTAL * 0.012, y_, f"{v:,}", va="center", fontsize=8)
ax.set_yticks(ypos); ax.set_yticklabels([NICE[k] for k in ORDER])
ax.set_xlabel("Cells")
ax.set_xlim(0, max(vals) * 1.30)
ax.set_title(f"B  Marker-score lineages ({N_TOTAL:,} cells)", loc="left")

# ───────── C：条件の内訳 ─────────
ax = axes[1, 0]
ct = pd.crosstab(lab.lineage_leiden, lab.diet, normalize="index").reindex(ORDER)
x = np.arange(len(ORDER))
ax.bar(x, ct["Lean"], color=C["blue"], width=0.66, label="Lean")
ax.bar(x, ct["Obese"], bottom=ct["Lean"], color=C["vermillion"], width=0.66, label="Obese")
ax.set_xticks(x); ax.set_xticklabels([NICE[k] for k in ORDER], rotation=38, ha="right")
ax.set_ylabel("Fraction of cells")
ax.set_ylim(0, 1.34)
ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
ax.legend(loc="upper center", ncol=2, frameon=False)
ax.set_title("C  Composition", loc="left")

# ───────── D：アンカー lncRNA の検出 ─────────
ax = axes[1, 1]
s = frac[ORDER].sum()
ax.bar(x, [s[k] for k in ORDER],
       color=[C["vermillion"] if k == "LSEC" else C["blue"] for k in ORDER], width=0.66)
for i, k in enumerate(ORDER):
    ax.text(i, s[k] + 0.11, f"{s[k]:.2f}", ha="center", fontsize=8)
ax.set_xticks(x); ax.set_xticklabels([NICE[k] for k in ORDER], rotation=38, ha="right")
ax.set_ylabel("Sum of detection fractions\n(8 anchor lncRNAs)", linespacing=1.7)
ax.set_ylim(0, 3.5)
ax.set_title("D  Anchor lncRNA detection", loc="left")

savefig(fig, f"{OUT}/Figure_1_atlas")
plt.close(fig)
print("系統内訳:", {k: int(n[k]) for k in ORDER})
print("アンカー合計:", {k: round(float(s[k]), 3) for k in ORDER})
