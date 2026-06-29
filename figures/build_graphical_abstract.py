"""Build the Graphical Abstract for the manuscript.

Figure formatting guidelines:
- Maximum width 11.7 cm, height 9.9 cm (must remain legible at this size)
- Minimum resolution 300 DPI
- Preferred file types: EPS, TIFF
- Single, concise visual summary of the main finding
- New figure (NOT reproduced from the manuscript)

This abstract communicates the central finding:
  Bulk RNA-seq cannot localize the cellular origin of liver disease-associated
  lncRNAs (KCNQ1OT1, MEG3). Single-cell mapping re-localizes both signals
  predominantly to LSECs, not hepatocytes, refining the mechanistic interpretation.
"""
import sys
sys.path.insert(0, '.')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Wedge, Rectangle
import numpy as np

from pipeline import config

# Figure size: 11.7 cm × 9.9 cm at 600 DPI (well above minimum 300 DPI)
# 11.7 cm = 4.606 inches, 9.9 cm = 3.898 inches
FIG_W = 11.7 / 2.54  # inches
FIG_H = 9.9 / 2.54

LSEC_COLOR = "#2ca02c"
HEPATO_COLOR = "#d62728"
BCELL_COLOR = "#e377c2"
OTHER_COLOR = "#bbbbbb"
ARROW_COLOR = "#444"

fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=600)
ax = fig.add_subplot(111)
ax.set_xlim(0, 100)
ax.set_ylim(0, 100)
ax.axis("off")

# ============================================================
# Title
# ============================================================
ax.text(50, 96,
        "Single-cell mapping re-localizes liver disease-associated\n"
        "lncRNAs KCNQ1OT1 and MEG3 to sinusoidal endothelial cells",
        ha="center", va="top", fontsize=8.5, fontweight="bold",
        family="DejaVu Sans")

# ============================================================
# LEFT: Bulk RNA-seq view
# ============================================================
# Tissue box
left_box = FancyBboxPatch((2, 30), 28, 38,
                          boxstyle="round,pad=0.5,rounding_size=2",
                          linewidth=0.8, facecolor="#fff3e0", edgecolor="#cc7700")
ax.add_patch(left_box)
ax.text(16, 64, "Bulk RNA-seq", ha="center", va="top",
        fontsize=7.5, fontweight="bold", color="#9c4f00")
ax.text(16, 60, "(whole-tissue)", ha="center", va="top",
        fontsize=6, color="#9c4f00", style="italic")

# Liver silhouette (simplified shape using polygon)
liver_xy = [(8, 48), (12, 51), (18, 52), (23, 51), (26, 48),
            (25, 43), (22, 40), (17, 39), (12, 40), (9, 43)]
liver_poly = plt.Polygon(liver_xy, facecolor="#c39bd3", edgecolor="#6a3d9a", linewidth=0.8)
ax.add_patch(liver_poly)
ax.text(16, 45.5, "Liver", ha="center", va="center", fontsize=6.5, fontweight="bold",
        color="#3a1f5d")

ax.text(16, 34, "↑ KCNQ1OT1, MEG3", ha="center", va="top",
        fontsize=6.5, fontweight="bold", color="#d62728")
ax.text(16, 31, "in NAFLD", ha="center", va="top",
        fontsize=6, color="#666", style="italic")

ax.text(16, 25, "Which cell type?", ha="center", va="center",
        fontsize=7, color="#444", style="italic",
        bbox=dict(boxstyle="round,pad=0.3", fc="#ffffff", ec="#888", lw=0.5))

# ============================================================
# Arrow
# ============================================================
arrow = FancyArrowPatch((31, 49), (42, 49),
                        arrowstyle="-|>", mutation_scale=14,
                        linewidth=2, color=ARROW_COLOR)
ax.add_patch(arrow)
ax.text(36.5, 54, "Single-cell\nre-analysis", ha="center", va="bottom",
        fontsize=6.5, color=ARROW_COLOR, fontweight="bold")
ax.text(36.5, 44, "152,513 cells\n7 lineages", ha="center", va="top",
        fontsize=5.5, color="#666", style="italic")

# ============================================================
# RIGHT: Cell-type distribution
# ============================================================
right_box = FancyBboxPatch((44, 8), 54, 64,
                           boxstyle="round,pad=0.5,rounding_size=2",
                           linewidth=0.8, facecolor="#e8f5e9", edgecolor="#2e7d32")
ax.add_patch(right_box)
ax.text(71, 68, "Single-cell view", ha="center", va="top",
        fontsize=7.5, fontweight="bold", color="#1b5e20")

# Two stacked bars: KCNQ1OT1 and MEG3
def stacked_bar(ax, y_base, height, gene_label, segments, lineage_pct):
    """
    segments: list of (label, fraction, color, show_pct_threshold)
    """
    ax.text(47, y_base + height/2, gene_label, ha="left", va="center",
            fontsize=6.5, fontweight="bold", color="#222")
    x_start = 58
    bar_total_w = 38
    cum = 0
    for label, pct, color, show in segments:
        w = bar_total_w * pct / 100
        rect = Rectangle((x_start + cum, y_base), w, height,
                         facecolor=color, edgecolor="white", linewidth=0.4)
        ax.add_patch(rect)
        if show and pct >= 5:
            ax.text(x_start + cum + w/2, y_base + height/2,
                    f"{label}\n{pct:.1f}%",
                    ha="center", va="center", fontsize=4.5, color="white",
                    fontweight="bold")
        cum += w
    # Frame
    frame = Rectangle((x_start, y_base), bar_total_w, height,
                      facecolor="none", edgecolor="#666", linewidth=0.5)
    ax.add_patch(frame)

# KCNQ1OT1 row — total expressing cells distribution
# Among detecting cells: LSEC 260 + Bcell 447 + Tcell 2328 + NK 1594 + Kupffer 1953 + Cholangiocyte 567 + Hepatocyte 13 = 7,162
# But the key % is "% of cells detecting" PER LINEAGE — let me show that as a horizontal bar.
# Actually for clarity, let me show: of all KCNQ1OT1+ cells, what fraction come from each lineage.
# Total detected = 7,162. LSEC contrib = 260, B = 447, T = 2328, NK = 1594, Kupffer = 1953, Chol = 567, Hep = 13
# Hmm but the manuscript reports "21.8% of LSECs" which is detection rate.
# For the abstract, the clearer message: in each lineage, what % of cells detect the lncRNA?

# Use detection rates (% of cells in each lineage that detect the lncRNA)
# KCNQ1OT1: LSEC 21.8%, Bcell 13.1%, Cholangiocyte 7.1%, Kupffer 5.4%, NK 3.8%, Tcell 3.9%, Hepatocyte 0.4%
# Show as horizontal lineage-wise % bars

def lineage_bar(ax, y_base, height, gene_label, rates):
    """
    rates: list of (lineage, pct, color, fontcolor)
    """
    ax.text(47, y_base + height/2, gene_label, ha="left", va="center",
            fontsize=7, fontweight="bold", color="#222")
    x_start = 58
    bar_width_max = 38
    # Normalize: show absolute % up to ~80% scale
    scale = bar_width_max / 80
    cum = 0
    # Show only top 4 lineages
    for lineage, pct, color in rates[:4]:
        w = pct * scale
        rect = Rectangle((x_start + cum, y_base), w, height,
                         facecolor=color, edgecolor="white", linewidth=0.4)
        ax.add_patch(rect)
        if pct >= 8:
            txt_color = "white"
            ax.text(x_start + cum + w/2, y_base + height/2,
                    f"{lineage}\n{pct:.1f}%",
                    ha="center", va="center", fontsize=4.8, color=txt_color,
                    fontweight="bold")
        else:
            ax.text(x_start + cum + w/2, y_base + height + 0.5,
                    f"{lineage} {pct:.1f}%",
                    ha="center", va="bottom", fontsize=4.2, color="#333")
        cum += w
    # Hepatocyte annotation (very small)
    hep_pct = next((p for lin, p, _ in rates if lin == "Hep"), 0)
    if hep_pct > 0:
        ax.text(x_start + cum + 0.5, y_base + height/2,
                f"Hep {hep_pct:.1f}%",
                ha="left", va="center", fontsize=4.5, color="#d62728",
                fontweight="bold", style="italic")

# KCNQ1OT1 detection rates (% of cells in each lineage detecting)
kcnq_rates = [
    ("LSEC", 21.8, LSEC_COLOR),
    ("Bcell", 13.1, BCELL_COLOR),
    ("Chol", 7.1, "#ff7f0e"),
    ("Kupffer", 5.4, "#1f77b4"),
    ("Hep", 0.4, HEPATO_COLOR),
]
lineage_bar(ax, 50, 8, "KCNQ1OT1", kcnq_rates)
ax.text(47, 60, "(% of cells per lineage detecting)", ha="left", va="bottom",
        fontsize=5, color="#666", style="italic")

# MEG3 detection rates
meg3_rates = [
    ("LSEC", 62.0, LSEC_COLOR),
    ("Chol", 4.9, "#ff7f0e"),
    ("Kupffer", 0.4, "#1f77b4"),
    ("Bcell", 0.2, BCELL_COLOR),
    ("Hep", 0.1, HEPATO_COLOR),
]
lineage_bar(ax, 33, 8, "MEG3", meg3_rates)

# ============================================================
# Bottom take-home message
# ============================================================
bottom_box = FancyBboxPatch((44, 11), 54, 14,
                            boxstyle="round,pad=0.4,rounding_size=2",
                            linewidth=0.8, facecolor="#e3f2fd", edgecolor="#1565c0")
ax.add_patch(bottom_box)
ax.text(71, 21, "Key finding", ha="center", va="top",
        fontsize=6.5, fontweight="bold", color="#0d47a1")
ax.text(71, 17, "LSECs — not hepatocytes —", ha="center", va="top",
        fontsize=6.5, fontweight="bold", color=LSEC_COLOR)
ax.text(71, 13.5, "host the bulk-level liver lncRNA signal",
        ha="center", va="top", fontsize=6, color="#0d47a1")

# Subtle footer with dataset attribution
ax.text(50, 3.5, "Guilliams 2022 Liver Cell Atlas re-analysis (GSE192742, 16 donors, Lean-to-Obese gradient)",
        ha="center", va="bottom", fontsize=5, color="#888", style="italic")

# ============================================================
# Save
# ============================================================
out_dir = config.FIGURES_DIR
out_dir.mkdir(parents=True, exist_ok=True)

base = out_dir / "graphical_abstract"
fig.savefig(str(base) + ".tif", dpi=600, bbox_inches="tight",
            pil_kwargs={"compression": "tiff_lzw"})
fig.savefig(str(base) + ".png", dpi=600, bbox_inches="tight")
fig.savefig(str(base) + ".pdf", bbox_inches="tight")
print(f"Saved: {base}.{{tif,png,pdf}}")

# Verify dimensions
from PIL import Image
img = Image.open(str(base) + ".png")
print(f"PNG size: {img.size[0]} × {img.size[1]} pixels (target: ≥442 × 374, actual ratio {img.size[0]/img.size[1]:.3f})")
print(f"Effective DPI: ~{img.size[0]/FIG_W:.0f}")
