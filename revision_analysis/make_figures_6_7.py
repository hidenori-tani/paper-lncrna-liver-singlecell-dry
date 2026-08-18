"""Figure 6・7 を Springer 規格（174 mm 幅・8 pt 以上・埋め込みフォント）で作り直す。

なぜ作り直すか（2026-08-19 実測）
--------------------------------
引き継いだ v1 の図は誌の規定に反する。F&IG の規定は逐語で
  "size figures to fit in the column width ... 174 mm"
  "lettering ... usually about 2-3 mm (8-12 pt)"
  "check that all lines and lettering within the figures are legible at final size"
実測は Fig6 355.6 mm / 最小 7.0 pt、Fig7 351.4 mm / 最小 7.0 pt。174 mm に縮めると
3.4-3.5 pt になり、規定を満たさない。加えて Type3 フォント（DejaVu）を埋め込んでいる。

中身は変えない。同じ一次データ・同じ量・同じ順序で、版面と字送りだけを規格に合わせる。
横長を縦方向の帯（細胞型を y 軸）に組み替えたのは、174 mm では細胞型名が
回転ラベルとして読めないため（規定の "legible at final size"）。
"""
import pathlib, sys, os, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "lib"))
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from figstyle import use_journal_style, new_figure, savefig, OKABE_ITO as C

use_journal_style("springer")
OUT = "submission/funig_r1/figures"
# Fig.6 は凡例が「cell types with >= 20 cells in a donor shown」と明記しているのでその規則。
# Fig.7 は凡例に規則が無く、凡例本文が「HSC は 2-24 cells per donor」と書いている＝
# 少数細胞のドナーも図に入っている。20 細胞で切ると 78 点中 28 点（36%）が消え、
# 凡例が言及している点そのものが図から無くなる。よって Fig.7 は全ドナーを出す。
MIN_CELLS = 20      # Fig.6
MIN_CELLS_7 = 1     # Fig.7（全ドナー）


def strip_panel(ax, df, gene, highlight, order, title, xlabel, min_cells=MIN_CELLS):
    sub = df[(df.gene == gene) & (df.n_cells >= min_cells)]
    ypos = {c: i for i, c in enumerate(order[::-1])}
    for c in order:
        v = sub.loc[sub.cell_type == c, "detection_rate"].values * 100
        if not len(v):
            continue
        y = ypos[c]
        col = C["vermillion"] if c in highlight else C["blue"]
        ax.scatter(v, np.full(len(v), y) + np.linspace(-0.16, 0.16, len(v)),
                   s=7, color=col, linewidths=0, zorder=3, clip_on=False)
        ax.plot([np.median(v)] * 2, [y - 0.34, y + 0.34], color="black",
                linewidth=1.1, zorder=4, solid_capstyle="butt")
    nd = (df[(df.gene == gene) & (df.n_cells >= min_cells)]
          .groupby("cell_type").donor.nunique())
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([f"{c.replace('_', ' ')} ({nd.get(c, 0)})" for c in order[::-1]])
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel(xlabel)
    ax.set_title(title, loc="left")
    ax.grid(axis="x", linewidth=0.3, color="#DDDDDD", zorder=0)
    ax.set_axisbelow(True)


def set_x(ax, df, gene, min_cells=MIN_CELLS):
    """データの最大値から x 範囲を決め、1点でも枠外なら落とす。

    🚨 2026-08-19、x 上限を手で 62 と書いたため KCNQ1OT1 の 74.2%（胆管細胞・
       ドナー1人分）が版面の外に落ちた。図ゲートも重なりゲートも検出しない
       （枠の中は正しく、外に出た点は「無い」ものとして扱われる）。
       上限はデータから決め、内蔵の検査で毎回確かめる。
    """
    v = df[(df.gene == gene) & (df.n_cells >= min_cells)].detection_rate.values * 100
    hi = float(np.ceil(v.max() / 10.0) * 10)
    ax.set_xlim(-hi * 0.03, hi * 1.03)
    ax.set_xticks(np.linspace(0, hi, 5) if hi % 4 == 0 else
                  np.arange(0, hi + 1, 20 if hi > 60 else 10))
    lo, up = ax.get_xlim()
    out = [x for x in v if not (lo <= x <= up)]
    assert not out, f"{gene}: {len(out)} 点が版面外（{out[:3]}） 上限 {up}"
    return hi


def order_by_median(df, gene, min_cells=MIN_CELLS):
    sub = df[(df.gene == gene) & (df.n_cells >= min_cells)]
    med = sub.groupby("cell_type").detection_rate.median().sort_values(ascending=False)
    return list(med.index)


# ────────────────────────── Figure 6 ──────────────────────────
d6 = pd.read_csv("submission/jgg_v1/donor_level/per_donor_detection.csv")
HL6 = {"Endothelial cells"}
order6 = order_by_median(d6, "MEG3")

fig, axes = new_figure("springer", "double", height_mm=104, nrows=1, ncols=2)
fig.subplots_adjust(wspace=0.10, left=0.275, right=0.99, bottom=0.13, top=0.90)
strip_panel(axes[0], d6, "KCNQ1OT1", HL6, order6,
            "A  KCNQ1OT1", "Per-donor detection (%)")
strip_panel(axes[1], d6, "MEG3", HL6, order6,
            "B  MEG3", "Per-donor detection (%)")
axes[1].set_yticklabels([])
set_x(axes[0], d6, "KCNQ1OT1"); set_x(axes[1], d6, "MEG3")
savefig(fig, f"{OUT}/Figure_6_donor_validation")
plt.close(fig)
s6 = d6[(d6.gene == "MEG3") & (d6.n_cells >= MIN_CELLS)]
print("Fig6 細胞型:", len(order6),
      "| MEG3 の LSEC 中央値: %.1f%%" %
      (s6.loc[s6.cell_type == "Endothelial cells", "detection_rate"].median() * 100),
      "| 次点: %.1f%%" %
      (s6[s6.cell_type != "Endothelial cells"].groupby("cell_type")
       .detection_rate.median().max() * 100))

# ────────────────────────── Figure 7 ──────────────────────────
d7 = pd.read_csv("replication_macparland/macparland_per_donor_detection.csv")
HL7 = {"Central_venous_LSECs", "Periportal_LSECs"}
order7 = order_by_median(d7, "MEG3", MIN_CELLS_7)

fig, axes = new_figure("springer", "double", height_mm=112, nrows=1, ncols=2)
fig.subplots_adjust(wspace=0.10, left=0.315, right=0.99, bottom=0.12, top=0.905)
strip_panel(axes[0], d7, "MEG3", HL7, order7,
            "A  MEG3", "Per-donor detection (%)", MIN_CELLS_7)
strip_panel(axes[1], d7, "KCNQ1OT1", HL7, order7,
            "B  KCNQ1OT1", "Per-donor detection (%)", MIN_CELLS_7)
axes[1].set_yticklabels([])
set_x(axes[0], d7, "MEG3", MIN_CELLS_7); set_x(axes[1], d7, "KCNQ1OT1", MIN_CELLS_7)
savefig(fig, f"{OUT}/Figure_7_replication")
plt.close(fig)
s7 = d7[(d7.gene == "MEG3") & (d7.n_cells >= MIN_CELLS_7)]
print("Fig7 細胞型:", len(order7),
      "| MEG3 の Central venous LSEC 中央値: %.1f%%" %
      (s7.loc[s7.cell_type == "Central_venous_LSECs", "detection_rate"].median() * 100))
