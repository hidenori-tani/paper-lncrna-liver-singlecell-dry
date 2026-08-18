#!/usr/bin/env python3
"""figstyle.py — 投稿可能な図を matplotlib で最初から作るための共有スタイル

なぜ必要か：既存の図生成スクリプトが matplotlib の既定のまま書き出しており、
`check_figure.py` で全 PDF が下記に該当した（2026-07-29 実測・18点中9 PDF 全部）。

  1. **Type3 フォント**（matplotlib 既定の `pdf.fonttype=3`）
     → Type1/TrueType を要求する誌の自動チェックで弾かれる
  2. **DejaVu Sans**（matplotlib 既定書体）
     → Arial/Helvetica/Times/Courier/Symbol しか許さない誌が多い
  3. **掲載サイズで作っていない**（幅 229 mm など）
     → 90 mm カラムに縮小されると 10 pt が 3.9 pt になり規定違反

このモジュールを図生成スクリプトの冒頭で呼べば、3つとも構造的に起こらない。

使い方
    import sys; sys.path.insert(0, "<repo>/research/scripts")
    from figstyle import use_journal_style, mm, new_figure, savefig, OKABE_ITO

    use_journal_style("elsevier")             # rcParams を投稿仕様に固定
    fig, ax = new_figure("elsevier", "single", height_mm=60)
    ax.plot(x, y, color=OKABE_ITO["blue"], label="control")
    savefig(fig, "figures/fig1")              # fig1.pdf と fig1.png を同時に出す

検証：`python3 check_figure.py figures/fig1.pdf --journal elsevier --column single`
      が FAIL 0 になることを確認済み（下部 `_selftest` を実行すると再現する）。
"""

from __future__ import annotations

import os

import matplotlib
import matplotlib.pyplot as plt

MM_PER_INCH = 25.4

# Okabe & Ito の色覚バリアフリー8色（Color Universal Design / jfly）
# Nature Methods (Wong 2011) が紹介して以来、カテゴリ色の事実上の標準。
OKABE_ITO = {
    "black": "#000000",
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
}
# 白背景で線に使う順（黄色は視認性が低いので後ろ）
OKABE_ITO_CYCLE = [
    OKABE_ITO["blue"],
    OKABE_ITO["vermillion"],
    OKABE_ITO["bluish_green"],
    OKABE_ITO["orange"],
    OKABE_ITO["sky_blue"],
    OKABE_ITO["reddish_purple"],
    OKABE_ITO["black"],
    OKABE_ITO["yellow"],
]

# check_figure.py の JOURNALS と同期させること
WIDTHS_MM = {
    "nature": {"single": 89.0, "double": 183.0},
    "elsevier": {"single": 90.0, "onehalf": 140.0, "double": 190.0},
    "cellpress": {"single": 85.0, "onehalf": 114.0, "double": 174.0},
    "plos": {"single": 66.8, "double": 190.5},
    "springer": {"single": 84.0, "double": 174.0},
    # ACS: 単段 ≤240pt(84.7mm) / 2段 300-504pt(105.8-177.8mm)。2段には**下限**がある
    "acs": {"single": 84.7, "double": 177.8},
    # T&F は幅を公表していない。値は掲載論文の実測（Crit Rev Anal Chem 2026・US Letter 2段）。
    # ⚠️ 誌が変われば測り直す（FIGURE_GUIDE §7.8b）。T&F は図を段幅に**引き伸ばさない**ので、
    #    小さめは安全・大きめだけが縮小されて pt 規定を割る＝迷ったら小さく作る。
    "tandf": {"single": 89.6, "double": 182.8},
    "generic": {"single": 85.0, "double": 174.0},
}
BASE_FONT_PT = {
    "nature": 7.0,
    "elsevier": 8.0,
    "cellpress": 9.0,
    "plos": 9.0,
    "springer": 9.0,
    "acs": 7.0,      # 誌の下限は 4.5pt。読みやすさを優先して 7pt を基準にする
    "generic": 8.0,
}


def mm(value_mm: float) -> float:
    """mm → inch（matplotlib の figsize 用）"""
    return value_mm / MM_PER_INCH


def _pick_sans() -> list:
    """Arial → Helvetica → 環境にある sans の順で使う。無ければ DejaVu に落ちる。"""
    available = {f.name for f in matplotlib.font_manager.fontManager.ttflist}
    order = ["Arial", "Helvetica", "Helvetica Neue", "Liberation Sans", "Nimbus Sans"]
    picked = [name for name in order if name in available]
    return picked + ["DejaVu Sans"]


def use_journal_style(journal: str = "generic", base_pt: float | None = None) -> None:
    """rcParams を投稿仕様に固定する。図を作る前に一度だけ呼ぶ。"""
    pt = base_pt if base_pt else BASE_FONT_PT.get(journal, 8.0)
    # 🚨 mathtext を放置すると、本文が Arial でも `$x^2$` の部分だけ DejaVu で埋め込まれる
    #    （2026-08-04 実測：ACS 図の唯一の書体違反が mathtext 由来だった）。
    #    fontset='custom' にして本文と同じ sans を数式にも使う。
    sans = _pick_sans()
    head = sans[0]
    plt.rcParams.update(
        {
            "mathtext.fontset": "custom",
            "mathtext.rm": head,
            "mathtext.it": f"{head}:italic",
            "mathtext.bf": f"{head}:bold",
            "mathtext.sf": head,
            "mathtext.tt": head,
            "mathtext.cal": f"{head}:italic",
            # ---- 最重要：フォントを TrueType(42) で埋め込む（既定の Type3 を回避）
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            # ---- 書体
            "font.family": "sans-serif",
            "font.sans-serif": sans,
            # ---- 文字サイズ（掲載サイズでの実寸 pt）
            "font.size": pt,
            "axes.titlesize": pt,
            "axes.labelsize": pt,
            "xtick.labelsize": pt - 0.5,
            "ytick.labelsize": pt - 0.5,
            "legend.fontsize": pt - 0.5,
            "figure.titlesize": pt + 1,
            # ---- 線・枠（細すぎると印刷で消える：0.5pt を下限に）
            "axes.linewidth": 0.6,
            "grid.linewidth": 0.5,
            "lines.linewidth": 1.0,
            "lines.markersize": 3.0,
            "patch.linewidth": 0.6,
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.minor.width": 0.5,
            "ytick.minor.width": 0.5,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            # ---- 体裁
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.prop_cycle": plt.cycler(color=OKABE_ITO_CYCLE),
            "legend.frameon": False,
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "standard",  # tight は寸法が崩れるので使わない
            "savefig.transparent": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def new_figure(
    journal: str = "generic",
    column: str = "single",
    height_mm: float = 60.0,
    **kwargs,
):
    """掲載サイズちょうどの Figure を作る（縮小されない＝フォント規定を割らない）。"""
    width_mm = WIDTHS_MM.get(journal, WIDTHS_MM["generic"])[column]
    return plt.subplots(figsize=(mm(width_mm), mm(height_mm)), **kwargs)


# Editorial Manager（Elsevier / Cell Press）が表示項目の上端に重ね刷りするスタンプの
# 占有域。JMB-D-26-00962_R2.pdf を実測して確定（2026-08-03）。EM は図を拡大縮小せず、
# 組み上がりのページ寸法はアップロードしたファイルと同一なので、この座標は絶対値。
# 詳細と設計規約は FIGURE_GUIDE.md §6.4、機械判定は check_figure.py の「EMスタンプ帯」。
EM_STAMP_Y_MM = (1.49, 8.55)     # 横はページ幅の中央から右端まで
EM_LABEL_MM = (1.76, 1.49, 30.4, 5.02)   # 左上の item 名ラベル x0,y0,x1,y1


def warn_if_em_stamp_collides(pdf_path: str) -> bool:
    """書き出した PDF が EM のスタンプ域に文字・図形を置いていないかを実測する。

    作った瞬間に気づけるのが一番安い。ここでは例外を投げず警告に留める
    （投稿しない図・EM 以外へ出す図もあるため）。投稿前の判定は
    check_figure.py --journal elsevier の「EMスタンプ帯」（FAIL）で行う。
    """
    try:
        import fitz
    except ImportError:
        return False
    to_pt = lambda v: v / 25.4 * 72
    page = fitz.open(pdf_path)[0]
    w, h = page.rect.width, page.rect.height
    stamp = fitz.Rect(w / 2, to_pt(EM_STAMP_Y_MM[0]), w, to_pt(EM_STAMP_Y_MM[1]))
    label = fitz.Rect(*(to_pt(v) for v in EM_LABEL_MM))

    def touches(r):
        if r.height <= 0 or r.width <= 0:
            r = fitz.Rect(r.x0 - .05, r.y0 - .05, r.x1 + .05, r.y1 + .05)
        return r.intersects(stamp) or r.intersects(label)

    hit = []
    for b in page.get_text("dict")["blocks"]:
        for line in b.get("lines", []):
            for s in line.get("spans", []):
                if s.get("text", "").strip() and touches(fitz.Rect(s["bbox"])):
                    hit.append(s["text"].strip()[:24])
    for it in page.get_drawings():
        r = it["rect"]
        if r.width >= w * .99 and r.height >= h * .99:
            continue
        if touches(r):
            hit.append("(図形)")
    if hit:
        print(f"  ⚠️ {os.path.basename(pdf_path)}: EM のスタンプ域に {len(hit)} 件"
              f"（{'／'.join(repr(t) for t in hit[:3])}）。Elsevier/Cell Press へ出すなら、"
              f"上端 y {EM_STAMP_Y_MM[1]} mm・ページ幅の中央から右を空けること"
              "（FIGURE_GUIDE.md §6.4）")
    return bool(hit)


def warn_if_text_offpage(fig) -> list:
    """紙の外に落ちた文字を書き出し前に検出する。

    🚨 なぜ必要か（2026-08-04 実測）：軸ラベルが下マージンに収まらないと、matplotlib は
    **警告を出さずに紙の外へ描く**。PDF に痕跡が残らないので `check_figure.py`（文字サイズ・
    書体）も `check_figure_overlap.py`（重なり・はみ出し）も**そもそも見えない**＝両方緑になる。
    「消えた要素」は緑のチェックでは絶対に捕まらないので、作った瞬間にここで捕まえる。
    """
    fig.canvas.draw()
    W, H = fig.bbox.x1, fig.bbox.y1
    lost = []
    for txt in fig.findobj(match=lambda o: hasattr(o, "get_text")):
        s = (txt.get_text() or "").strip()
        if not s or not txt.get_visible():
            continue
        try:
            bb = txt.get_window_extent()
        except Exception:
            continue
        if bb.x1 < 0 or bb.y1 < 0 or bb.x0 > W or bb.y0 > H:
            lost.append(s[:28])
    if lost:
        print(f"  🚨 紙の外に落ちた文字 {len(lost)} 件（{'／'.join(repr(t) for t in lost[:4])}）"
              "＝出力ファイルには存在しない。余白か図の高さを増やすこと")
    return lost


def savefig(fig, stem: str, dpi: int = 600, formats=("pdf", "png")) -> list:
    """PDF（本投稿用ベクタ）と PNG（確認・プレビュー用）を同じ寸法で書き出す。"""
    os.makedirs(os.path.dirname(stem) or ".", exist_ok=True)
    warn_if_text_offpage(fig)
    written = []
    for ext in formats:
        path = f"{stem}.{ext}"
        fig.savefig(path, format=ext, dpi=dpi)
        written.append(path)
        if ext == "pdf":
            warn_if_em_stamp_collides(path)
    return written


def cvd_safe_pair(name_a: str = "blue", name_b: str = "vermillion") -> tuple:
    """2群比較の既定色。赤緑ではなく青×朱で、P/D/T いずれでも分離する。"""
    return OKABE_ITO[name_a], OKABE_ITO[name_b]


def _selftest(outdir: str = "/tmp/figstyle_selftest") -> str:
    """スタイルが実際に規定を満たす図を作れるかを検証する（check_figure.py で確認）。"""
    import numpy as np

    use_journal_style("elsevier")
    fig, ax = new_figure("elsevier", "single", height_mm=55)
    x = np.linspace(0, 10, 100)
    for i, (label, shift) in enumerate([("control", 0.0), ("treated", 0.6)]):
        ax.plot(x, np.sin(x) + shift, label=label, linestyle=["-", "--"][i])
    ax.set_xlabel("time (h)")
    ax.set_ylabel("relative level")
    ax.legend()
    fig.tight_layout(pad=0.3)
    stem = os.path.join(outdir, "selftest")
    paths = savefig(fig, stem)
    plt.close(fig)
    return paths[0]


if __name__ == "__main__":
    print("wrote:", _selftest())
