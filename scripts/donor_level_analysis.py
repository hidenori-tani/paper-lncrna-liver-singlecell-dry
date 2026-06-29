#!/usr/bin/env python3
"""Donor-level validation of the lncRNA detection-rate claims (JGG revision, Path A).

Addresses the Codex pseudoreplication critique WITHOUT re-running the clustering pipeline,
by using the ORIGINAL publication's own cell-type annotation (annot_humanAll.csv: `annot`,
`patient`=donor, `diet`=Lean/Obese). All statistics are computed at the DONOR level
(n=16 donors), not the cell level.

No scanpy / harmony / leiden needed — only gzip, numpy, pandas, scipy, statsmodels, matplotlib.

Outputs (submission/jgg_v1/donor_level/):
  per_donor_detection.csv        donor x cell_type x gene detection rate + cell counts
  donor_lsec_enrichment.csv      paired Wilcoxon: LSEC vs pooled-other detection, per gene
  lean_vs_obese_donor.csv        Mann-Whitney across donors (Lean vs Obese) per gene x cell_type
  crosslineage_fdr.csv           existing screen re-FDR'd across ALL lineage x gene tests
  fig_donor_detection.pdf/png    per-donor detection-rate strip plots (KCNQ1OT1, MEG3)
"""
from __future__ import annotations
import gzip, re, sys, zipfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import wilcoxon, mannwhitneyu
from statsmodels.stats.multitest import multipletests
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PROJ = Path(__file__).resolve().parent.parent  # repo root
ANCHOR = PROJ / "data" / "GSE192742"
ZIP = ANCHOR / "rawData_human.zip"
EXTRACT = ANCHOR / "rawData_human_extracted"
OUTDIR = PROJ / "submission" / "jgg_v1" / "donor_level"
OUTDIR.mkdir(parents=True, exist_ok=True)

TARGETS = ["KCNQ1OT1", "MEG3", "ZFAS1", "LINC00996", "NEAT1", "MALAT1",
           "IDI2-AS1", "MIAT", "HULC", "LINC00261"]


def log(*a): print(*a, flush=True)


def extract():
    marker = EXTRACT / ".extracted"
    if marker.exists():
        log(f"  already extracted at {EXTRACT}"); return
    EXTRACT.mkdir(parents=True, exist_ok=True)
    log(f"  extracting {ZIP} -> {EXTRACT} ...")
    with zipfile.ZipFile(ZIP) as zf:
        zf.extractall(EXTRACT)
    marker.touch()


def find_rna_dir() -> Path:
    for d in EXTRACT.rglob("*"):
        if not d.is_dir() or "ADT" in d.name:
            continue
        files = {f.name for f in d.iterdir() if f.is_file()}
        if ({"matrix.mtx.gz", "matrix.mtx"} & files) and \
           ({"barcodes.tsv.gz", "barcodes.tsv"} & files) and \
           ({"features.tsv.gz", "features.tsv"} & files):
            return d
    raise RuntimeError("RNA matrix dir not found")


def _open(p):
    return gzip.open(p, "rt") if str(p).endswith(".gz") else open(p)


def parse_sample_order(path) -> dict:
    samples = []
    with open(path) as fh:
        for line in fh:
            m = re.match(r"^\s*-\s+([A-Za-z0-9]+):", line)
            if m: samples.append(m.group(1))
    return {i + 1: s for i, s in enumerate(samples)}


def load_target_counts(rna_dir: Path):
    """Stream the genes x cells mtx, keep only TARGET gene rows -> cells x targets DataFrame."""
    feat = rna_dir / ("features.tsv.gz" if (rna_dir / "features.tsv.gz").exists() else "features.tsv")
    bc = rna_dir / ("barcodes.tsv.gz" if (rna_dir / "barcodes.tsv.gz").exists() else "barcodes.tsv")
    mtx = rna_dir / ("matrix.mtx.gz" if (rna_dir / "matrix.mtx.gz").exists() else "matrix.mtx")

    with _open(feat) as fh:
        genes = [ln.rstrip("\n").split("\t")[0] for ln in fh]
    with _open(bc) as fh:
        barcodes = [ln.strip() for ln in fh]
    log(f"  genes={len(genes)} barcodes={len(barcodes)}")

    # map TARGET gene symbol -> 1-based mtx row (handle var_names_make_unique duplicates: take first)
    gene_to_row = {}
    for i, g in enumerate(genes):
        if g in TARGETS and g not in gene_to_row:
            gene_to_row[g] = i + 1  # 1-based
    present = list(gene_to_row.keys())
    missing = [g for g in TARGETS if g not in gene_to_row]
    log(f"  targets present: {present}")
    if missing: log(f"  targets MISSING from matrix: {missing}")
    target_rows = set(gene_to_row.values())
    row_to_gene = {v: k for k, v in gene_to_row.items()}

    # header
    with _open(mtx) as fh:
        line = fh.readline()
        while line.startswith("%"): line = fh.readline()
        n_genes, n_cells, nnz = map(int, line.split())
        log(f"  mtx header genes={n_genes} cells={n_cells} nnz={nnz}")

    # stream the body in chunks, keep only target rows
    kept = []
    with _open(mtx) as fh:
        line = fh.readline()
        while line.startswith("%"): line = fh.readline()
        # next reads are body
        reader = pd.read_csv(fh, sep=" ", header=None, names=["row", "col", "val"],
                             dtype={"row": "int32", "col": "int32", "val": "int32"},
                             engine="c", chunksize=20_000_000)
        for k, chunk in enumerate(reader):
            sub = chunk[chunk["row"].isin(target_rows)]
            if len(sub): kept.append(sub)
            log(f"    chunk {k}: {len(chunk)} rows, kept {len(sub)} (cum {sum(len(x) for x in kept)})")
    body = pd.concat(kept, ignore_index=True)
    log(f"  total target entries: {len(body)}")

    # build cells x targets dense (targets are few)
    counts = pd.DataFrame(0, index=range(n_cells), columns=present, dtype="int32")
    for g in present:
        r = gene_to_row[g]
        e = body[body["row"] == r]
        counts.loc[e["col"].values - 1, g] = e["val"].values
    counts.index = barcodes
    return counts, barcodes


def main():
    log("[1] extract"); extract()
    rna_dir = find_rna_dir(); log(f"  RNA dir: {rna_dir}")
    # annot lives at ANCHOR/ (separate download); sampleComp is inside the extracted zip
    annot_path = ANCHOR / "annot_humanAll.csv"
    if not annot_path.exists():
        annot_path = next(EXTRACT.rglob("annot_humanAll.csv"))
    sc_path = next(EXTRACT.rglob("sampleComp_humanAll.txt"))
    log(f"  annot: {annot_path}\n  sampleComp: {sc_path}")

    log("[2] load target counts from mtx (streaming)")
    counts, barcodes = load_target_counts(rna_dir)

    log("[3] join publication annotation (cell_type=annot, donor=patient, condition=diet)")
    split = pd.Series(barcodes).str.rsplit("-", n=1, expand=True)
    bc_root = split[0].values
    aggr_idx = pd.to_numeric(split[1], errors="coerce")
    sample_order = parse_sample_order(sc_path)
    sample = pd.Series(aggr_idx).map(sample_order).values

    annot = pd.read_csv(annot_path)
    log(f"  annot columns: {list(annot.columns)}")
    annot["barcode_root"] = annot["cell"].str.rsplit("-", n=1).str[0]
    annot_idx = annot.set_index(["sample", "barcode_root"])
    keys = pd.MultiIndex.from_arrays([sample, bc_root], names=["sample", "barcode_root"])
    aligned = annot_idx.reindex(keys)

    df = counts.reset_index(drop=True)
    df["cell_type"] = aligned["annot"].values
    df["donor"] = aligned["patient"].values
    df["condition"] = aligned["diet"].values
    df["assay"] = aligned["typeSample"].values if "typeSample" in aligned.columns else "NA"
    df = df[df["cell_type"].notna() & df["donor"].notna()].copy()
    log(f"  assay (typeSample) distribution: {df['assay'].value_counts(dropna=False).to_dict()}")
    log(f"  annotated cells: {len(df)}")
    log(f"  cell types: {sorted(df['cell_type'].dropna().unique())}")
    log(f"  donors: {df['donor'].nunique()} | conditions: {df['condition'].value_counts(dropna=False).to_dict()}")

    present = [g for g in TARGETS if g in df.columns]

    # identify LSEC label
    cts = df["cell_type"].dropna().unique()
    lsec_label = next((c for c in cts if re.search(r"lsec|sinusoidal|endothel", str(c), re.I)), None)
    log(f"  LSEC label resolved to: {lsec_label!r}")

    # ---- per-donor detection rate per cell_type per gene ----
    log("[4] per-donor detection rates")
    rows = []
    for (donor, ct), g in df.groupby(["donor", "cell_type"]):
        cond = g["condition"].iloc[0]
        for gene in present:
            n = len(g); ndet = int((g[gene] > 0).sum())
            rows.append(dict(donor=donor, condition=cond, cell_type=ct, gene=gene,
                             n_cells=n, n_detected=ndet, detection_rate=ndet / n,
                             mean_count=float(g[gene].mean())))
    perdonor = pd.DataFrame(rows)
    perdonor.to_csv(OUTDIR / "per_donor_detection.csv", index=False)
    log(f"  wrote per_donor_detection.csv ({len(perdonor)} rows)")

    # per-donor LSEC cell counts (Codex #8)
    if lsec_label:
        lsec_n = df[df["cell_type"] == lsec_label].groupby("donor").size().sort_values(ascending=False)
        log("  per-donor LSEC cell counts:"); log(lsec_n.to_string())

    # ---- donor-level LSEC enrichment: paired Wilcoxon LSEC vs pooled-other ----
    log("[5] donor-level LSEC enrichment (paired Wilcoxon across donors)")
    enr = []
    if lsec_label:
        for gene in present:
            sub = perdonor[perdonor["gene"] == gene]
            lsec = sub[sub["cell_type"] == lsec_label].set_index("donor")["detection_rate"]
            # pooled-other detection rate per donor (cells outside LSEC)
            other = []
            for donor, g in df[df["cell_type"] != lsec_label].groupby("donor"):
                other.append((donor, (g[gene] > 0).mean()))
            other = pd.Series(dict(other))
            common = lsec.index.intersection(other.index)
            l = lsec.reindex(common).values; o = other.reindex(common).values
            # require LSEC present in donor (>=10 LSEC cells) for a fair paired test
            lsec_counts = df[df["cell_type"] == lsec_label].groupby("donor").size().reindex(common).fillna(0)
            keep = lsec_counts.values >= 10
            l2, o2 = l[keep], o[keep]
            if len(l2) >= 5:
                try:
                    stat, p = wilcoxon(l2, o2, alternative="greater")
                except ValueError:
                    stat, p = np.nan, np.nan
            else:
                stat, p = np.nan, np.nan
            enr.append(dict(gene=gene, n_donors=int(keep.sum()),
                            median_LSEC_detection=float(np.median(l2)) if len(l2) else np.nan,
                            median_other_detection=float(np.median(o2)) if len(o2) else np.nan,
                            wilcoxon_stat=stat, p_greater=p))
        enr = pd.DataFrame(enr)
        if len(enr) and enr["p_greater"].notna().any():
            enr["q_bh"] = np.nan
            m = enr["p_greater"].notna()
            enr.loc[m, "q_bh"] = multipletests(enr.loc[m, "p_greater"], method="fdr_bh")[1]
        enr.to_csv(OUTDIR / "donor_lsec_enrichment.csv", index=False)
        log(enr.to_string(index=False))

    # ---- Lean vs Obese donor-level test (Mann-Whitney across donors) ----
    log("[6] Lean vs Obese donor-level test (Mann-Whitney across donors)")
    lvo = []
    for gene in present:
        for ct in cts:
            sub = perdonor[(perdonor["gene"] == gene) & (perdonor["cell_type"] == ct)]
            # require >=20 cells of this type in the donor for a stable per-donor rate
            sub = sub[sub["n_cells"] >= 20]
            lean = sub[sub["condition"].str.contains("lean", case=False, na=False)]["detection_rate"]
            obese = sub[sub["condition"].str.contains("obes", case=False, na=False)]["detection_rate"]
            if len(lean) >= 3 and len(obese) >= 3:
                stat, p = mannwhitneyu(lean, obese, alternative="two-sided")
                lvo.append(dict(gene=gene, cell_type=ct, n_lean=len(lean), n_obese=len(obese),
                                median_lean=float(lean.median()), median_obese=float(obese.median()),
                                mwu_stat=float(stat), p=float(p)))
    lvo = pd.DataFrame(lvo)
    if len(lvo):
        lvo["q_bh"] = multipletests(lvo["p"], method="fdr_bh")[1]
        lvo = lvo.sort_values("p")
        lvo.to_csv(OUTDIR / "lean_vs_obese_donor.csv", index=False)
        log(lvo.head(20).to_string(index=False))
    else:
        log("  (no gene x cell_type with >=3 Lean and >=3 Obese donors at >=20 cells)")

    # ---- cross-lineage FDR re-adjustment of the existing screen ----
    log("[7] cross-lineage FDR re-adjustment of the original screen")
    screen_csv = PROJ / "submission" / "hepres_v1" / "supplementary" / "progression_lncrna_table.csv"
    if screen_csv.exists():
        s = pd.read_csv(screen_csv)
        s = s[s["wald_p"].notna()].copy()
        s["q_within_lineage"] = np.nan
        for lin, g in s.groupby("lineage"):
            s.loc[g.index, "q_within_lineage"] = multipletests(g["wald_p"], method="fdr_bh")[1]
        s["q_crosslineage"] = multipletests(s["wald_p"], method="fdr_bh")[1]
        hits = s[s["gene"].isin(["MEG3", "ZFAS1", "LINC00996"]) &
                 s["lineage"].isin(["LSEC", "Cholangiocyte"])]
        cols = ["gene", "lineage", "spearman_rho", "wald_p", "q_within_lineage", "q_crosslineage"]
        keep = hits[hits.apply(lambda r: (r["gene"] in ["MEG3","ZFAS1"] and r["lineage"]=="LSEC") or
                                          (r["gene"]=="LINC00996" and r["lineage"]=="Cholangiocyte"), axis=1)]
        keep[cols].to_csv(OUTDIR / "crosslineage_fdr.csv", index=False)
        log("  three flagged hits, within-lineage vs cross-lineage FDR:")
        log(keep[cols].to_string(index=False))
        log(f"  (total tests across all lineage x gene: {len(s)})")

    # ---- figure: per-donor detection rate strip plots for KCNQ1OT1 & MEG3 ----
    log("[8] figure")
    main_genes = [g for g in ["KCNQ1OT1", "MEG3"] if g in present]
    # order cell types by overall LSEC-first, then by mean detection
    fig, axes = plt.subplots(1, len(main_genes), figsize=(7.0 * len(main_genes), 5.0), squeeze=False)
    for ax, gene in zip(axes[0], main_genes):
        sub = perdonor[(perdonor["gene"] == gene) & (perdonor["n_cells"] >= 20)]
        order = (sub.groupby("cell_type")["detection_rate"].median().sort_values(ascending=False).index.tolist())
        for i, ct in enumerate(order):
            vals = sub[sub["cell_type"] == ct]["detection_rate"].values
            x = np.random.default_rng(0).normal(i, 0.06, size=len(vals))
            color = "#d62728" if (lsec_label and ct == lsec_label) else "#1f77b4"
            ax.scatter(x, vals, s=18, alpha=0.7, color=color, edgecolors="none")
            ax.hlines(np.median(vals), i - 0.25, i + 0.25, color="black", lw=1.5)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order, rotation=45, ha="right", fontsize=7)
        ax.set_ylabel("per-donor detection rate")
        ax.set_title(f"{gene}  (each point = 1 donor)", fontsize=10)
        ax.set_ylim(-0.02, 1.02)
    fig.tight_layout()
    fig.savefig(OUTDIR / "fig_donor_detection.pdf")
    fig.savefig(OUTDIR / "fig_donor_detection.png", dpi=200)
    log(f"  wrote fig_donor_detection.{{pdf,png}}")
    log("DONE.")


if __name__ == "__main__":
    sys.exit(main())
