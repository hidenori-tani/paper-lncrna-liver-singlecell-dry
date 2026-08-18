"""R3-1 / R2-7(前半) / R2-1後半 / F4交絡 の実測.

入力: revision_analysis/target_gene_counts.tsv  (gene_row, cell_col, count)
      revision_analysis/cell_total_umi.tsv      (cell_col, total_umi)
出力: revision_analysis/ 配下の CSV と要約テキスト

検出率だけでなく「存在量」を出す (Reviewer 3-1):
  - normalized expression (CP10K) の平均
  - 検出細胞内での平均
  - ドナー疑似バルク (sum counts / sum total UMI)
すべて typeSample (scRnaSeq/citeSeq/nucSeq) で層別する (Reviewer 2-1, F4交絡).
"""
import gzip, re, sys
import numpy as np, pandas as pd
from scipy import stats

ROOT = "."
CT = f"{ROOT}/data/GSE192742/rawData_human_extracted/rawData_human/countTable_human"
SC = f"{ROOT}/data/GSE192742/rawData_human_extracted/rawData_human/sampleComp_humanAll.txt"
RA = f"{ROOT}/revision_analysis"

GENES = {586:"FCN3",2577:"PTPRC",7686:"ALB",8260:"LRAT",9793:"PDGFRB",10266:"HULC",
         13124:"LINC00996",16877:"IDI2-AS1",17478:"OIT3",18245:"KCNQ1OT1",19153:"NEAT1",
         19155:"MALAT1",21209:"DCN",21320:"STAB2",23355:"MEG3",27317:"COL1A1",
         28796:"LINC00261",29161:"ZFAS1",29708:"CLEC4G",31763:"MIAT"}

# --- 1. 細胞メタデータ（行列の列順） ---
with gzip.open(f"{CT}/barcodes.tsv.gz","rt") as fh:
    bc = [l.strip() for l in fh]
n_cells = len(bc)
print(f"barcodes: {n_cells:,}")

samples=[]
for line in open(SC):
    m = re.match(r"^\s*-\s+([A-Za-z0-9]+):", line)
    if m: samples.append(m.group(1))
sample_order = {i+1:s for i,s in enumerate(samples)}
print(f"samples in aggr order: {len(sample_order)}")

root, idx = zip(*(b.rsplit("-",1) for b in bc))
cells = pd.DataFrame({"col": np.arange(1, n_cells+1),
                      "barcode_root": root,
                      "aggr_idx": pd.to_numeric(idx, errors="coerce")})
cells["sample"] = cells.aggr_idx.map(sample_order)

annot = pd.read_csv(f"{ROOT}/data/GSE192742/annot_humanAll.csv")
annot["barcode_root"] = annot["cell"].str.rsplit("-", n=1).str[0]
cells = cells.merge(annot[["sample","barcode_root","annot","patient","typeSample","diet"]],
                    on=["sample","barcode_root"], how="left")
matched = cells.annot.notna().sum()
print(f"annot と一致した細胞: {matched:,} / {n_cells:,} ({100*matched/n_cells:.1f}%)")

# --- 2. 総UMI と 対象遺伝子カウント ---
tot = pd.read_csv(f"{RA}/cell_total_umi.tsv", sep="\t", header=None, names=["col","total_umi"])
cells = cells.merge(tot, on="col", how="left")
cells["total_umi"] = cells.total_umi.fillna(0)

tg = pd.read_csv(f"{RA}/target_gene_counts.tsv", sep="\t", header=None,
                 names=["gene_row","col","count"])
tg["gene"] = tg.gene_row.map(GENES)
wide = tg.pivot_table(index="col", columns="gene", values="count", aggfunc="sum").fillna(0)
cells = cells.merge(wide, on="col", how="left")
for g in GENES.values():
    if g not in cells.columns: cells[g] = 0.0
    cells[g] = cells[g].fillna(0.0)

d = cells[cells.annot.notna() & (cells.total_umi > 0)].copy()
print(f"解析対象（アノテ有・総UMI>0）: {len(d):,} cells")
d.to_parquet(f"{RA}/cells_annotated.parquet", index=False)

# --- 3. 細胞型 × typeSample の存在量指標 ---
rows=[]
for g in ["MEG3","KCNQ1OT1","ZFAS1","LINC00996","MALAT1","NEAT1"]:
    cp10k = d[g] / d.total_umi * 1e4
    for (ct, ts), sub in d.groupby(["annot","typeSample"], observed=True):
        c = sub[g].values; t = sub.total_umi.values
        det = c > 0
        rows.append(dict(gene=g, cell_type=ct, typeSample=ts, n_cells=len(sub),
            detection_rate=det.mean(),
            mean_cp10k=(c/t*1e4).mean(),
            mean_cp10k_in_detected=((c[det]/t[det]*1e4).mean() if det.any() else np.nan),
            pseudobulk_cp10k=(c.sum()/t.sum()*1e4) if t.sum()>0 else np.nan,
            n_donors=sub.patient.nunique()))
ab = pd.DataFrame(rows)
ab.to_csv(f"{RA}/abundance_by_celltype_and_assay.csv", index=False)
print(f"\n書き出し: abundance_by_celltype_and_assay.csv ({len(ab)} 行)")

# --- 4. MEG3: 内皮 vs 他 を測定法ごとに（ドナー水準・対応あり Wilcoxon） ---
print("\n=== MEG3 内皮 vs 非内皮（ドナー水準・測定法で層別）===")
out=[]
for ts, sub in d.groupby("typeSample", observed=True):
    pb = (sub.assign(is_endo=sub.annot.eq("Endothelial cells"))
             .groupby(["patient","is_endo"])
             .apply(lambda x: pd.Series({
                 "detection_rate": (x.MEG3>0).mean(),
                 "pseudobulk_cp10k": x.MEG3.sum()/x.total_umi.sum()*1e4,
                 "n_cells": len(x)}), include_groups=False)
             .reset_index())
    piv_d = pb.pivot(index="patient", columns="is_endo", values="detection_rate").dropna()
    piv_p = pb.pivot(index="patient", columns="is_endo", values="pseudobulk_cp10k").dropna()
    n = len(piv_d)
    line = {"typeSample": ts, "n_donors_paired": n}
    if n >= 3:
        for nm, piv in (("detection", piv_d), ("pseudobulk", piv_p)):
            try:
                w = stats.wilcoxon(piv[True], piv[False], alternative="greater")
                line[f"{nm}_median_endo"]=float(np.median(piv[True]))
                line[f"{nm}_median_other"]=float(np.median(piv[False]))
                line[f"{nm}_p_onesided"]=float(w.pvalue)
                line[f"{nm}_n_donors_endo_higher"]=int((piv[True]>piv[False]).sum())
            except Exception as e:
                line[f"{nm}_error"]=str(e)
    out.append(line)
res = pd.DataFrame(out)
res.to_csv(f"{RA}/meg3_endothelial_by_assay.csv", index=False)
print(res.to_string(index=False))
print("\n完了")
