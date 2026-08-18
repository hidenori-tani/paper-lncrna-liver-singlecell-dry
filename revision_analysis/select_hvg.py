"""gene_stats.tsv から HVG を選び、抽出する遺伝子行番号を書き出す。"""
import sys, gzip; sys.path.insert(0,"revision_analysis")
import numpy as np, pandas as pd
from hvg_from_stats import hvg_from_stats

N = 152559
CT = "data/GSE192742/rawData_human_extracted/rawData_human/countTable_human"
with gzip.open(f"{CT}/features.tsv.gz","rt") as fh:
    genes = np.array([l.split("\t")[0].strip() for l in fh])
NG = len(genes)

st = pd.read_csv("revision_analysis/gene_stats.tsv", sep="\t", header=None,
                 names=["g","cs","csq","s","sq","rs","rsq","n"])
cs = np.zeros(NG); csq = np.zeros(NG); rs = np.zeros(NG); rsq = np.zeros(NG); nd = np.zeros(NG)
idx = st.g.values - 1
cs[idx]=st.cs.values; csq[idx]=st.csq.values; rs[idx]=st.rs.values; rsq[idx]=st.rsq.values; nd[idx]=st.n.values

mean_cp = cs/N
var_cp  = (csq/N - mean_cp**2) * N/(N-1)

# 投稿版は QC で遺伝子を 32,738 → 24,189 に絞ってから HVG を選んでいる。
# scanpy の既定 filter_genes(min_cells=3) を再現し、同じ土俵で選ぶ。
expressed = nd >= 3
print(f"検出細胞>=3 の遺伝子: {expressed.sum():,}（投稿版の報告値 24,189 と比較）")
mean_f = np.where(expressed, mean_cp, 0.0)
var_f  = np.where(expressed, var_cp, 0.0)
hv, df = hvg_from_stats(mean_f, var_f, 4000)
hv &= expressed
print(f"HVG(top4000): {hv.sum()}")

tab = pd.read_csv("release_public/results/progression_lncrna_table.csv")
LNC = set(tab.gene.unique())
is_lnc = np.isin(genes, list(LNC))
print(f"lncRNA（投稿版と同一の863）が行列に存在: {is_lnc.sum()}")

forced = hv | is_lnc
print(f"HVG + 強制 lncRNA = {forced.sum()}（投稿版の報告値 4,751 と比較）")

# 系統割り当ての再現に要るマーカー遺伝子も抽出対象に含める（PCA には使わない）
LINEAGE_MARKERS = {
    "Hepatocyte": ["ALB","TF","APOB","CYP3A4","HNF4A"],
    "HSC": ["LRAT","ACTA2","COL1A1","COL3A1","PDGFRB"],
    "Kupffer": ["CD68","MARCO","VSIG4","CD163","TIMD4"],
    "LSEC": ["PECAM1","STAB2","CLEC4G","OIT3","FCN3"],
    "Cholangiocyte": ["KRT19","KRT7","EPCAM","SOX9","CFTR"],
    "Tcell": ["CD3D","CD3E","CD8A","CD4"],
    "Bcell": ["CD79A","MS4A1","IGHM"],
    "NK": ["KLRD1","NKG7","GNLY"],
}
allm = sorted({m for v in LINEAGE_MARKERS.values() for m in v})
is_marker = np.isin(genes, allm)
missing = [m for m in allm if m not in set(genes)]
print(f"マーカー {len(allm)} 件中、行列に存在 {is_marker.sum()} 件"
      + (f" / 不在: {missing}" if missing else ""))
print(f"  うち HVG に既に含まれる: {int((is_marker & forced).sum())} 件"
      f" / 追加抽出が要る: {int((is_marker & ~forced).sum())} 件")

extract = forced | is_marker
np.savetxt("revision_analysis/extract_rows.txt", np.where(extract)[0]+1, fmt="%d")
np.savetxt("revision_analysis/hvg_rows_forced.txt", np.where(forced)[0]+1, fmt="%d")
np.savetxt("revision_analysis/hvg_rows_standard.txt", np.where(hv)[0]+1, fmt="%d")
pd.DataFrame({"gene":genes,"is_marker":is_marker,"mean_cp10k":mean_cp,"var_cp10k":var_cp,
              "mean_raw":rs/N,"var_raw":(rsq/N-(rs/N)**2)*N/(N-1),
              "n_detected":nd,"is_lncRNA":is_lnc,
              "hvg_standard":hv,"hvg_forced":forced}).to_csv(
    "revision_analysis/gene_level_stats.csv", index=False)
print("書き出し: gene_level_stats.csv / hvg_rows_forced.txt / hvg_rows_standard.txt")
