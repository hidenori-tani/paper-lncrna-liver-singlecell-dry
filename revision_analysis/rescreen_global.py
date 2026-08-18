"""投稿版と同じ「アトラス全体の埋め込み」を再構築し、系統ラベルだけを差し替えて比較する。

比較の設計（交絡を切る）:
  埋め込み = アトラス全体の HVG→PCA(zero_center=False)→Harmony(donor_id)  ← 投稿版と同一
  ラベル   = (a) Leiden + マーカー idxmax（投稿版の再現＝陽性対照）
             (b) Guilliams 2022 の公表アノテ（R2-3 が求めるもの）
  これで「hit の消長」がラベルの差だけに帰属できる。
"""
import sys, gzip, warnings, time
warnings.filterwarnings("ignore")
sys.path.insert(0, "pipeline"); sys.path.insert(0, "revision_analysis")
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad
from pipeline.screening_core import screen_lineage

SEED = 42; np.random.seed(SEED)
RA = "revision_analysis"
CT = "data/GSE192742/rawData_human_extracted/rawData_human/countTable_human"
T0 = time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

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
PUB2LIN = {"Endothelial cells":"LSEC","Macrophages":"Kupffer",
           "Mono+mono derived cells":"Kupffer","Cholangiocytes":"Cholangiocyte",
           "Hepatocytes":"Hepatocyte","Fibroblasts":"HSC"}

# ---------- 1. 遺伝子・細胞 ----------
with gzip.open(f"{CT}/features.tsv.gz","rt") as fh:
    all_genes = np.array([l.split("\t")[0].strip() for l in fh])
rows = np.loadtxt(f"{RA}/extract_rows.txt", dtype=np.int64)      # 1-based
hvg_rows = set(np.loadtxt(f"{RA}/hvg_rows_forced.txt", dtype=np.int64).tolist())
gpos = {r:i for i,r in enumerate(rows)}
gene_names = all_genes[rows-1]
log(f"抽出遺伝子 {len(rows):,}（うち HVG/強制 {len(hvg_rows):,}）")

cells = pd.read_parquet(f"{RA}/cells_annotated.parquet")
cells = cells.sort_values("col").reset_index(drop=True)
cpos = {c:i for i,c in enumerate(cells.col.values)}
NC = len(cells); log(f"細胞 {NC:,}")

# ---------- 2. 疎行列（事前確保して1回で埋める） ----------
import subprocess
NNZ = int(subprocess.run(["wc","-l",f"{RA}/hvg_submatrix.tsv"],capture_output=True,text=True).stdout.split()[0])
log(f"非ゼロ {NNZ:,} を読み込み")
ri = np.empty(NNZ, np.int32); ci = np.empty(NNZ, np.int32); dv = np.empty(NNZ, np.float32)
o = 0
for ch in pd.read_csv(f"{RA}/hvg_submatrix.tsv", sep="\t", header=None,
                      names=["g","c","v"], dtype={"g":np.int32,"c":np.int32,"v":np.float32},
                      chunksize=20_000_000):
    n = len(ch)
    ri[o:o+n] = ch.c.map(cpos).values
    ci[o:o+n] = ch.g.map(gpos).values
    dv[o:o+n] = ch.v.values
    o += n; log(f"  {o:,}/{NNZ:,}")
X = sp.csr_matrix((dv,(ri,ci)), shape=(NC,len(rows)))
del ri, ci, dv
log(f"行列 {X.shape} 構築完了")

A = ad.AnnData(X=X, obs=cells.reset_index(drop=True),
               var=pd.DataFrame({"hvg": [r in hvg_rows for r in rows]}, index=gene_names))
A.var_names_make_unique()
del X

# ---------- 3. 正規化（全遺伝子の総UMIで割る＝投稿版と同じ）----------
inv = sp.diags((1e4/A.obs.total_umi.values).astype(np.float32))
A.X = inv @ A.X
sc.pp.log1p(A)
A.var["highly_variable"] = A.var.hvg.values
log("正規化・log1p 完了")

# ---------- 4. 全体の埋め込み ----------
sc.tl.pca(A, n_comps=50, use_highly_variable=True, zero_center=False, random_state=SEED)
log("PCA 完了")
import harmonypy
ho = harmonypy.run_harmony(A.obsm["X_pca"], A.obs, ["patient"], max_iter_harmony=10)
Z = np.asarray(ho.Z_corr); Z = Z.T if Z.shape[0]!=A.n_obs else Z
assert Z.shape==(A.n_obs,50), Z.shape
A.obsm["X_pca_harmony"] = Z; del ho, Z
log("Harmony 完了")
sc.pp.neighbors(A, use_rep="X_pca_harmony", n_neighbors=15)
log("neighbors 完了")
sc.tl.leiden(A, resolution=1.0, key_added="leiden", random_state=SEED)  # 投稿版と同一（flavor 既定=leidenalg）
log(f"Leiden 完了 clusters={A.obs.leiden.nunique()}")

# ---------- 5. ラベル2種 ----------
for name, mk in LINEAGE_MARKERS.items():
    valid = [m for m in mk if m in A.var_names]
    if valid: sc.tl.score_genes(A, valid, score_name=f"score_{name}")
sc_cols = [f"score_{n}" for n in LINEAGE_MARKERS if f"score_{n}" in A.obs.columns]
c2l = A.obs.groupby("leiden", observed=True)[sc_cols].mean().idxmax(axis=1).str.replace("score_","")
A.obs["lineage_leiden"] = A.obs.leiden.map(c2l).astype(str)
A.obs["lineage_pub"] = A.obs.annot.map(PUB2LIN).fillna("other").astype(str)
log("系統割り当て完了")
print("\n[Leiden+マーカー]\n", A.obs.lineage_leiden.value_counts().to_string())
print("\n[公表アノテ]\n", A.obs.lineage_pub.value_counts().to_string())

# ---------- 6. 系統ごとに擬時間＋スクリーニング ----------
tab = pd.read_csv("release_public/results/progression_lncrna_table.csv")
LNC = [g for g in sorted(tab.gene.unique()) if g in set(A.var_names)]
log(f"スクリーニング対象 lncRNA {len(LNC)}")

out = []
for key in ["lineage_leiden","lineage_pub"]:
    for lin in ["Hepatocyte","HSC","Kupffer","LSEC","Cholangiocyte"]:
        mask = (A.obs[key]==lin).values
        n = int(mask.sum())
        if n < 200:
            print(f"  {key}/{lin}: {n} cells → skip（投稿版と同じ 200 未満の除外）"); continue
        sub = A[mask].copy()
        sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_neighbors=15)
        sc.tl.diffmap(sub)
        lean = np.where(sub.obs.diet.astype(str).str.lower().str.contains("lean"))[0]
        sub.uns["iroot"] = int(np.random.RandomState(SEED).choice(lean)) if len(lean) else 0
        sc.tl.dpt(sub)
        scr = sub[:, LNC].copy()
        scr.obs["lineage"] = lin; scr.obs["pseudotime"] = sub.obs.dpt_pseudotime.values
        try:
            df = screen_lineage(scr, lin).reset_index()
        except ValueError as e:
            print(f"  {key}/{lin}: {e}"); continue
        df["label_set"] = key; df["lineage_name"] = lin; df["n_cells"] = n
        out.append(df)
        hits = df[df.is_progression_coupled]
        log(f"  {key}/{lin}: {n:,} cells → hit {len(hits)} 件 {list(hits.gene)}")

pd.concat(out).to_csv(f"{RA}/rescreen_global_table.csv", index=False)
A.obs[["col","patient","diet","annot","typeSample","leiden","lineage_leiden","lineage_pub"]].to_csv(
    f"{RA}/global_labels.csv", index=False)
log("書き出し完了: rescreen_global_table.csv / global_labels.csv")
