"""Figure 4 の B・C パネル用に、アンカー lncRNA の擬時間軌跡を「再実行」で取り直す。

なぜ必要か
----------
v1 の `anchor_trajectory.csv` は投稿版の実行から出ている。改訂本文は
「All numbers in this revision come from the re-execution」と述べているので、
図だけ投稿版の数値を出すわけにいかない。

やること＝`rescreen_global.py` と同じ手順を、同じ種で再現する。ただし
Leiden（12 分）は `global_labels.csv` に保存済みのラベルを読み込んで省略する
（クラスタ番号も系統も保存されているので、埋め込みさえ作り直せば同一になる）。
擬時間の根の選び方・10 分位のビン分けも `screening_core._bin_means` と同一。
"""
import sys, gzip, warnings, time; warnings.filterwarnings("ignore")
sys.path.insert(0, ".")
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad
import subprocess

SEED = 42
RA = "revision_analysis"
CT = "data/GSE192742/rawData_human_extracted/rawData_human/countTable_human"
GENES = ["MEG3", "NEAT1", "MALAT1", "KCNQ1OT1"]
LINEAGES = ["Hepatocyte", "Kupffer", "LSEC", "Cholangiocyte"]
T0 = time.time()


def log(m):
    print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)


with gzip.open(f"{CT}/features.tsv.gz", "rt") as fh:
    all_genes = np.array([l.split("\t")[0].strip() for l in fh])
rows = np.loadtxt(f"{RA}/extract_rows.txt", dtype=np.int64)
hvg_rows = set(np.loadtxt(f"{RA}/hvg_rows_forced.txt", dtype=np.int64).tolist())
gpos = {r: i for i, r in enumerate(rows)}
gene_names = all_genes[rows - 1]

cells = pd.read_parquet(f"{RA}/cells_annotated.parquet").sort_values("col").reset_index(drop=True)
cpos = {c: i for i, c in enumerate(cells.col.values)}
NC = len(cells)
log(f"細胞 {NC:,} / 抽出遺伝子 {len(rows):,}")

NNZ = int(subprocess.run(["wc", "-l", f"{RA}/hvg_submatrix.tsv"],
                         capture_output=True, text=True).stdout.split()[0])
ri = np.empty(NNZ, np.int32); ci = np.empty(NNZ, np.int32); dv = np.empty(NNZ, np.float32)
o = 0
for ch in pd.read_csv(f"{RA}/hvg_submatrix.tsv", sep="\t", header=None,
                      names=["g", "c", "v"], dtype={"g": np.int32, "c": np.int32, "v": np.float32},
                      chunksize=20_000_000):
    n = len(ch)
    ri[o:o+n] = ch.c.map(cpos).values
    ci[o:o+n] = ch.g.map(gpos).values
    dv[o:o+n] = ch.v.values
    o += n
X = sp.csr_matrix((dv, (ri, ci)), shape=(NC, len(rows)))
del ri, ci, dv
A = ad.AnnData(X=X, obs=cells.reset_index(drop=True),
               var=pd.DataFrame({"hvg": [r in hvg_rows for r in rows]}, index=gene_names))
A.var_names_make_unique()
del X
log("行列構築完了")

inv = sp.diags((1e4 / A.obs.total_umi.values).astype(np.float32))
A.X = inv @ A.X
sc.pp.log1p(A)
A.var["highly_variable"] = A.var.hvg.values
sc.tl.pca(A, n_comps=50, use_highly_variable=True, zero_center=False, random_state=SEED)
log("PCA 完了")
import harmonypy
ho = harmonypy.run_harmony(A.obsm["X_pca"], A.obs, ["patient"], max_iter_harmony=10)
Z = np.asarray(ho.Z_corr); Z = Z.T if Z.shape[0] != A.n_obs else Z
assert Z.shape == (A.n_obs, 50), Z.shape
A.obsm["X_pca_harmony"] = Z; del ho, Z
sc.pp.neighbors(A, use_rep="X_pca_harmony", n_neighbors=15)
log("Harmony + neighbors 完了")

# Leiden は再計算せず、保存済みの再実行ラベルを使う（同じ埋め込み・同じ種＝同一）
lab = pd.read_csv(f"{RA}/global_labels.csv").set_index("col")
assert lab.index.is_unique and len(lab) == NC, (len(lab), NC)
A.obs["lineage_leiden"] = A.obs.col.map(lab.lineage_leiden).astype(str).values
got = A.obs.lineage_leiden.value_counts().to_dict()
log(f"系統ラベル読み込み: {got}")

out = []
for lin in LINEAGES:
    mask = (A.obs.lineage_leiden == lin).values
    n = int(mask.sum())
    sub = A[mask].copy()
    sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.diffmap(sub)
    lean = np.where(sub.obs.diet.astype(str).str.lower().str.contains("lean"))[0]
    sub.uns["iroot"] = int(np.random.RandomState(SEED).choice(lean)) if len(lean) else 0
    sc.tl.dpt(sub)
    pt = sub.obs.dpt_pseudotime.values
    valid = ~np.isnan(pt)
    bins = pd.qcut(pt[valid], 10, labels=False, duplicates="drop")
    g = [x for x in GENES if x in sub.var_names]
    M = sub[valid, g].X
    M = M.toarray() if hasattr(M, "toarray") else np.asarray(M)
    df = pd.DataFrame(M, columns=g)
    df["_bin"] = bins
    bm = df.groupby("_bin").mean()
    bc = df.groupby("_bin").size()
    for gene in g:
        for b in bm.index:
            out.append({"gene": gene, "lineage": lin, "bin": int(b),
                        "mean_expr": float(bm.loc[b, gene]), "n_cells_bin": int(bc[b]),
                        "n_cells": n})
    log(f"  {lin}: {n:,} cells / 有効 {valid.sum():,} / ビン {len(bm)}")

pd.DataFrame(out).to_csv(f"{RA}/anchor_trajectory_reexec.csv", index=False)
log("書き出し: anchor_trajectory_reexec.csv")
