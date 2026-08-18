"""原著アノテ基準での再スクリーニング（R2-3）＋ 根の感度解析（R1-6）.

統計核は pipeline/pipeline/screening_core.py をそのまま使う（比較可能性のため）。
系統は Leiden 由来ではなく Guilliams 2022 の公表アノテを使う。
"""
import gzip, sys, warnings
import numpy as np, pandas as pd, scipy.sparse as sp
import scanpy as sc, anndata as ad
warnings.filterwarnings("ignore")
sys.path.insert(0, "pipeline")
from pipeline.screening_core import screen_lineage

SEED = 42
np.random.seed(SEED)
RA = "revision_analysis"
CT = "data/GSE192742/rawData_human_extracted/rawData_human/countTable_human"

# --- 遺伝子名 ---
with gzip.open(f"{CT}/features.tsv.gz","rt") as fh:
    genes = np.array([l.split("\t")[0].strip() for l in fh])
print(f"genes: {len(genes):,}")

# --- スクリーニング対象 lncRNA（投稿版と同一の 863 遺伝子）---
tab = pd.read_csv("release_public/results/progression_lncrna_table.csv")
LNC = sorted(tab.gene.unique())
print(f"lncRNA (投稿版と同一): {len(LNC)}")

# --- 部分行列 ---
cells = pd.read_csv(f"{RA}/screen_cells_pubannot.csv")
cols = np.sort(cells.col.values)
col_pos = {c:i for i,c in enumerate(cols)}
print(f"cells: {len(cols):,}  読み込み中...")
sm = pd.read_csv(f"{RA}/submatrix_pubannot.tsv", sep="\t", header=None,
                 names=["g","c","v"], dtype={"g":np.int32,"c":np.int32,"v":np.float32})
print(f"  非ゼロ: {len(sm):,}")
X = sp.csr_matrix((sm.v.values,
                   (sm.c.map(col_pos).values.astype(np.int32), sm.g.values-1)),
                  shape=(len(cols), len(genes)))
del sm
obs = cells.set_index("col").loc[cols].reset_index()
obs.index = obs.col.astype(str)
A = ad.AnnData(X=X, obs=obs, var=pd.DataFrame(index=genes))
A.var_names_make_unique()
print(f"AnnData: {A.shape}")

sc.pp.normalize_total(A, target_sum=1e4); sc.pp.log1p(A)

lnc_present = [g for g in LNC if g in set(A.var_names)]
print(f"  行列に存在する lncRNA: {len(lnc_present)}/{len(LNC)}")

results, sens = [], []
for lin in ["LSEC","Kupffer"]:
    sub = A[A.obs.lineage_pub == lin].copy()
    print(f"\n=== {lin}: {sub.n_obs:,} cells / {sub.obs.patient.nunique()} donors ===")
    sc.pp.highly_variable_genes(sub, n_top_genes=4000)
    forced = sub.var_names.isin(lnc_present)
    sub.var["highly_variable"] = sub.var.highly_variable.values | forced
    print(f"  HVG {int(sub.var.highly_variable.sum()):,}（うち強制 lncRNA {int(forced.sum())}）")
    emb = sub[:, sub.var.highly_variable].copy()
    sc.tl.pca(emb, n_comps=50, zero_center=False, random_state=SEED)
    import harmonypy
    ho = harmonypy.run_harmony(emb.obsm["X_pca"], emb.obs, ["patient"], max_iter_harmony=10)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != emb.n_obs:   # 版によって (PC, cell) / (cell, PC) が入れ替わる
        Z = Z.T
    assert Z.shape == (emb.n_obs, 50), Z.shape
    sub.obsm["X_pca_harmony"] = Z
    sc.pp.neighbors(sub, use_rep="X_pca_harmony", n_neighbors=15)
    sc.tl.diffmap(sub)

    lean_idx = np.where(sub.obs.diet.astype(str).str.lower().str.contains("lean"))[0]
    print(f"  Lean 細胞 {len(lean_idx):,} → 根の候補")

    # --- 主解析（投稿版と同じ：Lean からランダム1点・seed 42）---
    rng = np.random.RandomState(SEED)
    sub.uns["iroot"] = int(rng.choice(lean_idx))
    sc.tl.dpt(sub)
    pt_main = sub.obs["dpt_pseudotime"].values.copy()

    scr = sub[:, lnc_present].copy()
    scr.obs["lineage"] = lin; scr.obs["pseudotime"] = pt_main
    df = screen_lineage(scr, lin).reset_index(); df["lineage"] = lin
    results.append(df)
    hit = df[df.is_progression_coupled]
    print(f"  → progression-coupled: {len(hit)} 件 {list(hit.gene)}")

    # --- R1-6: 根の感度解析（Lean から 20 通り・行列演算）---
    from scipy.stats import rankdata
    G = scr.X.toarray() if hasattr(scr.X, "toarray") else np.asarray(scr.X)
    Rg = np.apply_along_axis(rankdata, 0, G)                  # 遺伝子の順位は1回だけ
    Rg = Rg - Rg.mean(0); Rg_n = np.linalg.norm(Rg, axis=0); Rg_n[Rg_n == 0] = np.nan
    rho_mat = []
    for k in range(20):
        sub.uns["iroot"] = int(np.random.RandomState(1000+k).choice(lean_idx))
        sc.tl.dpt(sub)
        rp = rankdata(sub.obs["dpt_pseudotime"].values).astype(float)
        rp = rp - rp.mean()
        rho_mat.append(pd.Series((Rg * rp[:, None]).sum(0) / (Rg_n * np.linalg.norm(rp)),
                                 index=lnc_present))
    R = pd.DataFrame(rho_mat)
    s = pd.DataFrame({"gene":lnc_present, "lineage":lin,
                      "rho_main":df.set_index("gene").loc[lnc_present,"spearman_rho"].values,
                      "rho_roots_median":R.median().values,
                      "rho_roots_min":R.min().values, "rho_roots_max":R.max().values,
                      "sign_consistent_frac":(np.sign(R)==np.sign(R.median())).mean().values})
    sens.append(s)
    print(f"  根20通りの感度解析 完了")

pd.concat(results).to_csv(f"{RA}/rescreen_pubannot_table.csv", index=False)
pd.concat(sens).to_csv(f"{RA}/root_sensitivity.csv", index=False)
print("\n書き出し: rescreen_pubannot_table.csv / root_sensitivity.csv")
