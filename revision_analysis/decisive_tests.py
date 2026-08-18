"""R2-4の決定的検定 / R1-6 根の感度 / R3-3 擬時間とdietの関連 / R1-4 HVGベンチマーク."""
import sys, gzip, warnings, time, subprocess
warnings.filterwarnings("ignore")
sys.path.insert(0,"pipeline"); sys.path.insert(0,"revision_analysis")
import numpy as np, pandas as pd, scipy.sparse as sp, scanpy as sc, anndata as ad
from scipy.stats import rankdata, spearmanr, mannwhitneyu
from sklearn.metrics import adjusted_rand_score, adjusted_mutual_info_score
from pipeline.screening_core import screen_lineage
import harmonypy

SEED=42; np.random.seed(SEED); RA="revision_analysis"
CT="data/GSE192742/rawData_human_extracted/rawData_human/countTable_human"
T0=time.time()
def log(m): print(f"[{time.time()-T0:7.1f}s] {m}", flush=True)

with gzip.open(f"{CT}/features.tsv.gz","rt") as fh:
    all_genes=np.array([l.split("\t")[0].strip() for l in fh])
rows=np.loadtxt(f"{RA}/extract_rows.txt",dtype=np.int64)
hvg_f=set(np.loadtxt(f"{RA}/hvg_rows_forced.txt",dtype=np.int64).tolist())
hvg_s=set(np.loadtxt(f"{RA}/hvg_rows_standard.txt",dtype=np.int64).tolist())
gpos={r:i for i,r in enumerate(rows)}; gene_names=all_genes[rows-1]

lab=pd.read_csv(f"{RA}/global_labels.csv")           # 前回の leiden 結果を再利用
cells=pd.read_parquet(f"{RA}/cells_annotated.parquet").sort_values("col").reset_index(drop=True)
cells=cells.merge(lab[["col","leiden","lineage_leiden","lineage_pub"]],on="col",how="left")
cpos={c:i for i,c in enumerate(cells.col.values)}; NC=len(cells)

NNZ=int(subprocess.run(["wc","-l",f"{RA}/hvg_submatrix.tsv"],capture_output=True,text=True).stdout.split()[0])
ri=np.empty(NNZ,np.int32); ci=np.empty(NNZ,np.int32); dv=np.empty(NNZ,np.float32); o=0
for ch in pd.read_csv(f"{RA}/hvg_submatrix.tsv",sep="\t",header=None,names=["g","c","v"],
                      dtype={"g":np.int32,"c":np.int32,"v":np.float32},chunksize=20_000_000):
    n=len(ch); ri[o:o+n]=ch.c.map(cpos).values; ci[o:o+n]=ch.g.map(gpos).values; dv[o:o+n]=ch.v.values; o+=n
A=ad.AnnData(X=sp.csr_matrix((dv,(ri,ci)),shape=(NC,len(rows))),obs=cells,
             var=pd.DataFrame({"hvg_forced":[r in hvg_f for r in rows],
                               "hvg_standard":[r in hvg_s for r in rows]},index=gene_names))
del ri,ci,dv; A.var_names_make_unique()
A.X = sp.diags((1e4/A.obs.total_umi.values).astype(np.float32)) @ A.X
sc.pp.log1p(A); log(f"行列 {A.shape} 準備完了")

def embed(adata, hvg_col):
    adata.var["highly_variable"]=adata.var[hvg_col].values
    sc.tl.pca(adata,n_comps=50,use_highly_variable=True,zero_center=False,random_state=SEED)
    ho=harmonypy.run_harmony(adata.obsm["X_pca"],adata.obs,["patient"],max_iter_harmony=10)
    Z=np.asarray(ho.Z_corr); Z=Z.T if Z.shape[0]!=adata.n_obs else Z
    adata.obsm["X_pca_harmony"]=Z; return adata

A=embed(A,"hvg_forced"); log("PCA+Harmony（強制あり）完了")

GENES=["MEG3","ZFAS1","LINC00996"]
res_contam=[]; res_root=[]; res_diet=[]

for key,lin in [("lineage_leiden","LSEC"),("lineage_pub","LSEC"),
                ("lineage_leiden","Kupffer"),("lineage_pub","Kupffer")]:
    m=(A.obs[key]==lin).values
    if m.sum()<200: continue
    sub=A[m].copy()
    sc.pp.neighbors(sub,use_rep="X_pca_harmony",n_neighbors=15); sc.tl.diffmap(sub)
    lean=np.where(sub.obs.diet.astype(str).str.lower().str.contains("lean"))[0]
    sub.uns["iroot"]=int(np.random.RandomState(SEED).choice(lean)) if len(lean) else 0
    sc.tl.dpt(sub); pt=sub.obs.dpt_pseudotime.values.copy()
    log(f"{key}/{lin}: n={sub.n_obs:,} 擬時間 完了")

    # --- R2-4 決定的検定: 同じ擬時間のまま、公表アノテで内皮/該当細胞だけに絞る ---
    if lin=="LSEC":
        pure=(sub.obs.annot=="Endothelial cells").values
        for g in GENES:
            x=np.asarray(sub[:,g].X.todense()).ravel()
            r_all=spearmanr(x,pt)[0]
            r_pure=spearmanr(x[pure],pt[pure])[0]
            contam=~pure
            r_contam=spearmanr(x[contam],pt[contam])[0] if contam.sum()>10 else np.nan
            res_contam.append(dict(label_set=key,gene=g,n_all=len(x),n_pure=int(pure.sum()),
                rho_all=r_all,rho_pure_endothelial_only=r_pure,rho_contaminants_only=r_contam))
            log(f"  {g}: rho(全{len(x)})={r_all:+.3f} → rho(内皮{int(pure.sum())}のみ)={r_pure:+.3f}")

    # --- R1-6 根の感度（20通り）---
    scr=sub[:,GENES].copy()
    G=np.asarray(scr.X.todense()); Rg=np.apply_along_axis(rankdata,0,G)
    Rg=Rg-Rg.mean(0); Rn=np.linalg.norm(Rg,axis=0); Rn[Rn==0]=np.nan
    rr=[]
    for k in range(20):
        sub.uns["iroot"]=int(np.random.RandomState(1000+k).choice(lean)) if len(lean) else 0
        sc.tl.dpt(sub)
        rp=rankdata(sub.obs.dpt_pseudotime.values).astype(float); rp-=rp.mean()
        rr.append((Rg*rp[:,None]).sum(0)/(Rn*np.linalg.norm(rp)))
    R=pd.DataFrame(rr,columns=GENES)
    for g in GENES:
        res_root.append(dict(label_set=key,lineage=lin,gene=g,
            rho_main=spearmanr(np.asarray(sub[:,g].X.todense()).ravel(),pt)[0],
            rho_median=R[g].median(),rho_min=R[g].min(),rho_max=R[g].max(),
            sign_consistent=float((np.sign(R[g])==np.sign(R[g].median())).mean())))
    log(f"  根20通り 完了")

    # --- R3-3 ドナー水準: 擬時間の代表値 と diet の関連 ---
    dd=pd.DataFrame({"patient":sub.obs.patient.values,"diet":sub.obs.diet.values,"pt":pt})
    per=dd.groupby(["patient","diet"]).pt.agg(["mean","median","size"]).reset_index()
    lean_v=per[per.diet.str.lower()=="lean"]["mean"]; ob_v=per[per.diet.str.lower()=="obese"]["mean"]
    u=mannwhitneyu(ob_v,lean_v,alternative="two-sided") if len(ob_v)>=2 and len(lean_v)>=2 else None
    res_diet.append(dict(label_set=key,lineage=lin,n_donors=len(per),
        n_lean=len(lean_v),n_obese=len(ob_v),
        mean_pt_lean=lean_v.mean(),mean_pt_obese=ob_v.mean(),
        mannwhitney_p=(u.pvalue if u else np.nan)))
    per.to_csv(f"{RA}/donor_pseudotime_{key}_{lin}.csv",index=False)

pd.DataFrame(res_contam).to_csv(f"{RA}/R2-4_contamination_test.csv",index=False)
pd.DataFrame(res_root).to_csv(f"{RA}/R1-6_root_sensitivity.csv",index=False)
pd.DataFrame(res_diet).to_csv(f"{RA}/R3-3_donor_pseudotime_vs_diet.csv",index=False)
log("R2-4 / R1-6 / R3-3 書き出し完了")

# --- R1-4: 強制取り込みが系統分離に与えた影響（ARI/AMI で比較）---
log("R1-4 ベンチマーク開始（標準HVGで埋め込み→Leiden）")
B=A.copy(); B=embed(B,"hvg_standard")
sc.pp.neighbors(B,use_rep="X_pca_harmony",n_neighbors=15)
sc.tl.leiden(B,resolution=1.0,key_added="leiden_std",random_state=SEED)
log(f"標準HVG Leiden 完了 clusters={B.obs.leiden_std.nunique()}")
truth=A.obs.annot.values
bench=pd.DataFrame([
  dict(feature_selection="forced lncRNA included (n_HVG=4750)",
       n_clusters=A.obs.leiden.nunique(),
       ARI_vs_published=adjusted_rand_score(truth,A.obs.leiden.values),
       AMI_vs_published=adjusted_mutual_info_score(truth,A.obs.leiden.values)),
  dict(feature_selection="standard HVG only (n_HVG=4000)",
       n_clusters=B.obs.leiden_std.nunique(),
       ARI_vs_published=adjusted_rand_score(truth,B.obs.leiden_std.values),
       AMI_vs_published=adjusted_mutual_info_score(truth,B.obs.leiden_std.values)),
])
bench.to_csv(f"{RA}/R1-4_hvg_benchmark.csv",index=False)
print(bench.to_string(index=False))
B.obs[["col","leiden_std"]].to_csv(f"{RA}/leiden_standard_hvg.csv",index=False)
log("R1-4 書き出し完了")
