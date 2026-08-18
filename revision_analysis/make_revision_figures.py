"""改訂で新設・差し替えとなる図を作る（Springer 寸法・figstyle 経由）。

Fig 2 差し替え: 系統割り当ての監査（R2-3）
Fig 3 差し替え: スクリーニングが偽陽性だった機序（R2-4・R1-6・R3-3）
Fig 5 差し替え: 存在量指標（R3-1）
"""
import pathlib, sys, warnings; warnings.filterwarnings("ignore")
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent / "lib"))
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from figstyle import use_journal_style, new_figure, savefig, OKABE_ITO as C

use_journal_style("springer")
RA="revision_analysis"; OUT="submission/funig_r1/figures"
import os; os.makedirs(OUT, exist_ok=True)

lab=pd.read_csv(f"{RA}/global_labels.csv")
cells=pd.read_parquet(f"{RA}/cells_annotated.parquet")
annot=pd.read_csv("data/GSE192742/annot_humanAll.csv")
annot["barcode_root"]=annot["cell"].str.rsplit("-",n=1).str[0]
u=cells[["col","sample","barcode_root"]].merge(
    annot[["sample","barcode_root","UMAP_1","UMAP_2"]],on=["sample","barcode_root"],how="left")
d=lab.merge(u[["col","UMAP_1","UMAP_2"]],on="col",how="left").dropna(subset=["UMAP_1"])
print(f"UMAP 座標つき細胞: {len(d):,}")

# ───────── Fig 2: 系統割り当ての監査 ─────────
fig,axes=new_figure("springer","double",height_mm=150,nrows=2,ncols=2)
KEY={"Endothelial cells":C["blue"],"Cholangiocytes":C["vermillion"],
     "Fibroblasts":C["bluish_green"],"Hepatocytes":"#000000"}
ax=axes[0,0]
ax.scatter(d.UMAP_1,d.UMAP_2,s=0.05,c="#D9D9D9",rasterized=True,linewidths=0)
for lab_,col in KEY.items():
    m=(d.annot==lab_).values
    if m.sum(): ax.scatter(d.UMAP_1[m],d.UMAP_2[m],s=1.6,c=col,rasterized=True,linewidths=0,label=f"{lab_} ({m.sum():,})")
ax.set_title("A  Published annotation",loc="left"); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("UMAP 1"); ax.set_ylabel("UMAP 2"); ax.legend(markerscale=3,loc="lower left")

ax=axes[0,1]
KEY2={"LSEC":C["blue"],"Cholangiocyte":C["vermillion"],"Hepatocyte":"#000000"}
ax.scatter(d.UMAP_1,d.UMAP_2,s=0.05,c="#D9D9D9",rasterized=True,linewidths=0)
for lab_,col in KEY2.items():
    m=(d.lineage_leiden==lab_).values
    if m.sum(): ax.scatter(d.UMAP_1[m],d.UMAP_2[m],s=1.6,c=col,rasterized=True,linewidths=0,label=f'"{lab_}" ({m.sum():,})')
ax.set_title("B  Marker-score lineage (this study)",loc="left"); ax.set_xticks([]); ax.set_yticks([])
ax.set_xlabel("UMAP 1"); ax.legend(markerscale=3,loc="lower left")

s2=pd.read_csv("submission/funig_r1/supplementary/TableS2b_lineage_purity_recall.csv")
s2=s2[s2.purity_pct.notna()].sort_values("n_cells",ascending=False)
ax=axes[1,0]
x=np.arange(len(s2)); w=0.36
ax.bar(x-w/2,s2.purity_pct,w,color=C["blue"],label="Purity")
ax.bar(x+w/2,s2.recall_pct,w,color=C["vermillion"],hatch="///",edgecolor="white",linewidth=0,label="Recall")
ax.set_xticks(x); ax.set_xticklabels(s2.assigned_lineage,rotation=20,ha="right")
ax.set_ylabel("Agreement with published\nannotation (%)"); ax.set_ylim(0,108)
ax.legend(loc="upper right"); ax.set_title("C  Purity and recall",loc="left")
for xi,(pu,re_) in enumerate(zip(s2.purity_pct,s2.recall_pct)):
    if pu==0: ax.text(xi,4,"0 / 0",ha="center",va="bottom",fontsize=8.5)

ax=axes[1,1]
cm=pd.read_csv("submission/funig_r1/supplementary/TableS2_lineage_contingency.csv",index_col=0)
SHADES=["#0072B2","#7FB3D5","#BFD7EA","#EEEEEE"]
for lin,y in [("Cholangiocyte",1),("Hepatocyte",0)]:
    row=cm.loc[lin].sort_values(ascending=False)[:4]; left=0
    for i,(k,v) in enumerate(row.items()):
        ax.barh(y,v,left=left,color=SHADES[i],edgecolor="white",linewidth=0.5)
        if v>700: ax.text(left+v/2,y,f"{k}\n{v:,}",ha="center",va="center",fontsize=8,
                          color=("white" if i==0 else "black"))
        left+=v
ax.set_yticks([0,1]); ax.set_yticklabels(['Assigned\n"Hepatocyte"','Assigned\n"Cholangiocyte"'])
ax.set_xlabel("Cells (by published label)")
ax.set_title("D  What the two failed lineages contain",loc="left")
fig.subplots_adjust(left=0.115,right=0.97,top=0.955,bottom=0.10,wspace=0.30,hspace=0.42)
savefig(fig,f"{OUT}/Figure_2_lineage_audit"); plt.close(fig)
print("Fig 2 完了")

# ───────── Fig 3: 偽陽性の機序 ─────────
ct=pd.read_csv(f"{RA}/R2-4_contamination_test.csv"); ct=ct[ct.label_set=="lineage_leiden"]
rs=pd.read_csv(f"{RA}/R1-6_root_sensitivity.csv")
fig,axes=new_figure("springer","double",height_mm=78,ncols=3)
ax=axes[0]
genes=["MEG3","ZFAS1"]; x=np.arange(len(genes)); w=0.26
for i,(col,lbl,c) in enumerate([("rho_all","All 1,190 cells",C["vermillion"]),
                                ("rho_pure_endothelial_only","920 endothelial only",C["blue"]),
                                ("rho_contaminants_only","270 contaminants only",C["bluish_green"])]):
    v=[ct[ct.gene==g][col].iloc[0] for g in genes]
    ax.bar(x+(i-1)*w,v,w,color=c,label=lbl)
ax.axhline(0,color="black",lw=0.6); ax.axhline(0.3,ls=":",lw=0.6,color="#666666")
ax.axhline(-0.3,ls=":",lw=0.6,color="#666666")
ax.set_xticks(x); ax.set_xticklabels(genes,style="italic")
ax.set_ylabel("Spearman ρ vs pseudotime"); ax.legend(loc="lower left")
ax.set_title("A  Same pseudotime, purer cells",loc="left")

ax=axes[1]
for i,g in enumerate(genes):
    r=rs[(rs.label_set=="lineage_leiden")&(rs.lineage=="LSEC")&(rs.gene==g)].iloc[0]
    ax.plot([i,i],[r.rho_min,r.rho_max],color=C["blue"],lw=1.6,solid_capstyle="round")
    ax.plot(i,r.rho_median,"o",color=C["blue"],ms=4)
    ax.plot(i,r.rho_main,"D",color=C["vermillion"],ms=4)
    ax.text(i+0.13,r.rho_max,f"{100*r.sign_consistent:.0f}%\nsame sign",fontsize=8,va="top")
ax.axhline(0,color="black",lw=0.6)
ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes,style="italic"); ax.set_xlim(-0.5,1.7)
ax.set_ylabel("ρ across 20 alternative roots")
ax.set_title("B  Root sensitivity",loc="left")

ax=axes[2]
per=pd.read_csv(f"{RA}/donor_pseudotime_lineage_pub_Kupffer.csv")
for i,(dt,c) in enumerate([("Lean",C["blue"]),("Obese",C["vermillion"])]):
    v=per[per.diet.str.lower()==dt.lower()]["mean"]
    ax.scatter(np.random.RandomState(i).normal(i,0.06,len(v)),v,s=14,color=c,zorder=3)
    ax.plot([i-0.2,i+0.2],[v.mean()]*2,color="black",lw=1.0)
ax.set_xticks([0,1]); ax.set_xticklabels(["Lean\n(n = 12)","Obese\n(n = 3)"])
ax.set_ylabel("Donor mean pseudotime"); ax.set_title("C  No donor-level association",loc="left")
ax.text(0.5,0.95,"P = 0.73",transform=ax.transAxes,ha="center",va="top",fontsize=8)
fig.subplots_adjust(left=0.095,right=0.975,top=0.92,bottom=0.16,wspace=0.46)
savefig(fig,f"{OUT}/Figure_3_why_the_screen_failed"); plt.close(fig)
print("Fig 3 完了")

# ───────── Fig 5: 存在量 ─────────
ab=pd.read_csv(f"{RA}/abundance_by_celltype_and_assay.csv")
m=ab[(ab.gene=="MEG3")&(ab.typeSample=="scRnaSeq")&(ab.n_cells>=100)].sort_values("pseudobulk_cp10k",ascending=False)
fig,axes=new_figure("springer","double",height_mm=88,ncols=2)
ax=axes[0]
y=np.arange(len(m))[::-1]
ax.barh(y,m.pseudobulk_cp10k,color=[C["blue"] if c=="Endothelial cells" else "#BBBBBB" for c in m.cell_type])
ax.set_yticks(y); ax.set_yticklabels(m.cell_type)
ax.set_xscale("log"); ax.set_xlim(1e-3,40); ax.set_xticks([1e-3,1e-2,1e-1,1,10]); ax.set_xticklabels(['0.001','0.01','0.1','1','10']); ax.minorticks_off(); ax.set_xlabel("Donor pseudobulk expression (CP10K)")
ax.set_title("A  MEG3 abundance by cell type",loc="left")
for yi,v in zip(y,m.pseudobulk_cp10k):
    ax.text(v*1.25,yi,f"{v:.3g}",va="center",fontsize=8)

ax=axes[1]
mg=pd.read_csv(f"{RA}/meg3_endothelial_by_assay.csv")
mg=mg[mg.n_donors_paired>=3]
xs=np.arange(len(mg))
for i,(_,r) in enumerate(mg.iterrows()):
    ax.plot([i-0.16,i+0.16],[r.pseudobulk_median_other,r.pseudobulk_median_endo],
            "-o",color=C["blue"],ms=4,lw=1.0)
    ax.text(i,r.pseudobulk_median_endo*1.5,f"{int(r.pseudobulk_n_donors_endo_higher)}/{int(r.n_donors_paired)}\nP = {r.pseudobulk_p_onesided:.3g}",
            ha="center",fontsize=8)
ax.set_xticks(xs); ax.set_xticklabels([f"{t}\n(n = {int(n)} donors)" for t,n in zip(mg.typeSample,mg.n_donors_paired)])
ax.set_yscale("log"); ax.set_ylim(1e-3,40); ax.set_yticks([1e-3,1e-2,1e-1,1,10]); ax.set_yticklabels(['0.001','0.01','0.1','1','10']); ax.minorticks_off()
ax.set_ylabel("Donor pseudobulk MEG3 (CP10K)")
ax.set_xlim(-0.5,len(mg)-0.5)
ax.set_title("B  Endothelial vs non-endothelial",loc="left")
ax.text(0.02,0.5,"non-endothelial",transform=ax.transAxes,fontsize=8,color="#666666")
fig.subplots_adjust(left=0.27,right=0.955,top=0.93,bottom=0.14,wspace=0.60)
savefig(fig,f"{OUT}/Figure_5_abundance"); plt.close(fig)
print("Fig 5 完了")
