# lncRNA-aware single-cell re-analysis of the human Liver Cell Atlas

Analysis code and result tables for the manuscript:

> **An lncRNA-aware single-cell framework with donor-level validation identifies reproducible MEG3 enrichment in human liver sinusoidal endothelial cells**
> Hidenori Tani — Department of Health Pharmacy, Yokohama University of Pharmacy
> ORCID: [0000-0001-6390-4136](https://orcid.org/0000-0001-6390-4136)

Single-author computational study. This repository contains the custom analysis
code, the intermediate result tables, the donor-level validation analysis, and the
independent-cohort replication. It does **not** redistribute the primary data, which
are publicly available from GEO.

## Summary

GENCODE v45 long non-coding RNAs (lncRNAs) were mapped across seven hepatic cell
lineages of a published human Liver Cell Atlas (GEO accession **GSE192742**;
Guilliams et al., 2022) using an lncRNA-aware highly variable gene selection that
force-retains all detectable lncRNAs. The central cell-type claims were then tested
at the **donor level** (per-donor detection rates; paired Wilcoxon across donors)
using the original publication's own cell-type annotation, independently of our
clustering.

Headline result: *MEG3* is reproducibly enriched in liver sinusoidal endothelial
cells (LSECs) — per-donor median detection 78% in LSECs versus ≤ 2.9% in every other
cell type with at least three informative donors; two-sided paired Wilcoxon
q = 0.0098; higher in all nine informative donors. The enrichment also holds as
abundance rather than detection frequency (donor pseudobulk 6.27 versus 1.05 CP10K
for the next-ranked cell type).

> ### ⚠ What the 2026-08 revision withdraws
>
> Three claims made in the version first submitted are **withdrawn** after an audit
> of the lineage assignment against the source publication's own annotation:
>
> 1. **The three "progression-coupled" lncRNAs** (*MEG3*, *ZFAS1*, *LINC00996*).
>    Holding the pseudotime fixed and restricting to correctly annotated endothelial
>    cells, *MEG3* goes from ρ = −0.315 to +0.057 and *ZFAS1* from +0.324 to −0.091;
>    both reverse sign. The correlations came from contaminating cells.
> 2. **The *KCNQ1OT1* hepatocyte-depletion claim.** The lineage it was measured in
>    contains no annotated hepatocytes (95% of it is neutrophils); the matrix holds
>    only 24 annotated hepatocytes in total.
> 3. **The *LINC00996* cholangiocyte result.** The lineage it was measured in
>    contains no annotated cholangiocytes.
>
> The cause is a marker-score `idxmax` lineage assignment with no acceptance
> threshold. Purity against the published annotation is 77.3% for the LSEC lineage
> but **0%** for the lineages named "cholangiocyte" and "hepatocyte".
> `revision_analysis/` reproduces every one of these checks.

**Independent-cohort replication.** The *MEG3*-in-LSEC result reproduces in a fully
independent human liver atlas (MacParland et al. 2018; GEO **GSE115469**, 5 donors),
where per-donor LSEC detection exceeds pooled non-LSEC detection in all 5 donors
(median 48% vs 1.7%; one-sided paired Wilcoxon P = 0.031, the n = 5 minimum), with
central-venous LSECs the highest-detecting cell type. The honest controls also
reproduce (*KCNQ1OT1* not LSEC-specific; *ZFAS1* not LSEC-enriched). Script:
`scripts/replication_macparland.py`.

> Interpretation caveat: a 10X 3′ detection rate is not transcript localization or
> abundance. All cell-type signals reported here are candidate signals requiring
> orthogonal validation (e.g., RNAscope/FISH).

## Data availability

- **Primary data (not included here):** human Liver Cell Atlas, GEO **GSE192742**
  (Guilliams et al., 2022), available from <https://www.livercellatlas.org>.
- **Replication cohort (not included here):** MacParland et al. 2018 human liver
  atlas, GEO **GSE115469**.
- **Stability reference:** BRIC-seq half-life table (Tani et al., 2012, *Genome Res.*).
- **lncRNA annotation:** GENCODE v45 long non-coding RNA gene set.

`scripts/parallel_download.sh` fetches the raw atlas matrix into `data/`. The
replication cohort is fetched with:

```bash
mkdir -p data/GSE115469 && cd data/GSE115469
curl -sL -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/GSE115469_Data.csv.gz
curl -sL -O https://ftp.ncbi.nlm.nih.gov/geo/series/GSE115nnn/GSE115469/suppl/GSE115469_CellClusterType.txt.gz
```

## Repository layout

```
pipeline/      analysis pipeline (Stages 0-5): matrix build, QC, atlas,
               pseudotime, screening, focused deep-dive, stability, validation
figures/       figure-generation scripts
scripts/       donor_level_analysis.py (donor-level validation),
               replication_macparland.py (independent-cohort replication, GSE115469),
               make_fig6.py, parallel_download.sh (data fetch)
revision_analysis/  analyses added at the 2026-08 revision: the annotation audit,
               the re-screen under both label sets, the fixed-pseudotime
               contamination test, the root-sensitivity sweep, the abundance and
               dispersion measures, the feature-selection benchmark, and the
               scripts that draw every figure of the revised manuscript
lib/           figstyle.py (journal figure style used by the revision figures)
results/revision/   result tables produced by revision_analysis/
results/       intermediate result tables (CSV); results/donor_level/ holds the
               per-donor validation tables; results/replication_macparland/ holds
               the MacParland independent-cohort replication tables + figure
environment.yml, requirements.txt   reproducible environment
```

## Reproduction

```bash
conda env create -f environment.yml          # or: pip install -r requirements.txt
conda activate liver-singlecell-lncrna

# 1) fetch the public atlas matrix into data/
bash scripts/parallel_download.sh

# 2) run the pipeline (from the repository root)
python pipeline/00_build_anchor_h5ad.py
python pipeline/01_qc.py
python pipeline/02_atlas.py
python pipeline/03_pseudotime.py
python pipeline/04a_screening.py
python pipeline/04b_focused.py
python pipeline/04c_stability.py
python pipeline/05_validation.py

# 3) donor-level validation (clustering-independent)
python scripts/donor_level_analysis.py

# 4) independent-cohort replication (MacParland 2018, GSE115469; fetch first — see Data availability)
python scripts/replication_macparland.py

# 5) figures
python figures/build_all_figures.py
```

The random seed is fixed (42) throughout. Software: Python 3.13, scanpy 1.12,
anndata 0.12, harmonypy 2.0, statsmodels, scikit-learn, matplotlib, pandas.

## Result tables

| File | Contents |
|---|---|
| `results/progression_lncrna_table.csv` | per-lineage Spearman ρ, p/q, and log2 dynamic range for every lncRNA |
| `results/anchor_composition.csv` | anchor lncRNA cell-type detection fractions |
| `results/anchor_trajectory.csv` | anchor lncRNA per-bin pseudotime trajectories |
| `results/anchor_coexpression.csv` | top co-expressed genes per anchor lncRNA × lineage |
| `results/stability_x_pseudotime.csv` | flagged lncRNA peak bin × BRIC-seq stability class |
| `results/screening_sensitivity.csv` | screening threshold sensitivity sweep |
| `results/qc_report.csv` | cell/gene quality-control summary |
| `results/donor_level/per_donor_detection.csv` | per-donor × cell-type × lncRNA detection rate |
| `results/donor_level/donor_lsec_enrichment_twosided.csv` | two-sided paired Wilcoxon, LSEC vs pooled non-LSEC (q used in the manuscript) |
| `results/donor_level/lean_vs_obese_donor.csv` | donor-level Lean-vs-Obese Mann–Whitney U |
| `results/donor_level/crosslineage_fdr.csv` | cross-lineage FDR re-adjustment |
| `results/replication_macparland/macparland_lsec_enrichment.csv` | MacParland (GSE115469) per-gene LSEC vs non-LSEC, per-donor + one-sided paired Wilcoxon |
| `results/replication_macparland/macparland_per_donor_detection.csv` | MacParland per-donor × cell-type × lncRNA detection rate |
| `results/replication_macparland/macparland_celltype_ranking.csv` | MacParland per-gene cell-type ranking by median per-donor detection |

## Citation

If you use this code, please cite the associated manuscript and this archived
release. Zenodo concept DOI (resolves to the latest version):
[10.5281/zenodo.21025050](https://doi.org/10.5281/zenodo.21025050).

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
