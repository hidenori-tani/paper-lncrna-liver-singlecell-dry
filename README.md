# lncRNA-aware single-cell re-analysis of the human Liver Cell Atlas

Analysis code and result tables for the manuscript:

> **MEG3 is a donor-reproducible endothelial lncRNA replicated across two human liver single-cell atlases**
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

Headline result: *MEG3* is a donor-reproducible endothelial (liver sinusoidal
endothelial cell, LSEC) lncRNA — per-donor median detection 78% in LSECs versus
≤ 2.9% in every other cell type; two-sided paired Wilcoxon q = 0.0098; higher in all
nine informative donors. *KCNQ1OT1* is depleted from hepatocytes but is not
LSEC-specific. A within-lineage pseudotime screen flagged exploratory candidates
(*MEG3*, *ZFAS1*, *LINC00996*) that did not validate as Lean-versus-Obese
differences at the donor level.

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
release. Zenodo DOI: [10.5281/zenodo.21025051](https://doi.org/10.5281/zenodo.21025051)
(concept DOI — resolves to the latest version).

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
