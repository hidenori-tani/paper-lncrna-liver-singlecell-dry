# lncRNA-aware single-cell re-analysis of the human Liver Cell Atlas

Analysis code and result tables for the manuscript:

> **MEG3 is a donor-reproducible endothelial lncRNA in the human liver single-cell atlas**
> Hidenori Tani — Department of Health Pharmacy, Yokohama University of Pharmacy
> ORCID: [0000-0001-6390-4136](https://orcid.org/0000-0001-6390-4136)

Single-author computational study. This repository contains the custom analysis
code, the intermediate result tables, and the donor-level validation analysis. It
does **not** redistribute the primary data, which are publicly available from GEO.

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

> Interpretation caveat: a 10X 3′ detection rate is not transcript localization or
> abundance. All cell-type signals reported here are candidate signals requiring
> orthogonal validation (e.g., RNAscope/FISH).

## Data availability

- **Primary data (not included here):** human Liver Cell Atlas, GEO **GSE192742**
  (Guilliams et al., 2022), available from <https://www.livercellatlas.org>.
- **Stability reference:** BRIC-seq half-life table (Tani et al., 2012, *Genome Res.*).
- **lncRNA annotation:** GENCODE v45 long non-coding RNA gene set.

`scripts/parallel_download.sh` fetches the raw atlas matrix into `data/`.

## Repository layout

```
pipeline/      analysis pipeline (Stages 0-5): matrix build, QC, atlas,
               pseudotime, screening, focused deep-dive, stability, validation
figures/       figure-generation scripts
scripts/       donor_level_analysis.py (donor-level validation),
               make_fig6.py, parallel_download.sh (data fetch)
results/       intermediate result tables (CSV); results/donor_level/ holds the
               per-donor validation tables
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

# 4) figures
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

## Citation

If you use this code, please cite the associated manuscript and this archived
release. Zenodo DOI: _to be inserted on first release_.

## License

Released under the MIT License — see [`LICENSE`](LICENSE).
