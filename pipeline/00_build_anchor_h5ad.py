"""Stage 0: Build the anchor AnnData from Liver Cell Atlas raw downloads.

Inputs (downloaded by data/GSE192742/download.sh):
- data/GSE192742/rawData_human.zip       # 10X filtered_feature_bc_matrix per donor
- data/GSE192742/annot_humanAll.csv      # QC-retained cell metadata

Output:
- data/GSE192742/human_allcells.h5ad     # concatenated, annotated, ready for QC
"""
from __future__ import annotations

import sys
import zipfile
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc

from pipeline import config


def extract_if_needed(zip_path: Path, target_dir: Path) -> Path:
    """Unzip rawData_human.zip if not already extracted. Returns extracted root dir."""
    marker = target_dir / ".extracted"
    if marker.exists():
        print(f"  already extracted at {target_dir}")
        return target_dir
    target_dir.mkdir(parents=True, exist_ok=True)
    print(f"  extracting {zip_path} → {target_dir}")
    with zipfile.ZipFile(zip_path) as zf:
        zf.extractall(target_dir)
    marker.touch()
    return target_dir


def find_rna_matrix_dir(root: Path) -> Path:
    """Identify the RNA 10X matrix directory.

    Guilliams 2022 zip layout (confirmed 2026-05-25):
      rawData_human/
        countTable_human/             ← RNA (this one)
          barcodes.tsv.gz, matrix.mtx.gz, features.tsv.gz
        countTableADT_human_*/        ← ADT (CITE-seq antibodies; skip)

    Returns the RNA matrix dir. Raises if not found.
    """
    for d in root.rglob("*"):
        if not d.is_dir():
            continue
        if "ADT" in d.name:  # skip antibody-capture matrices
            continue
        files = {f.name for f in d.iterdir() if f.is_file()}
        has_matrix = ("matrix.mtx.gz" in files) or ("matrix.mtx" in files)
        has_barcodes = ("barcodes.tsv.gz" in files) or ("barcodes.tsv" in files)
        has_features = ("features.tsv.gz" in files) or ("features.tsv" in files)
        if has_matrix and has_barcodes and has_features:
            return d
    raise RuntimeError(
        f"No RNA 10X matrix directory found under {root}. "
        "Expected: rawData_human/countTable_human/{barcodes,matrix,features}."
    )


def load_rna_matrix(dir_path: Path) -> ad.AnnData:
    """Load the single combined RNA 10X matrix.

    Guilliams 2022 uses 10X v2 format (features.tsv.gz has only gene_symbol column,
    no Ensembl ID). scanpy.read_10x_mtx assumes v3 and fails with KeyError: 1.
    scipy.io.mmread on 351M nonzero entries is very slow (~30+ min). We use
    pandas.read_csv (C-based) on the streamed gzipped Matrix Market body,
    which is ~10-20x faster.

    Donor information comes from annot_humanAll.csv (joined later).
    """
    import gzip
    import scipy.sparse as sp

    barcodes_path = dir_path / "barcodes.tsv.gz"
    features_path = dir_path / "features.tsv.gz"
    matrix_path = dir_path / "matrix.mtx.gz"

    print(f"  reading {features_path.name}...")
    with gzip.open(features_path, "rt") as fh:
        gene_symbols = [line.strip().split("\t")[0] for line in fh]
    print(f"  genes: {len(gene_symbols)}")

    print(f"  reading {barcodes_path.name}...")
    with gzip.open(barcodes_path, "rt") as fh:
        barcodes = [line.strip() for line in fh]
    print(f"  barcodes: {len(barcodes)}")

    print(f"  reading {matrix_path.name} (pandas C-based, ~2-5 min for 351M entries)...")
    # Parse the Matrix Market header (% comment lines + 1 dimension line)
    with gzip.open(matrix_path, "rt") as fh:
        line = fh.readline()
        while line.startswith("%"):
            line = fh.readline()
        n_rows, n_cols, n_nonzero = map(int, line.strip().split())
        print(f"  header: rows={n_rows}, cols={n_cols}, nonzero={n_nonzero}")

        # Read remaining lines with pandas (C parser, vectorized int parsing)
        df = pd.read_csv(
            fh, sep=" ", header=None,
            names=["row", "col", "val"],
            dtype={"row": "int32", "col": "int32", "val": "int32"},
            engine="c",
        )
    print(f"  parsed {len(df)} entries")

    # Construct sparse matrix (Matrix Market is 1-indexed → subtract 1)
    mat = sp.coo_matrix(
        (df["val"].values, (df["row"].values - 1, df["col"].values - 1)),
        shape=(n_rows, n_cols),
    ).tocsr()
    print(f"  matrix shape (genes x cells): {mat.shape}")

    # AnnData convention: rows = cells, cols = genes → transpose
    a = ad.AnnData(
        X=mat.T.tocsr(),
        obs=pd.DataFrame(index=pd.Index(barcodes, name="barcode")),
        var=pd.DataFrame(index=pd.Index(gene_symbols, name="gene_symbol")),
    )
    # De-duplicate gene_symbol names (scanpy convention)
    a.var_names_make_unique()
    a.obs["matrix_source"] = dir_path.name
    return a


def annotate_with_metadata(adata: ad.AnnData, annot_csv: Path) -> ad.AnnData:
    """Join annot_humanAll.csv to adata.obs via barcode-level merge.

    The annot file has one row per QC-retained cell with a `cell` column that
    matches `<donor>_<barcode>` or a similar composite. We try several keying
    strategies in order of likelihood.
    """
    df = pd.read_csv(annot_csv)
    print(f"  annot rows: {len(df)}, columns: {list(df.columns)[:15]}")

    # Build adata barcode keys to try matching:
    adata_idx = adata.obs.index.astype(str)
    donor = adata.obs["donor_id"].astype(str)

    # Strategy A: annot has a single key column matching adata's obs.index directly
    obs_index_set = set(adata_idx)
    for col in df.columns:
        if df[col].astype(str).isin(obs_index_set).mean() > 0.5:
            print(f"  matched annot column '{col}' → obs.index (direct)")
            df = df.set_index(col)
            new_obs = adata.obs.join(df, how="left")
            adata.obs = new_obs
            return adata

    # Strategy B: composite "<donor>_<barcode>"
    composite = donor.astype(str) + "_" + adata_idx
    composite_set = set(composite)
    for col in df.columns:
        if df[col].astype(str).isin(composite_set).mean() > 0.5:
            print(f"  matched annot column '{col}' → <donor>_<barcode>")
            df_idx = df.set_index(col)
            # Construct lookup
            adata.obs["_lookup"] = composite.values
            joined = adata.obs.merge(df_idx, left_on="_lookup", right_index=True, how="left")
            joined.index = adata.obs.index
            adata.obs = joined.drop(columns=["_lookup"])
            return adata

    print("  WARNING: could not auto-match annot CSV to obs barcodes. "
          "Attaching the raw CSV as adata.uns['annot_humanAll'] for manual reconciliation.")
    adata.uns["annot_humanAll"] = df
    return adata


def main() -> int:
    anchor_dir = config.ANCHOR_DIR
    zip_path = anchor_dir / "rawData_human.zip"
    annot_path = anchor_dir / "annot_humanAll.csv"

    if not zip_path.exists():
        raise FileNotFoundError(
            f"{zip_path} not found. Run data/GSE192742/download.sh first."
        )
    if not annot_path.exists():
        raise FileNotFoundError(
            f"{annot_path} not found. Run data/GSE192742/download.sh first."
        )

    extracted_root = anchor_dir / "rawData_human_extracted"
    extract_if_needed(zip_path, extracted_root)

    print("Identifying RNA 10X matrix directory...")
    rna_dir = find_rna_matrix_dir(extracted_root)
    print(f"  RNA matrix dir: {rna_dir}")

    print("Loading combined RNA matrix...")
    combined = load_rna_matrix(rna_dir)
    print(f"  shape: {combined.shape}")

    # Annotation is handled by pipeline/00b_annotate_anchor.py (Cell Ranger aggr suffix mapping)
    out = anchor_dir / "human_allcells.h5ad"
    combined.write_h5ad(out)
    print(f"Saved {out}")
    print(f"Final shape: {combined.shape}")
    print(f"obs columns: {list(combined.obs.columns)}")
    print("\nNext step: python -m pipeline.00b_annotate_anchor")
    return 0


if __name__ == "__main__":
    sys.exit(main())
