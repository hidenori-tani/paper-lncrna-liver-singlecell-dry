"""GENCODE v45 lncRNA helpers.

Loads the GENCODE v45 long_noncoding_RNAs GTF, extracts gene IDs/names/biotypes.
Caches parsed result as `data/gencode/gencode_v45_lncrna_index.csv`.
"""
import gzip
import re
from functools import lru_cache
from pathlib import Path

import pandas as pd

from pipeline import config

LNCRNA_BIOTYPES = {
    "lncRNA",
    "lincRNA",
    "antisense_RNA",
    "antisense",
    "sense_intronic",
    "sense_overlapping",
    "processed_transcript",
    "macro_lncRNA",
    "bidirectional_promoter_lncRNA",
}

GTF_FILE = config.GENCODE_DIR / "gencode.v45.long_noncoding_RNAs.gtf.gz"
INDEX_CACHE = config.GENCODE_DIR / "gencode_v45_lncrna_index.csv"


def _parse_gtf_attribute(field: str, key: str) -> str | None:
    m = re.search(rf'{key} "([^"]+)"', field)
    return m.group(1) if m else None


def build_index() -> pd.DataFrame:
    """Parse GTF -> DataFrame with columns gene_id, gene_name, gene_type, chrom, start, end."""
    if not GTF_FILE.exists():
        raise FileNotFoundError(
            f"GENCODE v45 lncRNA GTF not found at {GTF_FILE}. "
            "Run the curl download from Task 0.5 Step 1."
        )

    rows = []
    with gzip.open(GTF_FILE, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = parts[8]
            gid = _parse_gtf_attribute(attrs, "gene_id")
            gname = _parse_gtf_attribute(attrs, "gene_name")
            gtype = (
                _parse_gtf_attribute(attrs, "gene_type")
                or _parse_gtf_attribute(attrs, "gene_biotype")
            )
            if gid is None or gtype not in LNCRNA_BIOTYPES:
                continue
            # Strip Ensembl version suffix
            gid_clean = gid.split(".")[0]
            rows.append(
                {
                    "gene_id": gid_clean,
                    "gene_name": gname,
                    "gene_type": gtype,
                    "chrom": parts[0],
                    "start": int(parts[3]),
                    "end": int(parts[4]),
                }
            )
    df = pd.DataFrame(rows)
    df.to_csv(INDEX_CACHE, index=False)
    return df


@lru_cache(maxsize=1)
def load_index() -> pd.DataFrame:
    if INDEX_CACHE.exists():
        return pd.read_csv(INDEX_CACHE)
    return build_index()


def get_lncrna_ids() -> list[str]:
    return load_index()["gene_id"].tolist()


def get_lncrna_name_mapping() -> dict[str, str]:
    """gene_id -> gene_name."""
    df = load_index()
    return dict(zip(df["gene_id"], df["gene_name"]))


def get_lncrna_name_to_id() -> dict[str, str]:
    """gene_name -> gene_id (first occurrence; gene_name is not always unique)."""
    df = load_index().dropna(subset=["gene_name"])
    return df.drop_duplicates("gene_name").set_index("gene_name")["gene_id"].to_dict()


def filter_anndata_to_lncrna(adata):
    """Return view of AnnData restricted to GENCODE v45 lncRNA biotype.

    AnnData var_names must be gene_ids (ENSG without version) or gene_names.
    Tries both lookup directions.
    """
    ids = set(get_lncrna_ids())
    names = set(get_lncrna_name_mapping().values())
    if any(v.startswith("ENSG") for v in adata.var_names[:5]):
        mask = adata.var_names.str.replace(r"\.\d+$", "", regex=True).isin(ids)
    else:
        mask = adata.var_names.isin(names)
    return adata[:, mask].copy()
