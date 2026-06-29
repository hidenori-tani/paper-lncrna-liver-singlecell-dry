"""Stage 4a: Unbiased lncRNA progression-coupled screening per lineage."""
from __future__ import annotations

import sys

import pandas as pd
import scanpy as sc

from pipeline import config, lncrna_utils, screening_core


def main():
    a = sc.read_h5ad(config.RESULTS_DIR / "liver_anchor_pseudotime.h5ad")
    a_lnc = lncrna_utils.filter_anndata_to_lncrna(a)
    print(f"lncRNA matrix shape: {a_lnc.shape}")

    all_results = []
    for lin in config.PSEUDOTIME["lineages"]:
        print(f"\nScreening {lin}...")
        try:
            df = screening_core.screen_lineage(
                a_lnc,
                lineage=lin,
                n_bins=config.SCREENING["n_bins"],
                spearman_abs_min=config.SCREENING["spearman_abs_min"],
                wald_q_max=config.SCREENING["wald_q_max"],
                log2fc_abs_min=config.SCREENING["log2fc_abs_min"],
                min_cells_per_bin=config.SCREENING["min_cells_per_bin"],
            )
        except ValueError as e:
            print(f"  skip: {e}")
            continue
        df["lineage"] = lin
        n_coupled = int(df["is_progression_coupled"].sum())
        print(f"  progression-coupled lncRNA: {n_coupled}")
        all_results.append(df.reset_index())

    out = pd.concat(all_results, ignore_index=True)
    out.to_csv(config.RESULTS_DIR / "progression_lncrna_table.csv", index=False)
    print(f"\nTotal progression-coupled lncRNA (across all lineages, deduplicated): "
          f"{out[out['is_progression_coupled']]['gene'].nunique()}")

    rows = []
    for rho_min in config.SCREENING_SENSITIVITY["spearman_abs_min"]:
        for lfc_min in config.SCREENING_SENSITIVITY["log2fc_abs_min"]:
            for lin in config.PSEUDOTIME["lineages"]:
                try:
                    d = screening_core.screen_lineage(
                        a_lnc, lineage=lin,
                        n_bins=config.SCREENING["n_bins"],
                        spearman_abs_min=rho_min,
                        wald_q_max=config.SCREENING["wald_q_max"],
                        log2fc_abs_min=lfc_min,
                        min_cells_per_bin=config.SCREENING["min_cells_per_bin"],
                    )
                    rows.append({
                        "lineage": lin, "rho_min": rho_min, "lfc_min": lfc_min,
                        "n_coupled": int(d["is_progression_coupled"].sum()),
                    })
                except ValueError:
                    pass
    pd.DataFrame(rows).to_csv(config.RESULTS_DIR / "screening_sensitivity.csv", index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
