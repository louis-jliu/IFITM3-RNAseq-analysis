"""03_extract_rosmap_microglia.py

Author : Jinrui (Louis) Liu, with Kushal Dey and Christina Leslie (Computational and Systems
         Biology Program, MSKCC). Contact: liuj11@mskcc.org

Extract microglia from the ROSMAP prefrontal cortex single-nucleus RNA-seq data
(Mathys et al., Nature 2019; AD Knowledge Portal, Synapse syn18485175) for the
IFITM3+ vs IFITM3- microglial program meta-analysis (Fig. 5d,e).

Input : data/rosmap/alzheimers.h5ad
        AnnData with 70,634 nuclei x 17,926 genes; raw counts in layers["counts"],
        log-normalised expression in X; obs columns include "broad.cell.type",
        "projid" and "DiseaseStatus" (no-, early-, late-pathology).
Output: output/rosmap_microglia_log_normalized.tsv   genes x cells
        output/rosmap_microglia_metadata.tsv         one row per cell

The ROSMAP data are distributed under the AD Knowledge Portal data use agreement.
Do not commit the input or the outputs of this script to a public repository.

Run from the repository root:  python scripts/03_extract_rosmap_microglia.py
"""

import argparse
import os

import anndata as ad
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--h5ad", default=os.path.join(
        os.environ.get("IFITM3_DATA_DIR", "data"), "rosmap", "alzheimers.h5ad"))
    parser.add_argument("--out", default=os.environ.get("IFITM3_OUT_DIR", "output"))
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    print(adata)

    mic = adata[adata.obs["broad.cell.type"] == "Mic", :].copy()
    print(f"microglia: {mic.n_obs} cells from {mic.obs['projid'].nunique()} individuals")
    print(mic.obs["DiseaseStatus"].value_counts())

    # log-normalised expression, genes x cells
    x = mic.X.toarray() if hasattr(mic.X, "toarray") else mic.X
    expr = pd.DataFrame(x, index=mic.obs_names, columns=mic.var_names).T
    expr.index.name = "gene"
    expr.to_csv(os.path.join(args.out, "rosmap_microglia_log_normalized.tsv"), sep="\t")

    meta = mic.obs.copy()
    meta.index.name = "TAG"
    meta.to_csv(os.path.join(args.out, "rosmap_microglia_metadata.tsv"), sep="\t")


if __name__ == "__main__":
    main()
