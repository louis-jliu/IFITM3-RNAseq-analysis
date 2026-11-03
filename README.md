# IFITM3-RNAseq-analysis

Analysis code and data accompanying the manuscript

> **The innate immunity protein IFITM3 modulates microglial states in Alzheimer's disease**
>
> Jain T.\*, Frost G.R.\*, Hur J.-Y., Chang J.C., Wong E., Liu C., Nikain C., Kehs Z., Cross A.C.,
> Li T., Jonas L.A., Liu L., Niu X., Schlachetzki J., Gorodetski L., Lim W.L., Yu G., Wang M.,
> Zhou X., Wang W., Li F., Sanford J., Khobrekar N.V., Evangelisti A., Sole Domenech S.,
> Miloushev V.Z., Funk C., Mazutis L., Lekaye C., Zou Y., Glass C.K., Keshari K.R., Cruchaga C.,
> Zhang B., Dey K., Leslie C., Studer L., Li Y.-M.
>
> *Nature Aging* (accepted, 2026). DOI: _to be added on publication_
>
> \* equal contribution. Corresponding authors: Yue-Ming Li (liy2@mskcc.org), Georgia R. Frost.

This repository contains the RNA-seq analysis code for the paper:

| Analysis | Authors | Panels | Code |
|---|---|---|---|
| Mouse bulk RNA-seq: Ifitm3-KO x 5xFAD interaction model and gene set enrichment | L. Liu, K. Dey, C. Leslie (MSKCC) | Fig. 2h-k; Extended Data Fig. 2g | `scripts/01`-`02`, `05` |
| Human microglia (ROSMAP snRNA-seq): IFITM3+ vs IFITM3- gene program meta-analysis | L. Liu, K. Dey, C. Leslie | Fig. 5e | `scripts/00`, `03`-`05` |
| Mouse snRNA-seq: clustering, annotation, MAST differential expression, CellChat | T. Jain | Fig. 3a-d; Fig. 4a-i; Extended Data Figs. 3-4 | `scripts/snRNA/` |

Every script names its author in its header. Questions about an analysis can go to its
author, to the corresponding author Yue-Ming Li (liy2@mskcc.org), or as a GitHub issue.

Repository maintainer: Jinrui (Louis) Liu (liuj11@mskcc.org),
ORCID [0000-0003-0778-4701](https://orcid.org/0000-0003-0778-4701)

## Data

- Mouse bulk RNA-seq (raw reads, raw counts): GEO [GSE353156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353156)
- Mouse snRNA-seq (raw reads, count matrices, cell annotations): GEO [GSE353157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353157)
- Human snRNA-seq: ROSMAP, AD Knowledge Portal (see `data/README.md`)

## Contents

```
data/       inputs (bulk counts, sample sheet, MSigDB gene sets, curated GO terms, gene programs)
results/    the exact result tables behind the manuscript figures
scripts/    bulk RNA-seq and human microglia scripts (run in order)
scripts/snRNA/  mouse snRNA-seq scripts (T. Jain)
```

| Script | What it does | Output |
|---|---|---|
| `scripts/00_build_gene_program_table.py` | Documents how `data/gene_programs/microglia_programs_Gan2024.tsv` was built from Gan et al. 2024 Table S1 (already done) | gene program table |
| `scripts/01_prepare_bulk_counts.R` | Gene x sample count matrix for the 15 mouse PFC samples | `output/bulk_*_N15.rds` |
| `scripts/02_bulk_deseq2_fgsea.R` | DESeq2 with an Ifitm3-KO x 5xFAD interaction term, and fgsea (MSigDB WikiPathways, GO) for every contrast | `output/DEres_*.csv`, `output/fgsea_*` |
| `scripts/03_extract_rosmap_microglia.py` | Extracts microglia from the ROSMAP snRNA-seq data | `output/rosmap_microglia_*.tsv` |
| `scripts/04_microglia_program_meta_analysis.R` | IFITM3+ vs IFITM3- microglia: per-gene effect sizes and random-effects meta-analysis per program and pathology stage | `output/microglia_program_meta_analysis_results.txt` |
| `scripts/05_make_figures.R` | Draws Fig. 2h-k, Ext. Data Fig. 2g and Fig. 5e | `output/figures/*.pdf` |

## Running

Run every script from the repository root.

Redraw the published panels from the exact result tables in `results/`:

```bash
Rscript scripts/05_make_figures.R
```

Re-run the full analysis, then draw the figures from the new results:

```bash
Rscript scripts/01_prepare_bulk_counts.R
Rscript scripts/02_bulk_deseq2_fgsea.R
python  scripts/03_extract_rosmap_microglia.py        # needs data/rosmap/alzheimers.h5ad (see data/README.md)
Rscript scripts/04_microglia_program_meta_analysis.R
IFITM3_FIG_INPUT=output Rscript scripts/05_make_figures.R
```

Set `IFITM3_DATA_DIR` or `IFITM3_OUT_DIR` to use other input or output folders.

### Mouse snRNA-seq (`scripts/snRNA/`, T. Jain)

These scripts were converted from T. Jain's notebooks (Jupyter and R Markdown) with only file
paths made relative and plots saved to files. They start from the processed AnnData object
(not included because of its size; see `data/README.md`).

| Script | What it does | Panels |
|---|---|---|
| `scripts/snRNA/01_snRNA_clustering_annotation_DE.py` | PhenoGraph clustering, cell type annotation, marker / IFITM3 / IL33 / IL34 plots, MAST DE (5xFAD;Ifitm3-/- vs 5xFAD) in microglia and oligodendrocytes | Fig. 3a, 3b, 3d, 4b-f, 4h, 4i; Ext. Data Figs. 3-4 |
| `scripts/snRNA/02_snRNA_volcano_plots.py` | Volcano plots of the MAST results | Fig. 3c, 4g |
| `scripts/snRNA/03_snRNA_cellchat.R` | CellChat "Secreted Signaling" analysis per genotype | Fig. 4a |

`01` uses T. Jain's MAST wrapper `run_MAST.r` and `02` her plotting helper `workshop_extra.py`;
these helpers are not part of this repository.

## Software

- R (>= 4.3) with DESeq2, fgsea, dplyr, ggplot2, ggrepel and rmeta
  (`environment/install_R_packages.R`). Scripts 02 and 04 save `sessionInfo()` to `output/`.
- Python (>= 3.9) with anndata and pandas (`environment/requirements.txt`).
- snRNA-seq: scanpy, phenograph, rpy2/anndata2ri with R MAST, scran and scuttle; CellChat (v1) and Seurat in R.

## Methods summary

**Bulk RNA-seq (Fig. 2h-k, Ext. Data Fig. 2g)**
- Mouse prefrontal cortex, 15 samples: WT n = 4, Ifitm3 KO n = 4, 5xFAD n = 4, 5xFAD;Ifitm3 KO n = 3.
- Filtering: mitochondrial and RIKEN (`Rik`) genes are removed, along with genes with <= 10
  total counts.
- Model: DESeq2 with `~ strain + condition + strain:condition` (reference: WT, Ifitm3 WT).
  Contrasts are listed in `02_bulk_deseq2_fgsea.R`.
- Differential expression: BH-adjusted p < 0.05 (with |log2FC| > log2(1.5) for Up/Down labels).
- GSEA: genes are ranked by log2 fold change and run through fgsea (minSize 15, maxSize 500)
  with MSigDB mouse v2024.1 WikiPathways and GO. Pathways with padj < 0.05 are reported.
- Ext. Data Fig. 2g shows a curated subset of synaptic and neurotransmitter GO terms
  (`data/EDFig2g_curated_GO_terms.txt`).
- In Fig. 2h-j, |log2FC| is capped at 5 for display and the Pearson r is computed on the capped values.

**Human microglia (Fig. 5e)**
- Data: ROSMAP snRNA-seq (Mathys et al., 2019), 1,920 microglia from 48 individuals.
- Microglia with log-normalised IFITM3 > 0 are IFITM3+.
- For each gene, the difference in mean expression (IFITM3+ minus IFITM3-, among cells that
  express the gene) and its SE are computed per pathology stage.
- Genes are combined per program with `rmeta::meta.summaries(method = "random")`.
- BH correction is applied across programs within each stage.
- tau < 0 means lower program expression in IFITM3+ microglia.

## Reproducibility notes

- fgsea estimates p values by sampling. `02_bulk_deseq2_fgsea.R` fixes the seed, but the
  submitted results were produced without one. A re-run therefore gives NES/p values within
  sampling error of `results/`; the figures are drawn from `results/` by default.
- ggrepel label positions are random; `05_make_figures.R` fixes the seed.

## Citation

If you use this code or data, please cite the manuscript above and the archived release of this
repository (Zenodo: https://doi.org/10.5281/zenodo.23196949, which always resolves to the latest version). Citation metadata is in `CITATION.cff`.

## License

Code: MIT (see `LICENSE`). Data files keep their original licences (see `data/README.md`).
