# Input data

| File | Description | Source / licence |
|---|---|---|
| `bulk/countData.csv` | Gene x sample count matrix, mouse prefrontal cortex bulk RNA-seq, 15 samples (TruSeq stranded mRNA, HiSeq 4000 PE50, MSK Integrated Genomics Operation). Column 1 is `GeneSymbol`; sample columns are named like `s_G_1__WT`. Raw reads and raw counts: GEO [GSE353156](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353156). | This study |
| `bulk/sample_metadata.csv` | Sample sheet: id, sex, strain (WT / 5xFAD), condition (noKO / KO), genotype | This study |
| `gene_sets/m2.cp.wikipathways.v2024.1.Mm.symbols.gmt` | MSigDB mouse WikiPathways collection, v2024.1.Mm | MSigDB (Broad Institute / UC San Diego), CC BY 4.0. https://www.gsea-msigdb.org |
| `gene_sets/m5.go.v2024.1.Mm.symbols.gmt` | MSigDB mouse Gene Ontology collection, v2024.1.Mm | as above |
| `EDFig2g_curated_GO_terms.txt` | GO terms shown in Extended Data Fig. 2g | This study |
| `gene_programs/microglia_programs_Gan2024.tsv` | Microglial gene programs (ATM, IRM, TGFbeta-dependent, DAM, LRRC33-dependent) taken from Table S1 of Gan et al., *Cancer Cell* 42, 1693-1712 (2024), doi:10.1016/j.ccell.2024.08.015, which compiles the original publications. Built with `scripts/00_build_gene_program_table.py`. | Gan et al. 2024 |

## Not included: ROSMAP single-nucleus RNA-seq

`scripts/03_extract_rosmap_microglia.py` expects `data/rosmap/alzheimers.h5ad`: the ROSMAP
prefrontal cortex snRNA-seq data from Mathys et al., *Nature* 570, 332-337 (2019) (70,634
nuclei, 48 individuals). Raw counts are in `layers["counts"]`, log-normalised values in `X`, and
cell annotations (`broad.cell.type`, `projid`, `DiseaseStatus`) in `obs`.

These data are distributed through the AD Knowledge Portal (Synapse syn18485175) under a data
use agreement and cannot be redistributed here. The summary statistics used for Fig. 5e are in
`results/microglia_program_meta_analysis_results.txt`.

## Not included: processed mouse snRNA-seq object

`scripts/snRNA/01_snRNA_clustering_annotation_DE.py` expects `data/snRNA/adata_snRNA_PFC.h5ad`,
the processed AnnData used for Fig. 3 and 4 (24,025 nuclei, 5xFAD and 5xFAD;Ifitm3-/-, 15 months).
It is not included because of its size (7.4 GB) and is available from the corresponding author
on request. Raw FASTQ files, per-sample SEQC count matrices and per-nucleus cell type
annotations are in GEO [GSE353157](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE353157).
