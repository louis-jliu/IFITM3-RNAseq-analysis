# 01_prepare_bulk_counts.R
#
# Author : Jinrui (Louis) Liu, with Kushal Dey and Christina Leslie (Computational and Systems
#          Biology Program, MSKCC). Contact: liuj11@mskcc.org
#
# Prepare the mouse prefrontal cortex (PFC) bulk RNA-seq count matrix for DESeq2.
# 15 samples: WT (n = 4), Ifitm3 KO (n = 4), 5xFAD (n = 4), 5xFAD;Ifitm3 KO (n = 3).
#
# Input : data/bulk/countData.csv        gene x sample count matrix (column 1 = GeneSymbol)
#         data/bulk/sample_metadata.csv  sample sheet (id, sex, strain, condition, genotype)
# Output: output/bulk_counts_N15.rds
#         output/bulk_metadata_N15.rds
#
# Run from the repository root:  Rscript scripts/01_prepare_bulk_counts.R

data_dir <- Sys.getenv("IFITM3_DATA_DIR", "data")
out_dir  <- Sys.getenv("IFITM3_OUT_DIR", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- sample sheet ----------------------------------------------------------
meta <- read.csv(file.path(data_dir, "bulk", "sample_metadata.csv"), stringsAsFactors = FALSE)
rownames(meta) <- meta$id
meta$strain    <- factor(meta$strain,    levels = c("WT", "5xFAD"))   # reference = WT
meta$condition <- factor(meta$condition, levels = c("noKO", "KO"))    # reference = Ifitm3 WT
meta$genotype  <- factor(meta$genotype,  levels = c("WT", "KO", "5x", "5x_KO"))

# ---- count matrix ----------------------------------------------------------
counts_df <- read.csv(file.path(data_dir, "bulk", "countData.csv"), check.names = FALSE)
# sample columns are named like "s_G_1__WT"; drop the genotype suffix
colnames(counts_df)[-1] <- sub("__[^_]+$", "", colnames(counts_df)[-1])

# collapse duplicated gene symbols, keeping the row with the highest total count
counts_df <- counts_df[order(rowSums(counts_df[, -1]), decreasing = TRUE), ]
counts_df <- counts_df[!duplicated(counts_df$GeneSymbol), ]
stopifnot(!any(duplicated(counts_df$GeneSymbol)))

# remove mitochondrial ("mt-") and RIKEN cDNA ("Rik") genes
counts_df <- counts_df[!grepl("mt-|Rik", counts_df$GeneSymbol), ]

counts <- as.matrix(counts_df[, -1])
rownames(counts) <- counts_df$GeneSymbol

# remove genes with <= 10 counts summed across samples; round to integers for DESeq2
counts <- round(counts[rowSums(counts) > 10, ])

# put samples in the same order as the sample sheet
stopifnot(setequal(colnames(counts), rownames(meta)))
counts <- counts[, rownames(meta)]

cat("genes:", nrow(counts), " samples:", ncol(counts), "\n")
print(table(meta$genotype))

saveRDS(counts, file.path(out_dir, "bulk_counts_N15.rds"))
saveRDS(meta,   file.path(out_dir, "bulk_metadata_N15.rds"))
