# 02_bulk_deseq2_fgsea.R
#
# Author : Jinrui (Louis) Liu (Computational and Systems Biology Program, MSKCC).
#          Contact: liuj11@mskcc.org
#
# Differential expression with an Ifitm3-KO x 5xFAD interaction model (DESeq2),
# followed by gene set enrichment analysis (fgsea) for each contrast.
#
# Design: ~ strain + condition + strain:condition
#   strain    : WT (reference) / 5xFAD
#   condition : noKO (reference) / KO (Ifitm3 knockout)
# Coefficients: Intercept, strain 5xFAD, condition KO, strain5xFAD:conditionKO
#
#   contrast      vector       meaning                                   used in
#   KO_WT         c(0,0,1,0)   KO effect in WT background                Fig. 2h, 2j
#   5x_KO_5x      c(0,0,1,1)   KO effect in 5xFAD background             Fig. 2h, 2i
#   interaction   c(0,0,0,1)   KO x 5xFAD interaction term               Fig. 2i, 2j, 2k
#   5x_WT         c(0,1,0,0)   5xFAD effect in Ifitm3-WT background      Extended Data Fig. 2g
#   5x_KO_KO      c(0,1,0,1)   5xFAD effect in Ifitm3-KO background      -
#
# Input : output/bulk_counts_N15.rds, output/bulk_metadata_N15.rds
#         data/gene_sets/m2.cp.wikipathways.v2024.1.Mm.symbols.gmt   (MSigDB mouse, WikiPathways)
#         data/gene_sets/m5.go.v2024.1.Mm.symbols.gmt                (MSigDB mouse, GO)
# Output: output/DEres_<contrast>.csv                 all tested genes
#         output/DEres.sig_<contrast>.csv             padj < 0.05 and |log2FC| > log2(1.5)
#         output/fgsea_<wiki|go>_<contrast>.rds       full fgsea results
#         output/fgsea_<wiki|go>_<contrast>_sig.rds   padj < 0.05 (also as .csv)
# The files used for the figures are also in results/ exactly as submitted:
#   DEres_KO_WT.csv, DEres_5x_KO_5x.csv, DEres_interaction.csv,
#   fgsea_wiki_interaction_sig.rds, fgsea_go_5x_WT_sig.rds
#
# Run from the repository root:  Rscript scripts/02_bulk_deseq2_fgsea.R

suppressPackageStartupMessages({
  library(DESeq2)
  library(fgsea)
})

data_dir  <- Sys.getenv("IFITM3_DATA_DIR", "data")
out_dir   <- Sys.getenv("IFITM3_OUT_DIR", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff   <- 0.05
log2FC_cutoff <- log2(1.5)
fgsea_seed    <- 42      # fgsea p-values are estimated by sampling; fix the seed

# ---- DESeq2 ----------------------------------------------------------------
counts <- readRDS(file.path(out_dir, "bulk_counts_N15.rds"))
meta   <- readRDS(file.path(out_dir, "bulk_metadata_N15.rds"))
stopifnot(identical(colnames(counts), rownames(meta)))

dds <- DESeqDataSetFromMatrix(counts, meta, design = ~ strain + condition + strain:condition)
dds <- DESeq(dds)
print(resultsNames(dds))
stopifnot(length(resultsNames(dds)) == 4)

contrasts <- list(
  KO_WT       = c(0, 0, 1, 0),
  `5x_KO_5x`  = c(0, 0, 1, 1),
  interaction = c(0, 0, 0, 1),
  `5x_WT`     = c(0, 1, 0, 0),
  `5x_KO_KO`  = c(0, 1, 0, 1)
)

make_res_df <- function(res, pcut = padj_cutoff, lfc_cut = log2FC_cutoff) {
  res <- as.data.frame(res)
  res <- cbind(gene = rownames(res), res)
  res$regulation <- "NS"
  res$regulation[res$padj < pcut & res$log2FoldChange >  lfc_cut] <- "Up"
  res$regulation[res$padj < pcut & res$log2FoldChange < -lfc_cut] <- "Down"
  # drop genes without a p value (all-zero / outliers) or adjusted p value (independent filtering)
  res <- res[!is.na(res$pvalue) & !is.na(res$padj), ]
  res[order(res$padj, abs(res$log2FoldChange), decreasing = c(FALSE, TRUE)), ]
}

de_results <- lapply(contrasts, function(cv) make_res_df(results(dds, contrast = cv)))

for (nm in names(de_results)) {
  res <- de_results[[nm]]
  write.csv(res, file.path(out_dir, paste0("DEres_", nm, ".csv")), row.names = FALSE)
  write.csv(res[res$regulation != "NS", ],
            file.path(out_dir, paste0("DEres.sig_", nm, ".csv")), row.names = FALSE)
  cat(sprintf("%-12s tested: %5d  up: %4d  down: %4d\n", nm, nrow(res),
              sum(res$regulation == "Up"), sum(res$regulation == "Down")))
}

# ---- fgsea -----------------------------------------------------------------
gene_set_files <- c(
  wiki = file.path(data_dir, "gene_sets", "m2.cp.wikipathways.v2024.1.Mm.symbols.gmt"),
  go   = file.path(data_dir, "gene_sets", "m5.go.v2024.1.Mm.symbols.gmt")
)

run_fgsea <- function(res_df, gmt_file, set_name, contrast_name) {
  # rank all tested genes by their DESeq2 log2 fold change
  ranks    <- setNames(res_df$log2FoldChange, res_df$gene)
  pathways <- gmtPathways(gmt_file)
  set.seed(fgsea_seed)
  fg <- as.data.frame(fgsea(pathways, ranks, minSize = 15, maxSize = 500))
  fg_sig <- fg[which(fg$padj < 0.05), ]

  prefix <- file.path(out_dir, paste0("fgsea_", set_name, "_", contrast_name))
  saveRDS(fg,     paste0(prefix, ".rds"))
  saveRDS(fg_sig, paste0(prefix, "_sig.rds"))
  write.csv(fg_sig[, c("pathway", "pval", "padj", "log2err", "ES", "NES", "size")],
            paste0(prefix, "_sig.csv"), row.names = FALSE)
  cat(sprintf("fgsea %-4s %-12s significant pathways: %d\n", set_name, contrast_name, nrow(fg_sig)))
  invisible(fg)
}

for (set_name in names(gene_set_files)) {
  for (nm in names(de_results)) {
    run_fgsea(de_results[[nm]], gene_set_files[[set_name]], set_name, nm)
  }
}

writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_02.txt"))
