# 04_microglia_program_meta_analysis.R
#
# Author : Jinrui (Louis) Liu, with Kushal Dey and Christina Leslie (Computational and Systems
#          Biology Program, MSKCC). Contact: liuj11@mskcc.org
#
# Human microglia (ROSMAP snRNA-seq, Mathys et al. 2019): compare IFITM3+ vs IFITM3- microglia
# for curated microglial gene programs, separately in each pathology stage (Fig. 5e).
#
#  1. Microglia are IFITM3+ if log-normalised IFITM3 expression > 0, otherwise IFITM3-.
#  2. For every gene in a program and every stage, using only cells that express the gene (> 0):
#       beta = mean(IFITM3+) - mean(IFITM3-)   (log-normalised scale)
#       se   = sqrt(var+/n+ + var-/n-)
#     Genes with < 3 expressing cells in either group are skipped.
#  3. Per program and stage, gene-level betas are combined with a random-effects meta-analysis
#     (rmeta::meta.summaries, method = "random"): tau, its SE, z = tau / SE, two-sided p value.
#  4. Benjamini-Hochberg correction across programs within each stage.
#  tau < 0 means the program has lower expression in IFITM3+ than in IFITM3- microglia.
#
# Input : output/rosmap_microglia_log_normalized.tsv, output/rosmap_microglia_metadata.tsv
#           (from 03_extract_rosmap_microglia.py)
#         data/gene_programs/microglia_programs_Gan2024.tsv
# Output: output/microglia_program_meta_analysis_results.txt   (Fig. 5e; same as
#           results/microglia_program_meta_analysis_results.txt for the submitted figure)
#         output/gene_effects_<program>.tsv                    per-gene effect sizes
#
# Run from the repository root:  Rscript scripts/04_microglia_program_meta_analysis.R

suppressPackageStartupMessages({
  library(rmeta)
  library(dplyr)
})

data_dir <- Sys.getenv("IFITM3_DATA_DIR", "data")
out_dir  <- Sys.getenv("IFITM3_OUT_DIR", "output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stages      <- c("no-pathology", "early-pathology", "late-pathology")
min_samples <- 3

# ---- gene programs (see scripts/00_build_gene_program_table.py) -----------
gp <- read.delim(file.path(data_dir, "gene_programs", "microglia_programs_Gan2024.tsv"),
                 stringsAsFactors = FALSE)
program_order <- c("ATM", "IRM", "TGFbeta-dependent", "DAM-program", "LRRC33-dependent")
programs <- lapply(setNames(program_order, program_order),
                   function(p) unique(gp$gene[gp$program == p]))

# ---- expression ------------------------------------------------------------
expr <- read.delim(file.path(out_dir, "rosmap_microglia_log_normalized.tsv"),
                   row.names = 1, check.names = FALSE)
meta <- read.delim(file.path(out_dir, "rosmap_microglia_metadata.tsv"))
expr <- expr[, colnames(expr) %in% meta$TAG]
expr <- expr[rowSums(expr) > 0, ]
stopifnot(identical(colnames(expr), meta$TAG))

ifitm3 <- as.numeric(expr["IFITM3", ])
meta$IFITM3_Group <- ifelse(ifitm3 > 0, "IFITM3+", "IFITM3-")
print(table(meta$DiseaseStatus, meta$IFITM3_Group))

gene_effect <- function(neg, pos) {
  if (length(neg) < min_samples || length(pos) < min_samples) {
    return(c(beta = NA, se = NA, pval = NA))
  }
  c(beta = mean(pos) - mean(neg),
    se   = sqrt(var(neg) / length(neg) + var(pos) / length(pos)),
    pval = wilcox.test(neg, pos)$p.value)   # per-gene Wilcoxon test, reported only
}

results <- list()
for (program in names(programs)) {
  genes <- intersect(programs[[program]], rownames(expr))
  cat("Program:", program, "-", length(genes), "genes in data\n")

  gene_stats <- do.call(rbind, lapply(genes, function(g) {
    x <- as.numeric(expr[g, ])
    expressing <- x > 0
    do.call(rbind, lapply(stages, function(s) {
      in_stage <- expressing & meta$DiseaseStatus == s
      st <- gene_effect(neg = x[in_stage & meta$IFITM3_Group == "IFITM3-"],
                        pos = x[in_stage & meta$IFITM3_Group == "IFITM3+"])
      data.frame(gene = g, DiseaseStatus = s, beta = st[["beta"]], se = st[["se"]], pval = st[["pval"]])
    }))
  }))
  write.table(gene_stats, file.path(out_dir, paste0("gene_effects_", program, ".tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)

  for (s in stages) {
    d <- gene_stats[gene_stats$DiseaseStatus == s & !is.na(gene_stats$beta) &
                      !is.na(gene_stats$se) & !is.na(gene_stats$pval), ]
    m <- meta.summaries(d$beta, d$se, method = "random")
    z <- m$summary / m$se.summary
    results[[length(results) + 1]] <- data.frame(
      Program = program, DiseaseStatus = s,
      tau = m$summary, tau_se = m$se.summary, z = z, p = 2 * pnorm(-abs(z)),
      n_genes = nrow(d))
  }
}

df_meta <- bind_rows(results) %>%
  filter(!is.na(tau)) %>%
  group_by(DiseaseStatus) %>%
  mutate(p_adj = p.adjust(p, method = "BH")) %>%
  ungroup()

write.table(df_meta, file.path(out_dir, "microglia_program_meta_analysis_results.txt"),
            sep = "\t", row.names = FALSE, quote = FALSE)
print(as.data.frame(df_meta))
writeLines(capture.output(sessionInfo()), file.path(out_dir, "sessionInfo_04.txt"))
