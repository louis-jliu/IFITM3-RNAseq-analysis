# 05_make_figures.R
#
# Author : Jinrui (Louis) Liu, with Kushal Dey and Christina Leslie (Computational and Systems
#          Biology Program, MSKCC). Contact: liuj11@mskcc.org
#
# Recreates the submitted panels: Fig. 2h, 2i, 2j, 2k, Extended Data Fig. 2g and Fig. 5e.
# This is figures.R from the submission with file paths made relative to the repository
# (and a fixed seed for ggrepel label placement).
#
# By default it plots from results/, which holds the exact tables behind the submitted figures.
# To plot from a fresh run of scripts 01-04 instead:
#   IFITM3_FIG_INPUT=output Rscript scripts/05_make_figures.R
#
# Output: output/figures/ (same file names as the submitted panels)
#
# Run from the repository root:  Rscript scripts/05_make_figures.R

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

data_dir  <- Sys.getenv("IFITM3_DATA_DIR", "data")
input_dir <- Sys.getenv("IFITM3_FIG_INPUT", "results")
fig_dir   <- file.path(Sys.getenv("IFITM3_OUT_DIR", "output"), "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
cat("Plotting from:", input_dir, "\n")

res.KO.5xFAD = read.csv(file.path(input_dir, 'DEres_5x_KO_5x.csv'))
res.KO.WT = read.csv(file.path(input_dir, 'DEres_KO_WT.csv'))
res.interaction = read.csv(file.path(input_dir, 'DEres_interaction.csv'))

# cap |log2FC| at 5 for display (the Pearson r below is computed on the capped values)
.set_max_abs_log2FC = function(res, max_abs_log2FC = 5){
  signFC = sign(res$log2FoldChange)
  res$log2FoldChange[abs(res$log2FoldChange) > max_abs_log2FC] = max_abs_log2FC*signFC[abs(res$log2FoldChange) > max_abs_log2FC]
  return(res)
}
res.KO.WT = .set_max_abs_log2FC(res.KO.WT)
res.KO.5xFAD = .set_max_abs_log2FC(res.KO.5xFAD)
res.interaction = .set_max_abs_log2FC(res.interaction)

.cor_label = function(x, y){
  cor_value = round(cor(x, y, method='pearson'),3)
  cor_pvalue = cor.test(x, y, method='pearson')$p.value
  stars = if(cor_pvalue < 0.001) '***' else if(cor_pvalue < 0.01) '**' else if(cor_pvalue < 0.05) '*' else 'NS'
  cat(sprintf("  Pearson r = %.3f (p = %.3g)\n", cor_value, cor_pvalue))
  sprintf("Pearson r = %.3f %s", cor_value, stars)
}

# Figure 2h: 4-way plot of 5x_KO/5x vs. KO/WT log2FoldChange
df.KO.5xFAD_KO.WT = inner_join(res.KO.5xFAD[, c('gene','log2FoldChange','padj')],
           res.KO.WT[, c('gene','log2FoldChange','padj')],
           by = 'gene') %>%
  mutate(
    significance = factor(
      case_when(
        padj.x < 0.05 & padj.y < 0.05 ~ 'Both Significant',
        padj.x < 0.05 ~ '5x_KO/5x Significant',
        padj.y < 0.05 ~ 'KO/WT Significant',
        TRUE ~ 'Not Significant'
        ),
      levels = c('Both Significant', '5x_KO/5x Significant', 'KO/WT Significant', 'Not Significant')
      ),
    label = ifelse((significance != 'Not Significant' &
                     abs(log2FoldChange.x) > log2(1.5)) |
                     (significance != 'Not Significant' &
                        abs(log2FoldChange.y) > log2(1.5)), gene, NA)
  )
cat("Fig. 2h\n")
cor_label = .cor_label(df.KO.5xFAD_KO.WT$log2FoldChange.x, df.KO.5xFAD_KO.WT$log2FoldChange.y)

p3 <- ggplot(df.KO.5xFAD_KO.WT, aes(x = log2FoldChange.x, y = log2FoldChange.y, color = significance)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 'dashed', color = '#999999') +
  geom_hline(yintercept = c(-log2(1.5), log2(1.5)), linetype = 'dashed', color = '#999999') +
  geom_point(data = . %>% filter(significance == 'Not Significant'), alpha = 0.5) +
  geom_point(data = . %>% filter(significance == '5x_KO/5x Significant'), alpha = 1) +
  geom_point(data = . %>% filter(significance == 'KO/WT Significant'), alpha = 1) +
  geom_point(data = . %>% filter(significance == 'Both Significant'), alpha = 1) +
  geom_text_repel(aes(label = label),
                  max.overlaps = 10,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  force = 10, color = 'black',
                  fontface = 'italic',
                  show.legend = FALSE, seed = 42) +
  scale_color_manual(values = c(
    'Not Significant' = '#999999',
    '5x_KO/5x Significant' = '#009E73',
    'KO/WT Significant' = '#F0E442',
    'Both Significant' = '#D55E00'
  )) +
  theme_minimal() +
  labs(x = '5x_KO/5x log2FoldChange',
       y = 'KO/WT log2FoldChange',
       title = '5x_KO/5x vs. KO/WT Effects',
       color = 'Significance:') +
  theme(
    legend.position = 'bottom',
    legend.box = 'horizontal'
  ) +
  annotate("text", x = -5, y = 5,
           label = cor_label,
           hjust = 0, vjust = 1, size = 3.5)
ggsave(file.path(fig_dir, '2h-4way-5x_KO_5x_vs_KO_WT.pdf'), plot = p3,
       device = 'pdf',
       height = 6, width = 8)

# Figure 2i: 4-way plot of 5x_KO/5x vs. interaction log2FoldChange
df.KO.5xFAD_interaction = inner_join(res.KO.5xFAD[, c('gene','log2FoldChange','padj')],
           res.interaction[, c('gene','log2FoldChange','padj')],
           by = 'gene') %>%
  mutate(
    significance = factor(
      case_when(
        padj.x < 0.05 & padj.y < 0.05 ~ 'Both Significant',
        padj.x < 0.05 ~ '5x_KO/5x Significant',
        padj.y < 0.05 ~ 'Interaction Significant',
        TRUE ~ 'Not Significant'
      ),
      levels = c('Not Significant', 'Interaction Significant', '5x_KO/5x Significant', 'Both Significant')
    ),
    label = ifelse((significance != 'Not Significant' &
                      abs(log2FoldChange.x) > log2(1.5)) |
                     (significance != 'Not Significant' &
                        abs(log2FoldChange.y) > log2(1.5)), gene, NA)
  )
cat("Fig. 2i\n")
cor_label = .cor_label(df.KO.5xFAD_interaction$log2FoldChange.x, df.KO.5xFAD_interaction$log2FoldChange.y)

p1 <- ggplot(data = df.KO.5xFAD_interaction,aes(x = log2FoldChange.x, y = log2FoldChange.y, color = significance)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 'dashed', color = '#999999') +
  geom_hline(yintercept = c(-log2(1.5), log2(1.5)), linetype = 'dashed', color = '#999999') +
  geom_point(data = . %>% filter(significance == 'Not Significant'), alpha = 0.5) +
  geom_point(data = . %>% filter(significance == '5x_KO/5x Significant'), alpha = 1) +
  geom_point(data = . %>% filter(significance == 'Interaction Significant'), alpha = 1) +
  geom_point(data = . %>% filter(significance == 'Both Significant'), alpha = 1) +
  geom_text_repel(aes(label = label),
                  max.overlaps = 10,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  force = 10, color = 'black',
                  fontface = 'italic',
                  show.legend = FALSE, seed = 42) +
  scale_color_manual(values = c(
    'Not Significant' = '#999999',
    '5x_KO/5x Significant' = '#009E73',
    'Interaction Significant' = '#56B4E9',
    'Both Significant' = '#D55E00'
  )) +
  theme_minimal() +
  labs(x = '5x_KO/5x log2FoldChange',
       y = 'Interaction log2FoldChange',
       title = '5x_KO/5x vs. Interaction Effects',
       color = 'Significance:') +
  theme(
    legend.position = 'bottom',
    legend.box = 'horizontal'
  ) +
  annotate("text", x = -5, y = 5,
           label = cor_label,
           hjust = 0, vjust = 1, size = 3.5)
ggsave(file.path(fig_dir, '2i-4way-5xKO_5x_vs_interaction.pdf'), plot = p1,
       device = 'pdf',
       height = 6, width = 8)

# Figure 2j: 4-way plot of KO/WT vs. interaction log2FoldChange
df.KO.WT_interaction = inner_join(res.KO.WT[, c('gene','log2FoldChange','padj')],
           res.interaction[, c('gene','log2FoldChange','padj')],
           by = 'gene') %>%
  mutate(
    significance = factor(
      case_when(
        padj.x < 0.05 & padj.y < 0.05 ~ 'Both Significant',
        padj.x < 0.05 ~ 'KO/WT Significant',
        padj.y < 0.05 ~ 'Interaction Significant',
        TRUE ~ 'Not Significant'
      ),
      levels = c('Not Significant', 'Interaction Significant', 'KO/WT Significant', 'Both Significant')
    ),
    label = ifelse((significance != 'Not Significant' &
                      abs(log2FoldChange.x) > log2(1.5)) |
                     (significance != 'Not Significant' &
                        abs(log2FoldChange.y) > log2(1.5)), gene, NA)
  )
cat("Fig. 2j\n")
cor_label = .cor_label(df.KO.WT_interaction$log2FoldChange.x, df.KO.WT_interaction$log2FoldChange.y)

p2 <- ggplot(df.KO.WT_interaction , aes(x = log2FoldChange.x, y = log2FoldChange.y, color = significance)) +
  geom_vline(xintercept = c(-log2(1.5), log2(1.5)), linetype = 'dashed', color = '#999999') +
  geom_hline(yintercept = c(-log2(1.5), log2(1.5)), linetype = 'dashed', color = '#999999') +
  geom_point(data = . %>% filter(significance == 'Not Significant'), alpha = 0.5) +
  geom_point(data = . %>% filter(significance == 'KO/WT Significant'), alpha = 1) +
  geom_point(data = . %>% filter(significance == 'Interaction Significant'), alpha = 1) +
  geom_point(data = . %>% filter(significance == 'Both Significant'), alpha = 1) +
  geom_text_repel(aes(label = label),
                  max.overlaps = 10,
                  box.padding = 0.5,
                  point.padding = 0.5,
                  force = 10, color = 'black',
                  fontface = 'italic',
                  show.legend = FALSE, seed = 42) +
  scale_color_manual(values = c(
    'Not Significant' = '#999999',
    'KO/WT Significant' = '#F0E442',
    'Interaction Significant' = '#56B4E9',
    'Both Significant' = '#D55E00'
  )) +
  theme_minimal() +
  labs(x = 'KO/WT log2FoldChange',
       y = 'Interaction log2FoldChange',
       title = 'KO/WT vs. Interaction Effects',
       color = 'Significance:') +
  theme(
    legend.position = 'bottom',
    legend.box = 'horizontal'
  ) +
  annotate("text", x = 5, y = 3.7,
           label = cor_label,
           hjust = 1, vjust = 1, size = 3.5)
ggsave(file.path(fig_dir, '2j-4way-KO_WT_vs_interaction.pdf'), plot = p2,
       device = 'pdf',
       height = 6, width = 8)

plot_scatter = function(fgsea,title,name) {
  fgsea = fgsea[!is.na(fgsea$pathway),]
  p <- ggplot(fgsea, aes(x=NES, y=reorder(pathway, NES))) +
    geom_point(aes(size=size, color=padj)) +
    scale_color_gradient(low = "red", high = "blue") +
    labs(size = 'pathway size',color='p adjusted', y='', x='NES',
         title = title) +
    theme_bw() +
    theme(axis.text.x = element_text(size=10),
          axis.text.y = element_text(size=10),
          plot.title = element_text(size=12))
  ggsave(file.path(fig_dir, paste0(name, '.pdf')), plot = p, width = 10, height = 8)
  cat(sprintf("%s: %d pathways\n", name, nrow(fgsea)))
}

# Figure 2k: GSEA WikiPathways scatter for the interaction term
fgsea.interaction_wiki = readRDS(file.path(input_dir, 'fgsea_wiki_interaction_sig.rds'))
plot_scatter(fgsea.interaction_wiki, 'GSEA-wiki: interaction', '2k-scatter-wiki_interaction')

# Extended Data Figure 2g: GSEA GO scatter of synaptic/neurotransmitter pathways for 5x vs. WT
# (manually curated list of terms: data/EDFig2g_curated_GO_terms.txt)
fgsea.5xFAD_go = readRDS(file.path(input_dir, 'fgsea_go_5x_WT_sig.rds'))
go_neuronal.paths = readLines(file.path(data_dir, 'EDFig2g_curated_GO_terms.txt'))
go_neuronal.paths = go_neuronal.paths[nzchar(go_neuronal.paths) & !startsWith(go_neuronal.paths, '#')]
missing = setdiff(go_neuronal.paths, fgsea.5xFAD_go$pathway)
if (length(missing) > 0) warning('Curated GO terms not significant in this input: ', paste(missing, collapse = ', '))
plot_scatter(fgsea.5xFAD_go[fgsea.5xFAD_go$pathway %in% go_neuronal.paths,], 'GSEA-GO: 5x vs. WT', 'Ex_2g-scatter-GO_5x_WT')

# Figure 5e: meta-analysis effect size of IFITM3+ vs. IFITM3- microglia programs by disease stage
df_meta = read.delim(file.path(input_dir, 'microglia_program_meta_analysis_results.txt'))
dplot <- df_meta[, c('Program', 'DiseaseStatus', 'p_adj')]
dplot$neg_log10_p <- -log10(dplot$p_adj)
dplot$DiseaseStatus <- factor(dplot$DiseaseStatus,
                              levels = c('no-pathology', 'early-pathology', 'late-pathology'))

dot <- ggplot(dplot, aes(x=DiseaseStatus, y=Program)) +
  geom_point(aes(size=neg_log10_p, color=neg_log10_p)) +
  scale_color_gradient(low="white", high="red", name="-log10(p_adj)") +
  scale_size_continuous(range = c(1, 8), name="-log10(p_adj)") +
  labs(y='', x='',
       title='Effect size estimation: IFITM3+ vs.IFITM3-') +
  guides(size = guide_legend(reverse=TRUE)) +
  theme_bw() +
  theme(
    legend.key=element_blank(),
    plot.title = element_text(size=12, face='bold'),
    axis.text.x = element_text(angle=45, hjust=1, size=12, face='bold'),
    axis.text.y = element_text(size=10, face='bold'),
    legend.text = element_text(size=10),
    legend.title = element_text(size=12),
    legend.key.size = unit(1.2, "cm"))
ggsave(plot = dot, device = "pdf",
       filename = file.path(fig_dir, '5e-kellis_microglia_meta_analysis_scatter_p_adj_V3.pdf'),
       width = 8, height = 8)
cat("Figures written to", fig_dir, "\n")
