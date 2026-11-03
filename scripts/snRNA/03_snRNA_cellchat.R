# 03_snRNA_cellchat.R
#
# CellChat cell-cell communication analysis ("Secreted Signaling" ligand-receptor pairs)
# on the mouse PFC snRNA-seq data, run separately per genotype (5xFAD, 5xFAD;Ifitm3-/-).
# Figure: Fig. 4a (chord diagram of CSF signaling, Pericytes/SM cells -> Microglia).
#
# Author : Tanya Jain (original analysis, R Markdown "Fig 4_Cell Chat", 2024-11-11)
# Script : converted from the knitted R Markdown by Jinrui (Louis) Liu. Code is kept as run;
#          file paths were made relative, and the genotype and signaling pathway are
#          command-line arguments (the notebook was re-run with the object and pathway edited).
# Contact: Yue-Ming Li (liy2@mskcc.org); Louis Liu (liuj11@mskcc.org) for repository issues.
#
# Input  : output/snRNA/adata_5x.h5ad, output/snRNA/adata_3KO.h5ad
#          per-genotype subsets of the annotated AnnData (obs["cluster_labels"]), written with
#            adata[adata.obs['batch_id'] == '5x_3KO'].write_h5ad('output/snRNA/adata_3KO.h5ad')
#          Gene symbols are upper case, so the human CellChatDB is used.
# Output : output/snRNA/cellchat/<genotype>_secreted_<pathway>.pdf, <genotype>_cellchat.rds
# Requires: CellChat (v1), Seurat, anndata (R), patchwork, NMF, ggalluvial
#
# Run from the repository root:
#   Rscript scripts/snRNA/03_snRNA_cellchat.R 3KO CSF
#   Rscript scripts/snRNA/03_snRNA_cellchat.R 5x  CSF

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(anndata)
  library(patchwork)
})
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
geno <- if (length(args) >= 1) args[1] else "3KO"
pathways.show <- if (length(args) >= 2) args[2] else "CSF"
out_dir <- file.path(Sys.getenv("IFITM3_OUT_DIR", "output"), "snRNA")
cc_dir <- file.path(out_dir, "cellchat")
dir.create(cc_dir, recursive = TRUE, showWarnings = FALSE)

# Load data
print("Convert from Scanpy to Seurat...")
data5x <- read_h5ad(file.path(out_dir, paste0("adata_", geno, ".h5ad")))
data5x <- CreateSeuratObject(counts = t(data5x$X), meta.data = data5x$obs)
print(str(data5x))

# Create a CellChat object
cellchat <- createCellChat(object = data5x, group.by = "cluster_labels")
cellchat <- setIdent(cellchat, ident.use = "cluster_labels")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  # use Secreted Signaling
cellchat@DB <- CellChatDB.use

# Preprocessing the expression data for cell-cell communication analysis
cellchat <- subsetData(cellchat)  # necessary even if using the whole database
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat, slot.name = "netP")
cellchat <- computeCommunProbPathway(cellchat)

pdf(file = file.path(cc_dir, paste0(geno, "_secreted_", pathways.show, ".pdf")), width = 11, height = 9)

# Aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge = F,
                 title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge = F,
                 title.name = "Interaction weights/strength")

# Signaling sent from each cell group
mat <- cellchat@net$weight
par(mfrow = c(3, 4), xpd = TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat),
                   title.name = rownames(mat)[i])
}

# All signaling pathways showing significant communications
print(cellchat@netP$pathways)

# Selected pathway: hierarchy, circle, chord (Fig. 4a) and heatmap
vertex.receiver <- seq(1, 5)
netVisual_aggregate(cellchat, signaling = pathways.show, vertex.receiver = vertex.receiver)
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
par(mfrow = c(1, 1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
par(mfrow = c(1, 1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Contribution of each ligand-receptor pair
netAnalysis_contribution(cellchat, signaling = pathways.show)
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
LR.show <- pairLR[1, ]  # show one ligand-receptor pair
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, vertex.receiver = vertex.receiver)
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")

# All significant pathways: network plots and L-R contributions
pathways.show.all <- cellchat@netP$pathways
levels(cellchat@idents)
vertex.receiver <- seq(1, 10)
for (i in seq_along(pathways.show.all)) {
  netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
  gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
  ggsave(filename = file.path(cc_dir, paste0(geno, "_", pathways.show.all[i], "_L-R_contribution.pdf")),
         plot = gg, width = 3, height = 2, units = "in", dpi = 300)
}

# Bubble and chord plots of significant L-R pairs
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:11), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = 3, targets.use = c(1:11), signaling = pathways.show, remove.isolate = FALSE)
pairLR.use <- extractEnrichedLR(cellchat, signaling = pathways.show)
netVisual_bubble(cellchat, sources.use = c(3, 4), targets.use = c(1:11), pairLR.use = pairLR.use, remove.isolate = TRUE)
netVisual_chord_gene(cellchat, sources.use = 4, targets.use = c(1:11), lab.cex = 0.5, legend.pos.y = 30)
netVisual_chord_gene(cellchat, sources.use = c(1, 2, 3, 4), targets.use = 8, legend.pos.x = 15)
netVisual_chord_gene(cellchat, sources.use = c(1, 2, 3, 4), targets.use = c(5:11), signaling = pathways.show,
                     legend.pos.x = 8)
netVisual_chord_gene(cellchat, sources.use = c(1, 2, 3, 4), targets.use = c(5:11), slot.name = "netP",
                     legend.pos.x = 10)
plotGeneExpression(cellchat, signaling = pathways.show)
plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = TRUE)

# Signaling roles (network centrality)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = pathways.show)
print(gg1 + gg2)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
print(ht1 + ht2)

# Communication patterns
library(NMF)
library(ggalluvial)
selectK(cellchat, pattern = "outgoing")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = 3)
netAnalysis_river(cellchat, pattern = "outgoing")
netAnalysis_dot(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = 4)
netAnalysis_river(cellchat, pattern = "incoming")
netAnalysis_dot(cellchat, pattern = "incoming")

# Functional and structural similarity of signaling networks
cellchat <- computeNetSimilarity(cellchat, type = "functional")
cellchat <- netEmbedding(cellchat, type = "functional")
cellchat <- netClustering(cellchat, type = "functional")
netVisual_embedding(cellchat, type = "functional", label.size = 3.5)
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
cellchat <- netClustering(cellchat, type = "structural")
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)

dev.off()
saveRDS(cellchat, file.path(cc_dir, paste0(geno, "_cellchat.rds")))
