"""01_snRNA_clustering_annotation_DE.py

Mouse prefrontal cortex single-nucleus RNA-seq (5xFAD vs 5xFAD;Ifitm3-/-, 15 months):
clustering, cell type annotation, marker/IFITM3/IL33/IL34 plots and MAST differential
expression in microglia and oligodendrocytes.

Figures : Fig. 3a, 3b, 3d; Fig. 4b, 4c, 4d, 4e, 4f, 4h, 4i; Extended Data Figs. 3-4
          (marker and IFITM3 expression panels).
          DE tables for Fig. 3c and 4g (plotted by 02_snRNA_volcano_plots.py).

Author  : Tanya Jain (original analysis, Jupyter notebook "Fig 3_4_Extended Data Fig 3_4")
Script  : converted from the notebook by Jinrui (Louis) Liu. Code cells are kept as run;
          only file paths were made relative, notebook-only lines (bare `adata` displays,
          %R magics) were translated, and plots are saved to output/snRNA/figures/.
Contact : Yue-Ming Li (liy2@mskcc.org); Louis Liu (liuj11@mskcc.org) for repository issues.

Input   : data/snRNA/adata_snRNA_PFC.h5ad
          Processed AnnData (24,025 nuclei; raw counts in layers["raw_counts"], log2
          library-size-normalised expression in X, per-gene z-scores in layers["zs_norm_log"],
          obs["batch_id"] = "5x" / "5x_3KO"). Produced from the SEQC count matrices
          (GEO GSE353157) after removing mitochondrial/ribosomal genes, library-size
          normalisation, log transformation and Scrublet doublet removal (see Methods).
          Not included in the repository because of its size (7.4 GB); available from the
          corresponding author on request.
Requires: scanpy, phenograph (sc.external.tl.phenograph), matplotlib, pandas;
          for the DE step rpy2 + anndata2ri and R packages scuttle, scran, MAST, data.table,
          plus Tanya Jain's helper run_MAST.r (pairwise_de), not included here.

Run from the repository root:  python scripts/snRNA/01_snRNA_clustering_annotation_DE.py
"""
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
from matplotlib.pyplot import rc_context

data_dir = os.environ.get("IFITM3_DATA_DIR", "data")
out_dir = os.path.join(os.environ.get("IFITM3_OUT_DIR", "output"), "snRNA")
fig_dir = os.path.join(out_dir, "figures")
de_dir = os.path.join(out_dir, "DE")
os.makedirs(fig_dir, exist_ok=True)
os.makedirs(de_dir, exist_ok=True)


def save(name):
    plt.savefig(os.path.join(fig_dir, name), bbox_inches="tight", dpi=300)
    plt.close("all")


adata = sc.read_h5ad(os.path.join(data_dir, "snRNA", "adata_snRNA_PFC.h5ad"))
print(adata)

# ---------------------------------------------------------------- clustering
# Clustering using PhenoGraph
sc.external.tl.phenograph(adata, clustering_algo='leiden', k=30, jaccard=True, primary_metric='euclidean',
                          resolution_parameter=1)
adata.obsp['pheno_jaccard_ig'] = scipy.sparse.csr_matrix(adata.obsp['pheno_jaccard_ig'])

sc.pl.umap(adata, color=['pheno_leiden'], size=20, show=False)
save("umap_pheno_leiden.png")

# ---------------------------------------------------------------- annotation
genes = {'Lymphocytes': ['CD3G'],
         'Inhibitory neurons/Interneurons': ['GAD1', 'GAD2'],
         'Excitatory neurons': ['SLC17A7', 'SLC17A6'],
         'Pyramidal': ['CUX1', 'RORB'],
         'Motor neurons': ['ISL1', 'CHAT', 'SOD1', 'OLIG2'],
         'Endothelial cells': ['PECAM1', 'CDH5', 'CLDN5'],
         'Pericytes/Smooth muscle cells': ['VTN', 'PDGFRB'],
         'Oligodendrocytes': ['MBP'],
         'Microglia': ['HEXB', 'TMEM119'],
         'OPCs': ['PDGFRA'],
         'Astrocytes': ['AQP4'],
         'BAMs': ['MRC1', 'CD163']}

sc.set_figure_params(scanpy=True, fontsize=14)
ac = sc.pl.matrixplot(adata, genes, layer='zs_norm_log', groupby='pheno_leiden',
                      show=False, figsize=(10, 9), standard_scale='var')
ac['mainplot_ax'].set_xlabel('Genes')
ac['mainplot_ax'].set_ylabel('Clusters')
save("matrixplot_markers_by_cluster.png")

# Cluster numbers refer to the PhenoGraph run in the submitted analysis (stored in
# obs["pheno_leiden"] of the input file). A new PhenoGraph run may number clusters differently.
cluster_labels = {0: 'Excitatory neurons', 1: 'Excitatory neurons', 2: 'Oligodendrocytes', 3: 'Excitatory neurons',
                  4: 'Oligodendrocytes', 5: 'Microglia', 6: 'Astrocytes', 7: 'OPCs', 8: 'Microglia',
                  9: 'Inhibitory neurons', 10: 'Inhibitory neurons', 11: 'Inhibitory neurons',
                  12: 'Inhibitory neurons', 13: 'Inhibitory neurons', 14: 'Excitatory neurons', 15: 'Microglia',
                  16: 'Excitatory neurons', 17: 'Endothelial cells', 18: 'Lymphocytes', 19: 'Oligodendrocytes',
                  20: 'Pericytes/SM cells', 21: 'Endothelial cells', 22: 'Pericytes/SM cells', 23: 'Microglia',
                  24: 'Endothelial cells'}
adata.obs['cluster_labels'] = adata.obs['pheno_leiden'].map(cluster_labels)

# Fig. 3a
sc.pl.umap(adata, color=['cluster_labels'], size=10, show=False)
plt.gcf().set_size_inches(3, 3)
save("Fig3a_umap_cluster_labels.png")

genes = ['AQP4', 'CDH5', 'SLC17A7', 'GAD1', 'CD3G', 'HEXB', 'PDGFRA', 'MBP', 'VTN']
sc.pl.stacked_violin(adata, genes, groupby='cluster_labels', swap_axes=True, dendrogram=False,
                     standard_scale='var', figsize=(6, 4), show=False)
save("stacked_violin_markers.png")

adata.layers['scaled'] = sc.pp.scale(adata, copy=True).X

# ---------------------------------------------------------------- IFITM3 expression
adata_5x = adata[~adata.obs['batch_id'].isin(['5x_3KO'])]
adata_3KO = adata[~adata.obs['batch_id'].isin(['5x'])]

genes = ['IFITM3']
sc.set_figure_params(scanpy=True, fontsize=14)
fig = plt.figure(figsize=(3 * 2, 4 * 1.5))
ax1 = fig.add_subplot(1, 1, 1)
sc.pl.dotplot(adata_5x, genes, groupby='cluster_labels', log=True, layer='zs_norm_log',
              show=False, ax=ax1, swap_axes=True)
save("dotplot_IFITM3_5x.png")
sc.pl.matrixplot(adata_5x, genes, 'cluster_labels', cmap='Blues', standard_scale='var',
                 colorbar_title='column scaled\nexpression', show=False)
save("matrixplot_IFITM3_5x_scaled.png")
sc.pl.matrixplot(adata_5x, genes, 'cluster_labels', dendrogram=False,
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r', show=False)
save("matrixplot_IFITM3_5x_zscore.png")

# ---------------------------------------------------------------- microglia states (Fig. 3b, 3d)
adata_microglia = adata[adata.obs['cluster_labels'] == 'Microglia', :]

genes = {'Homeostatic': ['TMEM119', 'HEXB', 'CX3CR1', 'P2RY12', 'CSF1R'],
         'DAM1': ['TYROBP', 'CTSB', 'CTSD', 'APOE', 'B2M', 'FTH1', 'LYZ2'],
         'DAM2': ['TREM2', 'AXL', 'CST7', 'CTSL', 'LPL', 'CD9', 'CSF1', 'CCL6', 'ITGAX', 'TIMP2']}
sc.set_figure_params(scanpy=True, fontsize=14)
fig = plt.figure(figsize=(5 * 2, 4 * 1.5))
ax1 = fig.add_subplot(1, 1, 1)
sc.pl.dotplot(adata_microglia, genes, groupby='batch_id', log=True, layer='zs_norm_log', show=False, ax=ax1)
save("dotplot_microglia_DAM.png")
sc.pl.matrixplot(adata_microglia, genes, 'batch_id', dendrogram=False,
                 colorbar_title='mean z-score', layer='zs_norm_log', vmin=-2, vmax=2, cmap='RdBu_r', show=False)
save("Fig3d_matrixplot_microglia_DAM_zscore.png")
sc.pl.matrixplot(adata_microglia, genes, 'batch_id', cmap='Blues', standard_scale='var',
                 colorbar_title='column scaled\nexpression', show=False)
save("matrixplot_microglia_DAM_scaled.png")
sc.set_figure_params(scanpy=True, fontsize=14)
sc.pl.matrixplot(adata_microglia, genes, layer='zs_norm_log', groupby='batch_id', show=False, cmap='RdBu_r')
save("matrixplot_microglia_DAM.png")

adata_peri = adata[adata.obs['cluster_labels'] == 'Pericytes/SM cells', :]
adata_oligo = adata[adata.obs['cluster_labels'] == 'Oligodendrocytes', :]
adata_NVU = adata[~adata.obs['cluster_labels'].isin(['Astrocytes', 'Excitatory neurons', 'Lymphocytes', 'Microglia',
                                                      'OPCs', 'Oligodendrocytes', 'Inhibitory neurons', 'BAMs'])]

# ---------------------------------------------------------------- MAST DE (Fig. 3c, 4g tables)
# Notebook %R magics translated to rpy2 calls. Requires run_MAST.r (Tanya Jain's helper,
# provides pairwise_de(df1, df2, name, folder)) in scripts/snRNA/.
import rpy2.robjects as ro
import anndata2ri

anndata2ri.activate()
ro.r('suppressPackageStartupMessages({library(scuttle); library(scran); library(MAST); library(data.table)})')
ro.r(f'source("{os.path.join("scripts", "snRNA", "run_MAST.r")}")')


def run_mast(adata_sub, name):
    # Since we want to compute differentially expressed genes comparing two genotypes,
    # we supply the cluster labels and genotype as part of sce_v2
    sce_v2 = sc.AnnData(X=adata_sub.X,
                        obs=pd.DataFrame({'phenograph_label': adata_sub.obs['cluster_labels'].astype('str'),
                                          'genotype': adata_sub.obs['batch_id'].astype('str')},
                                         index=adata_sub.obs.index),
                        var=pd.DataFrame(index=adata_sub.var.index))
    ro.globalenv['sce_v2'] = sce_v2
    ro.r('counts(sce_v2) <- assay(sce_v2)')
    ro.r("genoA <- which(colData(sce_v2)$genotype == '5x')")
    ro.r("genoB <- which(colData(sce_v2)$genotype == '5x_3KO')")
    ro.r('print(length(genoA)); print(length(genoB))')
    ro.r('df1 <- t(data.frame(counts(sce_v2)[, genoA]))')  # transpose because in sce genes are rows
    ro.r('df2 <- t(data.frame(counts(sce_v2)[, genoB]))')
    # Syntax: pairwise_de(dataframe1, dataframe2, 'output_filename', 'output_folder')
    ro.r(f"pairwise_de(df2, df1, '5x_3KO vs 5x', '{os.path.join(de_dir, name)}')")


run_mast(adata_microglia, 'microglia_ALL')  # -> microglia_ALL5x_3KO vs 5x.csv (Fig. 3c)
run_mast(adata_oligo, 'oligo')              # -> oligo5x_3KO vs 5x.csv (Fig. 4g); same call as microglia

# ---------------------------------------------------------------- IL34 / IL33 (Fig. 4)
adata_noNVU = adata[~adata.obs['cluster_labels'].isin(['Endothelial cells', 'Pericytes/SM cells'])]
adata_onlybraincells = adata[~adata.obs['cluster_labels'].isin(['Endothelial cells', 'Pericytes/SM cells', 'Lymphocytes'])]
adata_5X = adata[adata.obs['batch_id'] == '5x', :]
adata_3KO = adata[adata.obs['batch_id'] == '5x_3KO', :]
adata_5X_onlybraincells = adata_5X[~adata_5X.obs['cluster_labels'].isin(['Endothelial cells', 'Pericytes/SM cells', 'Lymphocytes'])]

# Fig. 4d
sc.pl.matrixplot(adata_onlybraincells, 'IL34', 'cluster_labels', dendrogram=False,
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r',
                 standard_scale='var', show=False)
save("Fig4d_right_matrixplot_IL34_braincells.png")
sc.pl.matrixplot(adata, 'IL34', 'cluster_labels', dendrogram=False,
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r',
                 standard_scale='var', show=False)
save("Fig4d_left_matrixplot_IL34_all.png")

sc.pl.matrixplot(adata_5X_onlybraincells, 'IFITM3', 'cluster_labels', dendrogram=False,
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r',
                 standard_scale='var', show=False)
save("matrixplot_IFITM3_5x_braincells.png")
sc.pl.matrixplot(adata_5X, 'IFITM3', 'cluster_labels', dendrogram=False,
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r',
                 standard_scale='var', show=False)
save("matrixplot_IFITM3_5x_all.png")

# Fig. 4e
fig = plt.figure(figsize=(3 * 1, 3 * 1))
ax1 = fig.add_subplot(1, 1, 1)
sc.pl.dotplot(adata_peri, 'IL34', groupby='batch_id', log=True, layer='zs_norm_log', show=False, ax=ax1, swap_axes=True)
save("Fig4e_dotplot_IL34_pericytes.png")
with rc_context({'figure.figsize': (3, 2.5)}):
    sc.pl.violin(adata_peri, 'IL34', groupby='batch_id', stripplot=False, standard_scale='var', show=False)
    save("violin_IL34_pericytes.png")

# Fig. 4h
fig = plt.figure(figsize=(3 * 1, 3 * 1))
ax1 = fig.add_subplot(1, 1, 1)
sc.pl.dotplot(adata_oligo, 'IL33', groupby='batch_id', log=True, layer='zs_norm_log', show=False, ax=ax1, swap_axes=True)
save("Fig4h_dotplot_IL33_oligodendrocytes.png")
with rc_context({'figure.figsize': (3, 2.5)}):
    sc.pl.violin(adata_oligo, 'IL33', groupby='batch_id', stripplot=False, standard_scale='var', show=False)
    save("violin_IL33_oligodendrocytes.png")

# Fig. 3b, 4c, 4f
sc.pl.umap(adata_microglia, color=['HEXB', 'TMEM119'], size=50, show=False)
plt.gcf().set_size_inches(6, 3)
save("Fig3b_umap_microglia_HEXB_TMEM119.png")
sc.pl.umap(adata_peri, color=['VTN'], size=50, show=False)
plt.gcf().set_size_inches(3, 3)
save("Fig4c_umap_pericytes_VTN.png")
sc.pl.umap(adata_oligo, color=['MBP'], size=50, show=False)
plt.gcf().set_size_inches(3, 3)
save("Fig4f_umap_oligodendrocytes_MBP.png")

# Fig. 4i
sc.pl.matrixplot(adata_3KO, 'IL33', 'cluster_labels', dendrogram=False,
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r',
                 standard_scale='var', show=False)
save("Fig4i_matrixplot_IL33_5x_3KO.png")
sc.pl.matrixplot(adata_5X, 'IL33', 'cluster_labels', dendrogram=False,
                 colorbar_title='mean z-score', layer='scaled', vmin=-2, vmax=2, cmap='RdBu_r', show=False)
save("matrixplot_IL33_5x.png")

# Fig. 4b and related violins
for g in ['TREM2', 'TYROBP', 'CSF1', 'IL33', 'IL34', 'CSF1R']:
    with rc_context({'figure.figsize': (9, 2)}):
        sc.pl.violin(adata_onlybraincells, g, groupby='cluster_labels', show=False)
        save(f"violin_{g}_by_celltype.png")
