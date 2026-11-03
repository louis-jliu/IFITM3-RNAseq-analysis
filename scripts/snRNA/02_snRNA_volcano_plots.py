"""02_snRNA_volcano_plots.py

Volcano plots of MAST differential expression, 5xFAD;Ifitm3-/- vs 5xFAD:
  Fig. 3c  microglia
  Fig. 4g  oligodendrocytes
Genes with |log2FC| > 0.1 and -log10(FDR) > 1.3 are coloured; selected genes are labelled.

Author  : Tanya Jain (original analysis, Jupyter notebook "Fig 3c_4g volcano plots")
Script  : converted from the notebook by Jinrui (Louis) Liu. Plot code is kept as run;
          only file paths were made relative and plots are saved to output/snRNA/figures/.
Contact : Yue-Ming Li (liy2@mskcc.org); Louis Liu (liuj11@mskcc.org) for repository issues.

Input   : output/snRNA/DE/microglia_ALL5x_3KO vs 5x.csv   (from 01_snRNA_clustering_annotation_DE.py)
          output/snRNA/DE/oligo5x_3KO vs 5x.csv
          MAST output tables with columns primerid, coef, fdr.
Requires: pandas, matplotlib, and Tanya Jain's plotting helper workshop_extra.py
          (volcano_plot), not included here; place it in scripts/snRNA/.

Run from the repository root:  python scripts/snRNA/02_snRNA_volcano_plots.py
"""
import os
import sys

import matplotlib.pyplot as plt
import pandas as pd

sys.path.append(os.path.join("scripts", "snRNA"))
import workshop_extra  # noqa: E402

out_dir = os.path.join(os.environ.get("IFITM3_OUT_DIR", "output"), "snRNA")
de_dir = os.path.join(out_dir, "DE")
fig_dir = os.path.join(out_dir, "figures")
os.makedirs(fig_dir, exist_ok=True)


def volcano(csv, genes_highlight, out_png):
    df = pd.read_csv(os.path.join(de_dir, csv), index_col=0)
    df_new = pd.DataFrame({'log2foldchange': df['coef'].values, 'pvalues': df['fdr'].values}, index=df['primerid'])

    # logfold_neg_thres = -0.1, logfold_pos_thres = 0.1 (all genes > 1.07 fold or < 0.9 fold),
    # pval_thres = 1.3 (-log10 of the FDR column)
    f = workshop_extra.volcano_plot(diff_expr_result=df_new,
                                    logfoldchange_label='log2foldchange',
                                    pvalue_label='pvalues',
                                    logfold_neg_thres=-0.1,
                                    logfold_pos_thres=0.1,
                                    pval_thres=1.3,
                                    s=10,
                                    genes_highlight=genes_highlight,
                                    genes_highlight_fontdict={'fontsize': 10, 'fontstyle': 'italic'},
                                    fig_size=(3, 3))
    f.legend(markerscale=1)
    f.axes.set_xlabel('log2(foldchange)', fontsize=12)
    f.axes.set_ylabel('-log10(pvalue)', fontsize=12)
    f.tick_params(labelsize=10)
    f.set_xlim([-2, 2])
    plt.savefig(os.path.join(fig_dir, out_png), bbox_inches="tight", dpi=300)
    plt.close("all")


# Fig. 3c
volcano("microglia_ALL5x_3KO vs 5x.csv",
        ['TREM2', 'TYROBP', 'CTSD', 'FTL1', 'C1QB', 'FTL1-PS1', 'MALAT1', 'B2M'],
        "Fig3c_volcano_microglia.png")

# Fig. 4g
volcano("oligo5x_3KO vs 5x.csv",
        ['APOD', 'IL33', 'MBP', 'MOBP', 'CNP', 'FTL1', 'FTH1', 'PTGDS', 'FTL1-PS1'],
        "Fig4g_volcano_oligodendrocytes.png")
