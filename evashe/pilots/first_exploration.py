# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/02/21
content:    First exploration. For the sake of simplicity, let's try following
            scanpy's tutorial and stray as little as possible.
'''
import os
import sys
import argparse
import numpy as np
import pandas as pd

import anndata
import scanpy as sc

import matplotlib as mpl
mpl.rcParams['font.size'] = 7
import matplotlib.pyplot as plt
import seaborn as sns


if __name__ == '__main__':

    print('Load data')
    data_fdn = '../../data/sequencing'
    fn_normalised = f'{data_fdn}/normalised.h5ad'
    adata = anndata.read_h5ad(fn_normalised)
    print('Loaded')

    #sc.pl.highest_expr_genes(adata, n_top=20)

    # No need to filter, the minimum is already 203
    #sc.pp.filter_cells(adata, min_genes=200)

    # So we skip this one too
    #sc.pp.filter_genes(adata, min_cells=3)

    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True,
        )
    #sc.pl.violin(
    #    adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
    #    jitter=0.4, multi_panel=True,
    #    )

    # Let's ignore filtering for now
    #adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    #adata = adata[adata.obs.pct_counts_mt < 5, :]

    # It's notmalized already
    #sc.pp.normalize_total(adata, target_sum=1e4)

    print('Log')
    sc.pp.log1p(adata)

    print('Feature selection')
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
    #sc.pl.highly_variable_genes(adata)
    adataf = adata[:, adata.var.highly_variable]

    # Some black magic with dynamic ranges
    sc.pp.regress_out(adataf, ['total_counts', 'pct_counts_mt'])
    sc.pp.scale(adataf, max_value=10)

    print('PCA')
    sc.tl.pca(adataf, svd_solver='arpack')
    #sc.pl.pca(adataf, color='Runx1')
    #sc.pl.pca_variance_ratio(adataf, log=True)

    print('Similarity graph')
    sc.pp.neighbors(adataf, n_neighbors=10, n_pcs=10)

    #print('Embedding')
    #sc.tl.umap(adataf)

    print('Clustering')
    sc.tl.leiden(adataf, resolution=0.1)

    sys.exit()

    print('Move back to full adata')
    adata.obs['leiden'] = adataf.obs['leiden']
    adata.obsm = adataf.obsm

    #sc.pl.umap(adata, color=['leiden', 'Runx1', 'Wnt2', 'Wnt5a', 'Notch1'])

    # Sort clusters by Runx1 expression
    df = adata.obs[['leiden']].copy()
    df['Runx1'] = adata[:, ['Runx1']].X
    cl_order = list(df.groupby('leiden').mean().sort_values('Runx1').index)
    new_cluster_names = ['Runx1_q'+str(cl_order.index(x) + 1) for x in adataf.obs['leiden'].cat.categories]
    adata.rename_categories('leiden', new_cluster_names)

    genes = ['Runx1', 'Ptprc', 'Cdh5', 'Pecam1', 'Sox17', 'Tal1', 'Gfi1', 'Gfi1b',
             'Gata2', 'Erg', 'Fli1', 'Lyl1', 'Lmo2']
    for x in ['Bmp', 'Notch', 'Wnt']:
        genes += list(adata.var_names[adata.var_names.str.startswith(x)])
    cat_order = ['Runx1_q'+str(i+1) for i in range(len(new_cluster_names))]
    sc.pl.dotplot(adata, genes, groupby='leiden', categories_order=cat_order,
            swap_axes=True)
    plt.subplots_adjust(bottom=0.2, top=1.0, left=0.3)

    genes_umap = ['Runx1', 'Ptprc', 'Cdh5', 'Bmp2k', 'Notch1', 'Bmpr1a', 'Notch4']
    sc.pl.umap(adata, color=['leiden'] + genes_umap)

