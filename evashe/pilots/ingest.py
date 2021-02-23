# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/02/21
content:    First, ingest the data into an h5ad file.
'''
import os
import sys
import argparse
import numpy as np
import pandas as pd

import anndata


if __name__ == '__main__':

    # Hash data is basically useless, so use only counts
    data_fdn = '../../data/sequencing'
    counts_fn = f'{data_fdn}/GSE163757_scRNAseq_9562cells_GE-rawcounts.txt'
    
    # Genes are rows
    print('Load raw counts from TSV')
    counts = pd.read_csv(counts_fn, sep='\t', index_col=0)


    # Convert gene ids into names
    conv_fn = '../../data/genome/mart_export.txt'
    conv = pd.read_csv(conv_fn, sep='\t')
    inter = set(counts.index) & set(conv['Gene stable ID'])
    missing = set(counts.index) - set(conv['Gene stable ID'])
    conv = conv.loc[conv['Gene stable ID'].isin(inter)]
    genes = conv['Gene name'].value_counts()
    multiplets = genes[genes > 1].index
    singlets = genes[genes == 1].index
    genes = np.sort(np.concatenate([genes.index, list(missing)]))
    conv_singl = conv.loc[conv['Gene name'].isin(singlets)].set_index('Gene name')['Gene stable ID']
    conv_multi = conv.loc[conv['Gene name'].isin(multiplets)].groupby('Gene name').apply(lambda x: x['Gene stable ID'].values)
    # Transpose too
    countsg = np.zeros((counts.shape[1], len(genes)), np.float32)
    for ig, gene in enumerate(genes):
        if gene in singlets:
            countsg[:, ig] = counts.loc[conv_singl[gene]]
            continue

        if gene in multiplets:
            countsg[:, ig] = counts.loc[conv_multi.loc[gene]].sum(axis=0)
            continue

        # missing, copy verbatim
        countsg[:, ig] = counts.loc[gene]

    print('Prepare for output')
    obs = pd.DataFrame([], index=counts.columns)
    obs['coverage'] = countsg.sum(axis=1)
    obs['ngenes'] = (countsg > 0).sum(axis=1)
    obs['Batch'] = obs.index.str.split('-').str.get(1)
    obs['CellBarcode'] = obs.index.str.split('-').str.get(0)
    var = pd.DataFrame([], index=genes)

    print('Normalize counts per 10000')
    countsg = (1e4 * countsg.T / obs['coverage'].values).T

    adata = anndata.AnnData(
        X=countsg,
        obs=obs,
        var=var,
        uns={
            'description': 'Endothelial differentiation in vitro',
            'normalisation': 'Counts per ten thousands',
            },
    )

    fn_out = f'{data_fdn}/normalised.h5ad'
    adata.write(fn_out, force_dense=True)
