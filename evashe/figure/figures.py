# vim: fdm=indent
'''
author:     Fabio Zanini
date:       19/02/21
content:    Make figure panels for the paper.
'''
import os
import sys
import argparse
import numpy as np
import pandas as pd

import matplotlib as mpl
#mpl.rcParams['font.size'] = 7
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.append('/home/fabio/university/postdoc/singlet')
import singlet



if __name__ == '__main__':

    fig_fdn = '../../figures'

    pa = argparse.ArgumentParser()
    pa.add_argument('--save', action='store_true')
    args = pa.parse_args()

    print('Load data')
    data_fdn = '../../data/sequencing'
    fn_normalised = f'{data_fdn}/normalised.h5ad'
    ds = singlet.Dataset(
        dataset={
            'path': fn_normalised,
        })
    print('Loaded')

    print('Load umap')
    fn_umap = f'{data_fdn}/umap.tsv'
    vs = pd.read_csv(fn_umap, sep='\t', index_col=0)
    vs = vs.loc[ds.samplenames]
    ds.obs['umap1'] = vs['umap1']
    ds.obs['umap2'] = vs['umap2']
    ds.obs['leiden'] = vs['leiden'].astype(str)

    print('Sort by average Runx1 expression')
    dsa = ds.average('samples', by='leiden')
    ds.obs['cluster_new'] = ds.obs['leiden'].map({
        '3': '0',
        '2': '1',
        '0': '2',
        '6': '3',
        '4': '4',
        '1': '5',
        '5': '6',
    })

    print('Load pathway genes from Vashe')
    fn_pway = '../../data/gene_lists/Signalling pathways_VC.csv'
    genes_pway = pd.read_csv(fn_pway)
    ngenes = [59, 26, 24]
    genes_pway = {col.split(' ')[0]: genes_pway[col].iloc[1:ngenes[i]].values for i, col in enumerate(genes_pway.columns)}
    trans = {'Wisp1': 'Ccn4', 'Wisp2': 'Ccn5', 'Wisp3': 'Ccn6', 'Fgfr1op': 'Cep43'}
    for pway, genes_pw in genes_pway.items():
        genes_pway[pway] = [trans.get(g, g) for g in genes_pw]

    print('Dot plot')
    genes1 = [
        'Pecam1',
        'Cdh5',
        'Sox17',
        #'Notch1', 'Notch4', 'Bmpr2',
        #'Bmp2k',
        'Ptprc',
        'Runx1',
        'Gfi1', 'Gfi1b',
        'Tal1',
        #'Gata2', 'Erg', 'Fli1', 'Lyl1', 'Lmo2', 'Etv6',
        ]
    genes = list(genes1)
    #for pway, genes_pw in genes_pway.items():
    #    genes += list(genes_pw)
    #genes += list(ds.featurenames[ds.featurenames.str.startswith('Col')])
    #genes += list(ds.featurenames[ds.featurenames.str.startswith('Act')])

    cmap = plt.cm.get_cmap('Reds')
    size_fun = lambda x: 2 + (x * 130)
    fig, ax = plt.subplots(figsize=(4.5, 3.5))
    cluster_order = ['0', '1', '2', '3', '4', '5', '6']
    ds.plot.dot_plot(
            group_axis='samples',
            group_by='cluster_new',
            group_order=cluster_order,
            ax=ax,
            plot_list=genes[::-1],
            threshold=0.5,
            size_function=size_fun,
            color_log=True,
            layout='vertical',
            vmin=-1, vmax=1,
            cmap=cmap,
    )
    ax.set_ylabel('Gene')
    ax.set_xlabel('Cluster')
    for tk in ax.get_xticklabels():
        tk.set_rotation(0)

    percentage = [0, 25, 50, 75, 100]
    labels = [str(x) for x in percentage]
    sizes = [size_fun(0.01 * x) for x in percentage]
    hs = [ax.scatter([], [], s=s, color='grey') for s in sizes]
    ax.legend(hs, labels, loc='upper left', bbox_to_anchor=(1.05, 1), bbox_transform=ax.transAxes, title='Fraction expressing:')

    ax2 = ax.twinx()
    ax2.set_axis_off()
    expr = [0.0, 0.5, 1.0]
    labels = ['$0$', '$1$', '$10$']
    cmap = ax._singlet_dotmap['level_color_map']
    hs = [ax.scatter([], [], s=30, color=cmap(x)) for x in expr]
    ax2.legend(hs, labels, loc='upper left', bbox_to_anchor=(1.05, 0.35), bbox_transform=ax.transAxes, title='Exp level [cptt]:      ')

    fig.tight_layout()
    if args.save:
        fig.savefig(f'{fig_fdn}/dotplot.png', dpi=300)
        fig.savefig(f'{fig_fdn}/dotplot.svg')

    print('UMAP by cluster')
    clusters = ds.obs['cluster_new'].unique()
    cmap = cmap_clusters = dict(zip(clusters, sns.color_palette('husl', n_colors=len(clusters))))
    fig, ax = plt.subplots(figsize=(3, 3))
    ds.plot.scatter_reduced(
            ('umap1', 'umap2'),
            color_by='cluster_new',
            color_log=None,
            ax=ax,
            alpha=0.3,
            cmap=cmap,
    )
    ax.grid(False)
    ax.set_axis_off()
    for cl in clusters:
        ind = ds.obs['cluster_new'] == cl
        xm, ym = ds.obs.loc[ind, ['umap1', 'umap2']].mean(axis=0)
        ax.text(xm, ym, cl, ha='center', va='center')
    fig.tight_layout()
    if args.save:
        fig.savefig(f'{fig_fdn}/UMAP_cluster.png', dpi=300)

    print('UMAPs for Runx1 and Sox17')
    fig, axs = plt.subplots(2, len(clusters), figsize=(8, 2.5))
    cluster_order = ['0', '1', '2', '3', '4', '5', '6']
    for gene, axr in zip(['Runx1', 'Sox17'], axs):
        axr[0].set_xlabel(gene, rotation=0, ha='right')
        for cl, ax in zip(cluster_order, axr):
            if axr[0] == axs[0, 0]:
                ax.set_title(cl)
            ind = ds.obs['cluster_new'] == cl
            ds.obs['is_cluster'] = ind
            dsp = ds.split('is_cluster')
            dsp[True].plot.scatter_reduced(
                    ('umap1', 'umap2'),
                    color_by=gene,
                    color_log=True,
                    ax=ax,
                    s=35,
                    alpha=0.02,
                    cmap='viridis',
                    vmin=-1, vmax=2
            )
            dsp[False].plot.scatter_reduced(
                    ('umap1', 'umap2'),
                    ax=ax,
                    s=15,
                    alpha=0.01,
            )
            ax.set_axis_off()
            ax.grid(False)
            if gene == 'Runx1':
                ax.add_artist(plt.Rectangle(
                    (0.38, 1.05), 0.24, 0.23, facecolor='none',
                    transform=ax.transAxes,
                    edgecolor=cmap_clusters[cl], lw=2, clip_on=False),
                    )
    fig.text(0.01, 0.7, 'Runx1')
    fig.text(0.01, 0.25, 'Sox17')
    fig.tight_layout(rect=(0.05, 0, 1, 1))
    if args.save:
        fig.savefig(f'{fig_fdn}/UMAP_genes.png', dpi=300)


    print('Distributions')
    genes = ['Notch1', 'Notch4', 'Bmpr2']
    colors = ['seagreen', 'tomato', 'navy']
    cluster_order = ['0', '1', '2', '3', '4', '5', '6']
    fig, axs = plt.subplots(len(clusters), len(genes), figsize=(3.8, 6), sharex=True, sharey=True)
    for cl, axr in zip(cluster_order, axs):
        dsi = ds.query_samples_by_metadata('cluster_new == @cl', local_dict=locals())
        axr[0].set_ylabel(cl, rotation=0, ha='right')
        for gene, ax, color in zip(genes, axr, colors):
            if axr[0] == axs[0, 0]:
                ax.set_title(gene)
            x = np.sort(dsi.counts.loc[gene])
            y = 1.0 - np.linspace(0, 1, len(x))
            ax.plot(x, y, lw=2, color=color)
            ax.set_xscale('log')
            if ax == axs[-1, 1]:
                ax.set_xlabel('Gene expression [cptt]')
    fig.text(0.02, 0.5, 'Cluster', rotation=90)
    fig.text(0.07, 0.35, 'Fraction of cells expressing > x', rotation=90)
    fig.tight_layout(rect=(0.1, 0, 1, 1), h_pad=0, w_pad=0)
    if args.save:
        fig.savefig(f'{fig_fdn}/distributions.png', dpi=300)
        fig.savefig(f'{fig_fdn}/distributions.svg')
