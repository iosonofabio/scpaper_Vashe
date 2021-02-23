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

import matplotlib.pyplot as plt


if __name__ == '__main__':

    print('Load data')
    data_fdn = '../../data/sequencing'
    fn_normalised = f'{data_fdn}/normalised.h5ad'
    adata = anndata.read_h5ad(fn_normalised)

    print('Plot coverage')
    fig, ax = plt.subplots(figsize=(3, 3))
    x = np.sort(adata.obs['coverage'])
    y = 1.0 - np.linspace(0, 1, len(x))
    ax.plot(x, y, lw=2, color='tomato')
    ax.set_xlabel('Coverage')
    ax.set_ylabel('Fraction of cells\nwith coverage > x')
    ax.set_xscale('log')
    fig.tight_layout()

    # Looks pretty good already, let's leave it for now

    plt.ion()
    plt.show()
