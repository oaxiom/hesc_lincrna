"""

Pack the scRNA-seq data using scanpy, prep for scran normalisation

"""

import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData
from matplotlib import rcParams
from matplotlib import colors
plt.rcParams['figure.figsize'] = (8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 10
sc.settings.autoshow = False

adata = sc.read('raw_data.h5ad')

#mito_genes = adata.var_names.str.startswith('mt-')
# for each cell compute fraction of counts in mito genes vs. all genes
# the `.A1` is only necessary as X is sparse (to transform to a dense array after summing)
#adata.obs['percent_mito'] = np.sum(adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
# add the total counts per cell as observations-annotation to adata
#adata.obs['n_counts'] = adata.X.sum(axis=1).A1

sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-pre-norm-replicates.pdf')

# Base filtering for QC failures:
sc.pp.filter_cells(adata, min_genes=1000)
sc.pp.filter_cells(adata, max_genes=8000)
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=100000)
sc.pp.filter_genes(adata, min_cells=100) # Only filter genes here;
#adata = adata[adata.obs['percent_mito'] < 0.2, :]

# Filter out the mt genes to stop them etting into the most variable
#mito_genes = adata.var_names.str.startswith('mt-')
#mask = np.isin(adata.var_names, mito_genes, invert=True, assume_unique=True)
#adata = adata[:, mask] # also slices .var[]
# The above causes trouble for some reason in scran

sc.pl.violin(adata, ['n_genes','n_counts'], groupby='cell_type', size=0, log=False, cut=0, show=False, save='qc1.pdf')
sc.pl.violin(adata, ['n_genes','n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-replicates.pdf')

print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

adata.write('./filtered.h5ad')

oh = open('gene_names.filtered.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()

