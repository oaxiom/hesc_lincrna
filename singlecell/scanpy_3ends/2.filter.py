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

sc.pl.violin(adata, ['n_genes', 'n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-pre-norm-replicates.pdf')

# Base filtering for QC failures:
sc.pp.filter_cells(adata, min_genes=1500)
sc.pp.filter_cells(adata, max_genes=8000)
sc.pp.filter_cells(adata, min_counts=3000)
sc.pp.filter_cells(adata, max_counts=50000)
sc.pp.filter_genes(adata, min_cells=100) # Only filter genes here;

sc.pl.violin(adata, ['n_genes','n_counts'], groupby='cell_type', size=0, log=False, cut=0, show=False, save='qc1.pdf')
sc.pl.violin(adata, ['n_genes','n_counts'], groupby='replicate', size=0, log=False, cut=0, show=False, save='qc1-replicates.pdf')

print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

adata.write('./filtered.h5ad')

oh = open('gene_names.filtered.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()


# rows = genes; cols=cells;
print('Save dense matrix')
data_mat = adata.X.T.tocsr()
# Save out a dense array of the sparse array:
oh = open('dense_array.tsv', 'w')
for i in range(data_mat.shape[0]):

    if i % 1000 == 0:
        print('{0}/{1}'.format(i, data_mat.shape[0]))

    oh.write('\t'.join([str(i) for i in data_mat.getrow(i).toarray()[0,:]]))
    oh.write('\n')
oh.close()
