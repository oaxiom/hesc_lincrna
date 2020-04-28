import logging, matplotlib, os, sys
import anndata
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=300)

adata = sc.read('normed.h5ad')
print(adata)

print('Number of cells: {:d}'.format(adata.n_obs))

sc.pp.highly_variable_genes(adata, flavor='cell_ranger', n_top_genes=4000)
sc.pl.highly_variable_genes(adata, show=False, save='highly_variable.pdf')

# Calculate the visualizations
sc.pp.pca(adata, n_comps=20, use_highly_variable=True, svd_solver='arpack') # PC=20 from Nature paper
sc.pp.neighbors(adata)
sc.tl.tsne(adata, n_jobs=3) #Note n_jobs works for MulticoreTSNE, but not regular implementation
sc.tl.umap(adata, min_dist=2.0)
#sc.tl.diffmap(adata)

sc.pl.pca_variance_ratio(adata, log=True, show=False, save='pca_variance.pdf')

# Perform clustering - using highly variable genes
res = [1.0, 0.6, 0.5, 0.4, 0.35, 0.3,
    0.25, 0.2, 0.15, 0.1
    ]
for r in res:
    sc.tl.leiden(adata, resolution=r, key_added='leiden_r{0:.2f}'.format(r))

adata.write('./learned.h5ad')

todraw = ['cell_type', 'replicate'] + ['leiden_r{0:.2f}'.format(r) for r in res]

#Visualize the clustering and how this is reflected by different technical covariates
sc.pl.tsne(adata, color=todraw, size=10, legend_loc='on data', show=False, save='tsne.pdf')
sc.pl.umap(adata, color=todraw, size=10, legend_loc='on data', show=False, save='umap.pdf')

