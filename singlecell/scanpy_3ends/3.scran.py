
import logging, os, sys
import rpy2
from rpy2 import robjects
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import anndata2ri
import anndata
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sb
#import rpy2.rinterface_lib.callbacks
import anndata2ri
from matplotlib import rcParams
from matplotlib import colors
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=300)
sc.settings.figdir = 'markers'
sc.settings.autoshow = False

# Ignore R warning messages
#Note: this can be commented out to get more verbose R output
#rpy2.rinterface_lib.callbacks.logger.setLevel(logging.ERROR)
# Automatically convert rpy2 outputs to pandas dataframes
pandas2ri.activate()
anndata2ri.activate()
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
#sc.set_figure_params(dpi=200, dpi_save=300)
sc.logging.print_versions()

adata = sc.read('filtered.h5ad')
print(adata)
print('Number of cells: {:d}'.format(adata.n_obs))

adata_pp = adata.copy()
sc.pp.normalize_per_cell(adata_pp, counts_per_cell_after=1e6)
sc.pp.log1p(adata_pp)
sc.pp.pca(adata_pp, n_comps=20, svd_solver='arpack')
sc.pp.neighbors(adata_pp)
sc.tl.leiden(adata_pp, key_added='groups', resolution=0.3)

input_groups = adata_pp.obs['groups']
data_mat = adata.X.T
data_mat = sp.sparse.csc_matrix(data_mat)

# -i data_mat -i input_groups -o size_factors
scran = importr('scran')
size_factors = scran.computeSumFactors(data_mat, clusters=input_groups, **{'min.mean': 0.1})
print(size_factors)

del adata_pp

adata.obs['size_factors'] = size_factors

#sc.pl.scatter(adata, 'size_factors', 'n_counts', save='-plot1.pdf')
#sc.pl.scatter(adata, 'size_factors', 'n_genes', save='-plot2.pdf')

#sb.distplot(size_factors, bins=50, kde=False)

adata.layers["counts"] = adata.X.copy()
adata.X /= adata.obs['size_factors'].values[:,None]
sc.pp.log1p(adata)
adata.X = np.clip(adata.X, 0, 100000000) # get rid of <0
adata.X = sp.sparse.csr_matrix(adata.X)
adata.raw = adata # You only need to do this if you do batch correction

# combat batch correction could go here:
sc.pp.combat(adata, key='replicate')

# resparsify:
adata.X = np.clip(adata.X, 0, 1e9) # get rid of <0
adata.X = sp.sparse.csr_matrix(adata.X)

adata.write('./normed.h5ad')

