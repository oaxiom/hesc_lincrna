import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from rpy2.robjects.packages import importr
#from gprofiler import gprofiler
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=300)

# Get back:
adata = sc.read('./filtered.h5ad')

size_factors = []

with open('size_factors.csv', 'rt') as oh:
    size_factors = [float(line.rstrip()) for line in oh if line and 'size_factors' not in line]

adata.obs['size_factors'] = size_factors

sc.pl.scatter(adata, 'size_factors', 'n_counts', show=False, save='size_factors_vs_ncounts.pdf')
sc.pl.scatter(adata, 'size_factors', 'n_genes', show=False, save='size_factors_vs_ncounts.pdf')
p = sb.distplot(size_factors, bins=50, kde=False)
p.get_figure().savefig('figures/size_factors.pdf')

#Normalize adata
adata.X /= adata.obs['size_factors'].values[:,None]
adata.X = np.clip(adata.X, 0, 1e10) # get rid of <0

sc.pp.log1p(adata)
adata.X = sp.sparse.csr_matrix(adata.X)

# Store the full data set in 'raw' as log-normalised data for statistical testing
adata.raw = adata

print('Combat')
sc.pp.combat(adata, key='replicate')

# resparsify:
print('Resparsify')
adata.X = np.clip(adata.X, 0, 1e10) # get rid of <0 to help sparsification:
adata.X = sp.sparse.csr_matrix(adata.X)

adata.write('./normed.h5ad')

