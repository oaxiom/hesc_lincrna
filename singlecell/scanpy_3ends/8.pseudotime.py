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

sc.settings.figdir = 'pseudotime'

adata = sc.read('de.h5ad')

# UMAP is already in the neighbours, so don't regenerate;
#sc.pp.neighbors(adata, n_neighbors=20, use_rep='X', method='umap')
sc.tl.diffmap(adata)
sc.tl.dpt(adata, n_branchings=1, n_dcs=10)
sc.pl.diffmap(adata, color=['dpt_pseudotime', 'dpt_groups', 'de_clusters'], show=False, save='-diffmap-umap.pdf')
sc.pl.dpt_timeseries(adata, show=False, save='-diffmap-umap.pdf')
sc.pl.dpt_timeseries(adata, show=False, save='-timeseries-umap.pdf')

sc.pp.neighbors(adata, n_neighbors=20, use_rep='X', method='gauss')
sc.tl.diffmap(adata)
sc.tl.dpt(adata, n_branchings=1, n_dcs=10)
sc.pl.diffmap(adata, color=['dpt_pseudotime', 'dpt_groups', 'de_clusters'], show=False, save='-diffmap-gauss.pdf')
sc.pl.dpt_timeseries(adata, show=False, save='-timeseries-gauss.pdf')
