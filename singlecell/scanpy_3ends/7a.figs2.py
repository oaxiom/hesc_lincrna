import logging, matplotlib, os, sys
import anndata
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=300)

adata = sc.read('learned.h5ad')

godsnot_64_noyell = [ # default_64
    # "#000000",  # remove the black, as often, we have black colored annotation
    #"#FFFF00", # 0 CD14+ mono
    #"#FF34FF", # "#1CE6FF", # 1 CD4 T
    '#FF4A46', # "#FF34FF", # 2 CD8 T
    #'#008941', # "#FF4A46", # 3 B cells
    '#1CE6FF', # "#008941", # 4 NK cells
    "#006FA6", # 5 FCGR+
    "#A30059", # 6 ?
    "#FFDBE5", # 7 DC
    "#7A4900", # 8 ?
    "#0000A6", #
    "#63FFAC", #
    "#B79762",]

#Visualize the clustering and how this is reflected by different technical covariates
todo = ['cell_type', 'replicate']
sc.pl.umap(adata, color=todo, size=20, legend_loc='right margin', show=False, save='umap-2.pdf')

todo = ['n_counts', 'n_genes']
sc.pl.umap(adata, color=todo, size=20, color_map='summer', legend_loc='right margin', show=False, save='umap-3.pdf')
