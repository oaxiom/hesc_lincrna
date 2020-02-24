import logging, matplotlib, os, sys, glob
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import pandas as pd
from glbase3 import genelist
plt.rcParams['figure.figsize']=(8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype']=42
matplotlib.rcParams['font.size']=10

from glbase3 import genelist, glload

sc.settings.figdir = 'diffexp'

[os.remove(f) for f in glob.glob('{0}/*.pdf'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.glb'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.tsv'.format(sc.settings.figdir))]

de_leiden = 'leiden_r0.70'

# Merge some clusters that have hardly any difference:


transcipt_ids = glload('../../transcript_assembly/packed/all_genes.glb')
transcipt_ids = {i['transcript_id']: i for i in transcipt_ids}

adata = sc.read('./learned.h5ad')

print(adata)

old_to_new = {
    '0': '0',
    '1': '0',
    '2': '1',
    '3': '0',
    '4': '0',
    '5': '2',
    '6': '0',
    '7': '3',
    '8': '4',
    '9': '5',
    }

adata.obs['de_clusters'] = (
    adata.obs['leiden_r0.70']
    .map(old_to_new)
    .astype('category')
    )

# I merge some unimportant clusters:
sc.pl.umap(adata, color=[de_leiden, 'de_clusters'], size=10, legend_loc='on data', show=False, save='-new-umap.pdf')

sc.tl.rank_genes_groups(adata, 'de_clusters', method='t-test_overestim_var', n_genes=10000)
#sc.tl.filter_rank_genes_groups(adata, min_fold_change=1)
adata.write('./de.h5ad')
