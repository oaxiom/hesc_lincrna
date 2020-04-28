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

sc.settings.figdir = 'diffexp-summaries'

[os.remove(f) for f in glob.glob('{0}/*.pdf'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.glb'.format(sc.settings.figdir))]
[os.remove(f) for f in glob.glob('gls/*.tsv'.format(sc.settings.figdir))]

de_leiden = 'leiden_r1.00'

# Merge some clusters that have hardly any difference:


transcipt_ids = glload('../../transcript_assembly/packed/all_genes.glb')
transcipt_ids = {i['transcript_id']: i for i in transcipt_ids}

adata = sc.read('./learned.h5ad')

de_clusters = de_leiden


# Merge clusters that have no DE genes:

print(adata)

old_to_new = {
    '0':  '0',
    '2':  '0',
    '13': '0',
    '3':  '0',
    '4':  '0',
    '7':  '0',
    '8':  '0',
    '9':  '0',
    '6':  '0',

    '1':  '1',
    '10': '1',

    '5':  '2',

    '11': '3', # keep
    '12': '4', # keep
    }

adata.obs['de_clusters'] = (
    adata.obs['leiden_r1.00']
    .map(old_to_new)
    .astype('category')
    )

adata.obs.to_csv('cell_data.csv')

de_clusters = 'de_clusters'


sc.pl.umap(adata, color=[de_leiden, de_clusters], size=10, legend_loc='on data', show=False, save='-new-umap.pdf')
sc.tl.rank_genes_groups(adata, de_clusters, method='t-test_overestim_var', n_genes=10000)
adata.write('./de.h5ad')
