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

transcript_id = glload('../../transcript_assembly/packed/all_genes.glb')
transcript_id = {i['transcript_id']: i for i in transcript_id}

de_leiden = 'de_clusters'

adata = sc.read('./de.h5ad')

sc.pl.rank_genes_groups(adata, n_genes=25, sharey=True, show=False, save='genes-top25.pdf')
sc.pl.rank_genes_groups(adata, key='rank_genes_groups', show=False, save='genes.pdf')
sc.pl.rank_genes_groups_dotplot(adata, key='rank_genes_groups', show=False, save='genes-top25.pdf')

topall = pd.DataFrame(adata.uns['rank_genes_groups']['names']) # get all;
fcs = pd.DataFrame(adata.uns['rank_genes_groups']['logfoldchanges'])
padj = pd.DataFrame(adata.uns['rank_genes_groups']['pvals_adj'])

# Matrix of DE genes:

groups = list(topall.columns.values)

newcols = {}

for group in groups:
    newcols[group] = []

    t = zip([i[group] for i in adata.uns['rank_genes_groups']['names']], [i[group] for i in adata.uns['rank_genes_groups']['logfoldchanges']], [i[group] for i in adata.uns['rank_genes_groups']['pvals_adj']])

    print('Group: {0}'.format(group))

    for item in t:
        if item[1] < 2: # fold change
            continue
        if item[2] > 0.01: # just in case
            continue

        newcols[group].append({'transcript_id': item[0], 'log2FC': item[1], 'q': item[2],
            'ensg': transcript_id[item[0]]['ensg'], 'enst': transcript_id[item[0]]['enst'], 'name': transcript_id[item[0]]['name']})

# join all and draw a dotplot:
for group in newcols:
    print('Top 10:\n', newcols[group][0:10])
    if newcols[group]:
        gl = genelist()
        gl.load_list(newcols[group])
        gl.saveTSV('gls/de_genes-grp{0}.tsv'.format(group))
        gl.save('gls/de_genes-grp{0}.glb'.format(group))

        genes = [i['transcript_id'] for i in newcols[group]]
        #sc.pl.dotplot(adata, genes, groupby=de_leiden, dot_max=0.7, dendrogram=True, standard_scale='var', show=False, save='de-grp{0}.pdf'.format(group))
        #sc.pl.matrixplot(adata, genes, groupby=de_leiden, dendrogram=True, standard_scale='var', show=False, save='de-grp{0}.pdf'.format(group))

for grp in newcols:
    if not newcols[grp]:
        continue
    for k in newcols[grp]:
        title = k['name']
        sc.pl.umap(adata, color=k['transcript_id'], size=30, legend_loc='on data',
            title=title,
            vmin=0, vmax=3,
            show=False, save='-markers-grp{0}-{1}-{2}.pdf'.format(grp, k['transcript_id'], k['name']))
        #sc.pl.violin(adata, [k], groupby='disease', size=0, log=False, cut=0, show=False, save='markers-{0}-disease.pdf'.format(k))
        #sc.pl.violin(adata, [k], groupby='cell_type', size=0, log=False, cut=0, show=False, save='markers-{0}-cell_type.pdf'.format(k))
