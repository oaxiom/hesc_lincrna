import os, sys, glob, numpy
import scanpy as sc
from glbase3 import genelist, glload
sc.settings.verbosity = 1
sc.set_figure_params(dpi=200, dpi_save=200)
sc.settings.figdir = 'markers'

adata = sc.read('../../learned.h5ad')

for filename in glob.glob('../te_containing/*.glb'):
    print(filename)
    # Just do a simple annotation of yh the TE-contianig DE genes from the scRNA-seq
    stub = os.path.split(filename)[1].split('.')[0].replace('de_genes-', '')
    de_list = glload(filename)

    sc.settings.figdir = stub

    for item in de_list:
        for te in item['doms']:
            sc.settings.figdir = '{0}/{1}'.format(stub, te['dom'])

            sc.pl.umap(adata,
                color=item['transcript_id'],
                size=10,
                legend_loc='on data',
                vmax=3,
                show=False,
                save='-{0}-markers-{1}-{2}.pdf'.format(te['dom'],  item['transcript_id'], item['name'])
                )
