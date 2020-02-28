import os, sys, glob, numpy
import scanpy as sc
from glbase3 import genelist, glload
sc.settings.verbosity = 1
sc.set_figure_params(dpi=200, dpi_save=200)
sc.settings.figdir = 'markers'

dfam = glload('../../../../te_discovery/dfam/dfam_annotation.glb') #, format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
dfam = {d['name']: '{0}:{1}'.format(d['type'], d['subtype'], d['name']) for d in dfam}

adata = sc.read('../../learned.h5ad')

for filename in glob.glob('../te_containing/*.glb'):
    print(filename)
    # Just do a simple annotation of yh the TE-contianig DE genes from the scRNA-seq
    stub = os.path.split(filename)[1].split('.')[0].replace('de_genes-', '')
    de_list = glload(filename)

    sc.settings.figdir = stub

    if stub == 'grp5':
        continue

    for item in de_list:
        unq_doms = set([i['dom'] for i in item['doms']]) # Only draw once per TE family;
        for te in unq_doms:
            type_subtype = dfam[te]

            sc.settings.figdir = '{0}/{1}'.format(stub, type_subtype)



            sc.pl.umap(adata,
                color=item['transcript_id'],
                size=10,
                legend_loc='on data',
                vmax=3,
                show=False,
                save='-{0}{1}-markers-{2}-{3}.pdf'.format(type_subtype, te, item['transcript_id'], item['name'])
                )
