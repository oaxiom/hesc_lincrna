import os, sys, glob, numpy
import matplotlib.pyplot as plot
from matplotlib import gridspec
from glbase3 import genelist, glload

TEs = glload('../../../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')

dfam = glload('../../../../te_discovery/dfam/dfam_annotation.glb') #, format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
dfam = {d['name']: '{0}:{1}:{2}'.format(d['type'], d['subtype'], d['name']) for d in dfam}

for filename in glob.glob('../../gls/*.glb'):
    print(filename)
    # Just do a simple annotation of yh the TE-contianig DE genes from the scRNA-seq
    stub = os.path.split(filename)[1].split('.')[0].replace('de_genes-', '')
    de_list = glload(filename)

    annot = de_list.map(genelist=TEs, key='transcript_id')

    for item in annot:
        item['TES'] = '; '.join(list(set([d['dom'] for d in item['doms']])))

    annot._optimiseData()

    annot = annot.getColumns(['ensg', 'transcript_id', 'name', 'gene_symbol', 'TES'])

    annot.save('{0}.glb'.format(stub))
    annot.saveTSV('{0}.tsv'.format(stub))
