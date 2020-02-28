
'''

Draw protein domain-style plots, but containing the locations of the TEs

The plots should be the mRNA, with a protein coding exon (if present), and the location of the TE;

This script collectes the data

'''

import glob, sys, os, gzip
from glbase3 import glload, utils, expression, genelist, genome_sql

sys.path.append('../')
import draw_domains_share
sys.path.append('../../../')
import shared

draw = 'pdf'
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
gencode_db = genome_sql(filename=os.path.expanduser('~/hg38/hg38_gencode_v29.sql'))
#all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')

for filename in glob.glob('../../CDS_insertions/*.glb'):
    typ = os.path.split(filename)[1].split('.')[0]

    if os.access('{0}'.format(typ), os.R_OK | os.W_OK):
        [os.remove(f) for f in glob.glob('{0}/*.pdf'.format(typ))]
    else:
        os.mkdir('{0}'.format(typ))

    tes = glload(filename)
    #tes = all_genes.map(genelist=tes, key='transcript_id')

    for n, gene in enumerate(tes):
        draw_domains_share.draw_domain(gene, '{0}/{1}_{2}_{3}.{4}'.format(typ, gene['name'], gene['transcript_id'], gene['enst'], draw), dfam)

