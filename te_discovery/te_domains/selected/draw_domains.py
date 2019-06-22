
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

draw = 'png'

#[os.remove(f) for f in glob.glob('%s/*.%s' % (draw, draw))]

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
gencode_db = genome_sql(filename=os.path.expanduser('~/hg38/hg38_gencode_v29.sql'))
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

genes = set(['SOX2', 'NANOG', 'SALL4', 'LIN28A', 'LIN28B', 'SALL1', 'POU5F1', 'BRCA1', 'BRCA2',
    'DPPA2', 'DPPA3', 'DPPA5', 'PRDM14', 'JARID2', 'SALL2', 'SALL3', 'TCF3',
    'DNMT3L', 'LEFTY2', 'FGF4', 'NODAL', 'CER1', 'NLRP7', 'SFRP1', 'ZIC2', 'KDR',
    'ZFP42', 'C9ORF135', 'ST6GAL1', 'LRP4', 'MSTO1', 'PRODH',# From Pontis et al., 2019 CSC
    'SFRP2', 'OTX2', 'KLF4', 'KLF5',
    'HNRNPK', 'HNRNPU'])

for n, gene in enumerate(doms):
    #print(gene['name'].split(' ')[0])
    if gene['name'].split(' ')[0] not in genes:
        continue

    draw_domains_share.draw_domain(gene, '%s/%s.%s.%s.%s' % (draw, gene['name'], gene['transcript_id'], gene['enst'], draw), gencode_db, dfam)

    if (n+1) % 1000 == 0:
        print('Processed: {:,} domains'.format(n+1))
        #break
