
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

#[os.remove(f) for f in glob.glob('%s/*.%s' % (draw, draw))]

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
CDSs = glload('../../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

doms = doms.map(genelist=CDSs, key='transcript_id')
print(doms)

for n, gene in enumerate(doms):
    found = False
    for d in gene['doms']:
        if gene['coding'] == 'noncoding':
            continue

        if 'HERVE' in d['dom']:
            found = True
            break

    if found:
        draw_domains_share.draw_domain(gene, '%s/%s.%s.%s.%s' % (draw, gene['name'], gene['transcript_id'], gene['enst'], draw), dfam)

    if (n+1) % 1000 == 0:
        print('Processed: {:,} domains'.format(n+1))
        #break
