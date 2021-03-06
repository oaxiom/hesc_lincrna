
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

[os.remove(f) for f in glob.glob('%s/*/*/*.%s' % (draw, draw))]
[os.remove(f) for f in glob.glob('%s/*/*.%s' % (draw, draw))]

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
CDSs = glload('../../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
CDSs = {i['transcript_id']: i for i in CDSs}
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

for n, gene in enumerate(doms):
    destination, alpha = shared.classify_transcript(gene['name'])
    if not os.access('%s/%s' % (draw, destination), os.R_OK | os.W_OK):
        os.mkdir('%s/%s' % (draw, destination))

    path = '%s/%s/%s' % (draw, destination, alpha)

    if not os.access(path, os.R_OK | os.W_OK):
        os.mkdir(path)

    if gene['transcript_id'] in CDSs:
        gene['cds_local_locs'] = CDSs[gene['transcript_id']]['cds_local_locs']

    draw_domains_share.draw_domain(gene, '%s/%s.%s.%s.%s' % (path, gene['name'], gene['transcript_id'], gene['enst'], draw), dfam)

    if (n+1) % 1000 == 0:
        print('Processed: {:,} domains'.format(n+1))
        #break
