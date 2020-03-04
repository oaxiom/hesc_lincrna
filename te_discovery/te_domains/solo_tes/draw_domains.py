
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

solo_tes = glload('../../solo_tes/solo_tes.glb')
print(solo_tes)

CDSs = glload('../../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
doms = solo_tes.map(genelist=CDSs, key='transcript_id')
print(doms)

tes_to_check = set(['HERVH', 'LTR7', 'FRAM', 'SVA', 'L2d', 'HERV-Fc2', 'HERVH48', 'AluJb', 'MLT1C2'])

#genes = set(solo_tes['transcript_id'])

for n, gene in enumerate(doms):
    #print(gene['name'].split(' ')[0])
    #print(gene)
    #if gene['transcript_id'].split(' ')[0] not in genes:
    #    continue

    all_TEs = set([g['dom'] for g in gene['doms']])
    tes = [te for te in all_TEs if te in tes_to_check]

    if not tes:
        continue

    path = '%s/%s/' % (draw, ';'.join(tes))
    if not os.access('%s' % (path), os.R_OK | os.W_OK):
        os.mkdir('%s' % (path))

    draw_domains_share.draw_domain(gene, '%s/%s.%s.%s.%s' % (path, gene['name'], gene['transcript_id'], gene['enst'], draw), dfam)


