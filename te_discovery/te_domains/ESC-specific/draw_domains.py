
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

[os.remove(f) for f in glob.glob('%s/*.%s' % (draw, draw))]

doms = glload('../transcript_table_HSC_SR_PB_merged.mapped.glb')
genes_to_do = set(genelist('../../../fantom5/som_results/H9 embryonic stem cells.tsv', format={'enst': 0})['enst']) # This is from the SOM, done on FANTOM5 data. See DPre
#gencode = glload(os.path.expanduser('~/hg38/hg38_gencode_v29.glb'))
gencode_db = genome_sql(filename=os.path.expanduser('~/hg38/hg38_gencode_v29.sql'))

for n, gene in enumerate(doms):
    # I have to pre-load the exonStarts:
    if gene['enst'] not in genes_to_do:
        continue

    draw_domains_share.draw_domain(gene,
        '%s/%s.%s.%s' % (draw, gene['name'], gene['enst'], draw),
        gencode_db)

    if (n+1) % 1000 == 0:
        print('Processed: {:,} domains'.format(n+1))
        #break
