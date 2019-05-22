
'''

Draw protein domain-style plots, but containing the locations of the TEs

The plots should be the mRNA, with a protein coding exon (if present), and the location of the TE;

This script collectes the data

Get a list of all transcripts that contain a TE, add their genome coordinates from the GTF

'''

import glob, sys, os, gzip
from glbase3 import glload, utils, expression, genelist

all_genes = glload('../../transcript_assembly/packed/all_genes.glb')
print(all_genes)
sam = '../hmmer_dfam/data/transcript_table_HSC_SR_PB_merged.transcripts.glb'
print('... %s' % sam)
doms = glload(sam)

doms = doms.map(all_genes, key='transcript_id')
print(doms)
doms.save('transcript_table_HSC_SR_PB_merged.mapped.glb') # fixed name version
doms.saveTSV('transcript_table_HSC_SR_PB_merged.mapped.tsv') # fixed name version
print(doms)
doms[0:100].save('transcript_table_HSC_SR_PB_merged.mapped.small.glb') # fixed name version
