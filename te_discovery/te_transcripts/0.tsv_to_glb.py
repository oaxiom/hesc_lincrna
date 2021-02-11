
'''

Build the final annotation tables;

'''

import glob, sys, os, gzip
from glbase3 import glload, utils, expression, genelist


f = {
    'transcript_id': 0, 	
    'doms': {'code': 'eval(column[1])'},
    'ensg': 2,
    'enst': 3,
    'name': 4,
    'gene_symbol': 5,
    'loc': 'location(loc=column[6])',
    'exonCounts': 7,
    'exonStarts': 8,
    'exonEnds': 9,
    'strand': 10,
    'tags': 11,
    'transcript_class': 12,
    'coding': 13,
    'expression': 14,
    'TPM': 15,
    }

doms = genelist('transcript_table_merged.mapped.tsv.gz', force_tsv=True, gzip=True,
    format=f)

print(doms)

doms.save('transcript_table_merged.mapped.glb') # fixed name version

