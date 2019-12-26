
'''

Get the CDS out of the GENCODE FASTA file;

'''

import glob, sys, os, gzip, numpy, math, re
import matplotlib.pyplot as plot
from glbase3 import glload, utils, expression, genelist, genome_sql, format
sys.path.append('../../')
import shared

fastas = genelist('gencode.v32.pc_transcripts.fa.gz', format=format.fasta, gzip=True)

newl = []
for fasta in fastas:
    enst = fasta['name'].split('|')
    cds = [i for i in enst if 'CDS:' in i][0].split(':')[1].split('-')
    newl.append({'enst': enst[0].split('.')[0], 'cds_local_locs': (int(cds[0]), int(cds[1])), 'tlength': len(fasta['seq'])})

newg = genelist()
newg.load_list(newl)
newg.save('gencode_cds.glb')
newg.saveTSV('gencode_cds.tsv')
