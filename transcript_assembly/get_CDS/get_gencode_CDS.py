
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
for f in fastas:
    f = f['name'].split('|')
    cds = [i for i in f if 'CDS:' in i][0].split(':')[1].split('-')
    newl.append({'enst': f[0], 'cds_local_locs': (cds[0], cds[1])})

newg = genelist()
newg.load_list(newl)
newg.save('gencode_cds.glb')
newg.saveTSV('gencode_cds.glb')
