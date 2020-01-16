
'''

Get the CDS out of the GENCODE FASTA file;

'''

import glob, sys, os, gzip, numpy, math, re
import matplotlib.pyplot as plot
from glbase3 import glload, utils, expression, genelist, genome_sql, format
sys.path.append('../../')
import shared

fastas = genelist('gencode.v32.pc_transcripts.fa.gz', format=format.fasta, gzip=True)

truncated = 0
newl = []
for fasta in fastas:
    enst = fasta['name'].split('|')
    cds = [i for i in enst if 'CDS:' in i][0].split(':')[1].split('-')
    cds = (int(cds[0])-1, int(cds[1]))
    newl.append({'enst': enst[0], 'cds_local_locs': (cds[0], cds[1]), 'tlength': len(fasta['seq'])})

    if fasta['seq'][cds[0]:cds[0]+3] != 'ATG':
        print('! Warning: {0} does not have an ATG'.format(fasta['name']))
        truncated += 1
        continue
    elif len(fasta['seq'][cds[1]:cds[1]+3]) != 3:
        print('! Warning: {0} does not have an in frame STOP'.format(fasta['name']))
        truncated += 1
        continue
    elif (cds[1] - cds[0]) % 3 != 0:
        print('! Warning, {0} not divisible by 3'.format(fasta['name']))
        truncated += 1
        continue

newg = genelist()
newg.load_list(newl)
newg.save('gencode_cds.glb')
newg.saveTSV('gencode_cds.tsv')

print('Truncated transcripts/CDS frames:', truncated)
