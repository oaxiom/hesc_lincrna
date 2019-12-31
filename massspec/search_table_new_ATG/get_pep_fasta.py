
import sys, os
from glbase3 import *
sys.path.append('../../')
import shared

dna_fasta = glload('../../transcript_assembly/fasta/transcripts.glb')
genes = glload('../../te_discovery/CDS_insertions/table_new_ATG.glb')
with_seq = genes.map(genelist=dna_fasta, key='transcript_id')

oh = open('result.fasta', 'w')

for gene in with_seq:
    cdsl = gene['cds_local_locs'][0]
    cdsr = gene['cds_local_locs'][1]
    CDS = gene['seq'][cdsl:cdsr]

    aa = shared.translateAA(shared.split3(CDS))
    if '_' not in aa:
        print('Warning {0} {1} has no STOP'.format(gene['transcript_id'], gene['name']))
    if aa.count('_') > 1:
        print('Warning {0} {1} has more than 1 STOP!'.format(gene['transcript_id'], gene['name']))
        1/0

    oh.write('>{0}|{1}|{2}\n'.format(gene['name'], gene['transcript_id'], gene['enst']))
    oh.write(''.join([a for a in aa if a != '_']))
    oh.write('\n')
oh.close()





