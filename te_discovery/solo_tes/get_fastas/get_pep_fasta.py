
import sys, os, glob
from glbase3 import *
from collections import defaultdict
sys.path.append('../../../')
import shared

dna_fasta = glload('../../../transcript_assembly/fasta/transcripts.glb')
CDSs = glload('../../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
CDSs = {i['transcript_id']: i for i in CDSs}

def d():
    return {'match': 0, 'no-match': 0}
res = defaultdict(d)

solotes = glload('../solo_tes.glb')

with_seq = solotes.map(genelist=dna_fasta, key='transcript_id')
stub = 'solo_tes'
oh = open('{0}.fasta'.format(stub), 'w')

for gene in with_seq:
    if gene['coding'] != 'coding':
        continue

    gene['cds_local_locs'] = CDSs[gene['transcript_id']]['cds_local_locs']

    cdsl = gene['cds_local_locs'][0]
    cdsr = gene['cds_local_locs'][1]
    CDS = gene['seq'][cdsl:cdsr]

    #print(shared.split3(CDS))
    try:
        aa = shared.translateAA(CDS)
    except KeyError:
        print('Warning {0} {1} malformed'.format(gene['transcript_id'], gene['name']))
        continue

    if '_' not in aa:
        print('Warning {0} {1} has no STOP'.format(gene['transcript_id'], gene['name']))
    if aa.count('_') > 1:
        print('Warning {0} {1} has more than 1 STOP!'.format(gene['transcript_id'], gene['name']))
        1/0

    aa = ''.join([a for a in aa if a != '_'])



    '''
    name = gene['name'].split(' ')[0]
    found = False
    if name in gencode_peptide_fastas_lookup:
        for seq in gencode_peptide_fastas_lookup[gene['name'].split(' ')[0]]:

            if seq == aa:
                found = True
                break
    else:
        print('{0} not found in gencode_peptide_fastas_lookup'.format(name))

    if found:
        res[stub]['match'] += 1 # don't add to the FASTA
        continue
    else:
        res[stub]['no-match'] += 1
    '''

    oh.write('>{0}|{1}|{2}\n'.format(gene['name'].replace(' ', ''), gene['transcript_id'], gene['enst']))
    oh.write(''.join([a for a in aa if a != '_']))
    oh.write('\n')

oh.close()

print()

print(stub)
for r in ['match', 'no-match']:
    print('  {0}: {1}'.format(r, res[stub][r]))




