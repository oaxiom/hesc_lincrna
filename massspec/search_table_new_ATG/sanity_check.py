
# Do a word match against the peeptides from the GENCODE FASTA

import sys, os
from glbase3 import *
sys.path.append('../../')
import shared

gencode_peptide_fastas = genelist('../../transcript_assembly/get_CDS/gencode.v32.pc_translations.fa.gz', format=format.fasta, gzip=True)

print(gencode_peptide_fastas)

gencode_peptide_fastas_lookup = {}
for gene in gencode_peptide_fastas:
    name = gene['name'].split('|')[6]
    if name not in gencode_peptide_fastas_lookup:
        gencode_peptide_fastas_lookup[name] = []
    gencode_peptide_fastas_lookup[name].append(gene['seq'])

dna_fasta = glload('../../transcript_assembly/fasta/transcripts.glb')
genes = glload('../../te_discovery/CDS_insertions/table_new_ATG.glb')
with_seq = genes.map(genelist=dna_fasta, key='transcript_id')

res = {'match': 0,
    'no-match': 0}

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

    aa = ''.join([a for a in aa if a != '_'])

    name = gene['name'].split(' ')[0]
    found = False
    if name in gencode_peptide_fastas_lookup:
        for seq in gencode_peptide_fastas_lookup[gene['name'].split(' ')[0]]:

            if seq == aa:
                print()
                print(gene)
                print(cdsl, cdsr, cdsr-cdsl)
                print(name)
                print(seq)
                print(aa)
                print()
                found = True
                break
    else:
        print('{0} not found in gencode_peptide_fastas_lookup'.format(name))

    if found:
        res['match'] += 1
    else:
        res['no-match'] += 1


print(res)

