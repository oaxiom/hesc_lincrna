
import sys, os, glob
from glbase3 import *
from collections import defaultdict
sys.path.append('../../')
import shared

dna_fasta = glload('../../transcript_assembly/fasta/transcripts.glb')

gencode_peptide_fastas = genelist('../../transcript_assembly/get_CDS/gencode.v32.pc_translations.fa.gz', format=format.fasta, gzip=True)
print(gencode_peptide_fastas)
gencode_peptide_fastas_lookup = {}
for gene in gencode_peptide_fastas:
    name = gene['name'].split('|')[6]
    if name not in gencode_peptide_fastas_lookup:
        gencode_peptide_fastas_lookup[name] = []
    gencode_peptide_fastas_lookup[name].append(gene['seq'])


def d():
    return {'match': 0, 'no-match': 0}
res = defaultdict(d)

for filename in glob.glob('../../te_discovery/CDS_insertions/*.glb'):
    if 'disrupt_coding' in filename or 'no_coding' in filename:
        continue # skip, as no CDS!
    stub = os.path.split(filename)[1].replace('.glb', '')

    genes = glload(filename)
    with_seq = genes.map(genelist=dna_fasta, key='transcript_id')

    oh = open('fasta/{0}.fasta'.format(stub), 'w')

    for gene in with_seq:
        cdsl = gene['cds_local_locs'][0]
        cdsr = gene['cds_local_locs'][1]
        CDS = gene['seq'][cdsl:cdsr]

        aa = shared.translateAA(CDS)
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
                    found = True
                    break
        else:
            print('{0} not found in gencode_peptide_fastas_lookup'.format(name))

        if found:
            res[stub]['match'] += 1 # don't add to the FASTA
            continue
        else:
            res[stub]['no-match'] += 1

        oh.write('>{0}|{1}|{2}\n'.format(gene['name'].replace(' ', ''), gene['transcript_id'], gene['enst']))
        oh.write(''.join([a for a in aa if a != '_']))
        oh.write('\n')

    oh.close()

print()
for filename in glob.glob('../../te_discovery/CDS_insertions/*.glb'):
    if 'disrupt_coding' in filename or 'no_coding' in filename:
        continue # skip, as no CDS!
    stub = os.path.split(filename)[1].replace('.glb', '')

    print(stub)
    for r in ['match', 'no-match']:
        print('  {0}: {1}'.format(r, res[stub][r]))




