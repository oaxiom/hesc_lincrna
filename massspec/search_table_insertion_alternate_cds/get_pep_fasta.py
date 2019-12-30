

from glbase3 import *

table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
    }

def split3(s):
    return [s[i:i+3] for i in range(0, len(s), 3)]

def translateAA(seq):
    return [table[codon.upper()] for codon in seq]

dna_fasta = glload('../../transcript_assembly/fasta/transcripts.glb')
genes = glload('../../te_discovery/CDS_insertions/table_insertion_alternate_cds.glb')
with_seq = genes.map(genelist=dna_fasta, key='transcript_id')

oh = open('result.fasta', 'w')

for gene in with_seq:
    cdsl = gene['cds_local_locs'][0]
    cdsr = gene['cds_local_locs'][1]
    CDS = gene['seq'][cdsl:cdsr]

    aa = translateAA(split3(CDS))
    if '_' not in aa:
        print('Warning {0} {1} has no STOP'.format(gene['transcript_id'], gene['name']))
    if aa.count('_') > 1:
        print('Warning {0} {1} has more than 1 STOP!'.format(gene['transcript_id'], gene['name']))
        1/0

    oh.write('>{0}|{1}|{2}\n'.format(gene['name'], gene['transcript_id'], gene['enst']))
    oh.write(''.join([a for a in aa if a != '_']))
    oh.write('\n')
oh.close()





