
import sys, os, glob
from glbase3 import *
from collections import defaultdict
sys.path.append('../../')
import shared

dna_fasta = glload('../../transcript_assembly/fasta/transcripts.glb')
pcat14 = dna_fasta.get(key='transcript_id', value='HPSCSR.312557.3')[0]['seq']

# get the longest ORFS:

# NEed a custom translateAA:
def translateAA_skip_incomplete_codons(seq):
    split = shared.split3(seq)
    split = [s for s in split if len(s) == 3]
    return [shared.table[codon.upper()] for codon in split]


frame1 = translateAA_skip_incomplete_codons(pcat14)
frame2 = translateAA_skip_incomplete_codons(pcat14[1:])
frame3 = translateAA_skip_incomplete_codons(pcat14[2:])

# rank all possible ORFs:
ORFs = []
for seq in [frame1, frame2, frame3]:
    corf = None
    for aa in seq:
        if corf is None and aa == 'M':
            corf = []
        elif corf and aa == '_':
            ORFs.append((len(corf), ''.join(corf)))
            corf = None

        if corf is not None:
            corf.append(aa)

# Not so many, so output all >50 aa
oh = open('all_orfs_gt50.fa', 'wt')
for idx, orf in enumerate(ORFs):
    print(ORFs)
    if orf[0] > 50:
        oh.write('>orf-{0}\n'.format(idx))
        oh.write('{0}\n\n'.format(orf[1]))
oh.close()
