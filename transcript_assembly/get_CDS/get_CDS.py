
'''

Get the CDS of all transcripts, and measure it against the GENCDOE for accuracy.

As the GTF is not guaranteed to match the genocode CDS locations, I need to generate
themselves for all transcripts, including the

'''

import glob, sys, os, gzip, numpy, math, re
import matplotlib.pyplot as plot
from glbase3 import glload, utils, expression, genelist, genome_sql

def find_cds(seq):
    '''
    Return the most likely CDS
    '''
    def split3(s):
        return ' '.join([s[i:i+3] for i in range(0, len(s), 3)])

    def in_same_frame(s, e):
        if e < s:
            return False
        l = e - s
        return l % 3 == 0

    # get all CDS:

    start = re.compile('ATG', re.IGNORECASE)
    ends = re.compile('TAA|TAG|TGA', re.IGNORECASE)

    all_starts = start.finditer(seq)
    all_stops = ends.finditer(seq)

    all_starts = [i.span()[0] for i in all_starts]
    all_ends = [i.span()[1] for i in all_stops]

    # All pairs of start/end;
    longest_transcripts = {}
    for s in all_starts:
        temp_orfs = {}
        for e in all_ends:
            # check they are in the same triplet frame
            if in_same_frame(s, e):
                temp_orfs[e-s] = (s, e)
        # get the shortest ATG -> STOP
        shortest_orf = sorted(temp_orfs.keys())[0]
        longest_transcripts[shortest_orf] = temp_orfs[shortest_orf]

    # I now need to check that
    print(longest_transcripts)
    for p in longest_transcripts.values():
        print(p, split3(seq[p[0]:p[1]]))

    1/0


    return cdsl, cdsr

#gencode = glload('../te_transcripts/transcript_table_gencode_pc.glb')
#gencode_sliced = gencode.getColumns(['cds_loc', 'transcript_id', 'loc'])
#gencode_sliced = gencode_sliced.renameKey('transcript_id', 'enst')

fastas = glload('../../transcript_assembly/fasta/transcripts.glb')

newl = []
for f in fastas:
    if f['coding'] == 'noncoding':
        continue
    if f['tags'][-1] == '~': # variant sequence
        cdsl, cdsr = find_cds(f['seq'])
        newl.append(f)
    elif f['tags'][-1] == '=': # can get this one from GENCODE
        pass
