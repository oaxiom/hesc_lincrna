
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
    all_ends.append(len(seq)) # Sometimes not STOP, as possible truncated, so add a hypothetical STOP in all three codon frames
    all_ends.append(len(seq)-1)
    all_ends.append(len(seq)-2)

    if not all_starts: # There is no ATG, so best to just make no prediction
        return False, -1, -1

    # All pairs of start/end;
    longest_transcripts = {}
    for s in all_starts:
        temp_orfs = {}
        for e in all_ends:
            # check they are in the same triplet frame
            if in_same_frame(s, e):
                temp_orfs[e-s] = (s, e)
        # get the shortest ATG -> STOP
        if temp_orfs:
            shortest_orf = sorted(temp_orfs.keys())[0]
            longest_transcripts[shortest_orf] = temp_orfs[shortest_orf]

    lengths = sorted(list(longest_transcripts.keys()), reverse=True)
    print(lengths)
    cdsl, cdsr = longest_transcripts[lengths[0]][0], longest_transcripts[lengths[0]][1]
    print(lengths[0], split3(seq[cdsl:cdsr]))

    # best guess is the longest intact transcript
    return True, longest_transcripts[lengths[0]][0], longest_transcripts[lengths[0]][1]

#gencode = glload('../te_transcripts/transcript_table_gencode_pc.glb')
#gencode_sliced = gencode.getColumns(['cds_loc', 'transcript_id', 'loc'])
#gencode_sliced = gencode_sliced.renameKey('transcript_id', 'enst')

fastas = glload('../../transcript_assembly/fasta/transcripts.glb')

no_prediction = 0
newl = []
for f in fastas:
    print(f['enst'])
    if f['coding'] == 'noncoding':
        continue
    if f['tags'][-1] == '~': # variant sequence
        status, cdsl, cdsr = find_cds(f['seq'])
        newl.append(f)
    elif f['tags'][-1] == '=': # can get this one from GENCODE
        status, cdsl, cdsr = find_cds(f['seq'])

    if status:
        f['cds_local_locs'] = (cdsl, cdsr)
        del f['seq']
        newl.append(f)
    else:
        no_prediction += 1
        print('Cound not predict for {0}'.format(f['name']))

print('Cound not make a ORF prediction for {0}'.format(no_prediction))

newd = genelist()
newd.load_list(newl)
newd.save('coding_genes_with_local_CDS.glb')
newd.saveTSV('coding_genes_with_local_CDS.tsv')
