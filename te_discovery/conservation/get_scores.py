import sys, os, itertools
import matplotlib.pyplot as plot
from glbase3 import *

# collect three things:
# 1. The total PhyloP score of the transcript
# 2. score for TE-containing bits
# 3. score for non-TE containign bits;

transcripts = glload('../te_transcripts/transcript_table_merged.mapped.glb')

phylo = flat_track(filename=os.path.expanduser('~/hg38/goldenPath/hg38.phyloP17way.flat'))

print(transcripts)

res = []

setter = itertools.count(-50)

for t in transcripts:
    # just extract all the per-base exon scores; then distribute them to the TE or not-TE categories

    all = 0.0
    TE = 0.0
    notTE = 0.0

    all_scores = phylo.get(t['loc'])

    new_scores = []

    # get only the exons;
    for es, ee in zip(t['exonStarts'], t['exonEnds']):
        ls = es - t['loc']['left']
        le = ee - t['loc']['left']+1
        new_scores += list(all_scores[ls:le])

    if t['strand'] == '-':
        new_scores = new_scores[::-1]

    te_scores = []
    for te in t['doms']:
        te_scores += new_scores[te['span'][0]:te['span'][1]]

    # To get the non-TE parts, I think it's just easiest to set all non-TE parts to an invalid value
    # and then trim them.
    for te in t['doms']:
        for i in range(te['span'][0], te['span'][1]-1):
            new_scores[i] = -50

    notte_scores = [i for i in new_scores if i > -50]


    result = {'name': t['name'],
        'transcript_id': t['transcript_id'],
        'phyloP_all': sum(new_scores) / len(new_scores),
        'phyloP_tes': sum(te_scores) / len(te_scores),
        'phyloP_nottes': sum(notte_scores) / len(notte_scores),
        'TPM': t['TPM'],
        }

    res.append(result)

gl = genelist()
gl.load_list(res)
gl.saveTSV('phyloP_conservation_table.tsv')
gl.save('phyloP_conservation_table.glb')
