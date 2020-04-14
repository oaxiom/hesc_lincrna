import sys, os, itertools
from collections import defaultdict
import matplotlib.pyplot as plot
from glbase3 import *

# collect three things:
# 1. The total PhyloP score of the transcript
# 2. score for TE-containing bits
# 3. score for non-TE containign bits;

transcripts = glload('../te_transcripts/transcript_table_merged.mapped.glb')

phylo = flat_track(filename=os.path.expanduser('~/hg38/goldenPath/hg38.phyloP17way.flat'))

print(transcripts)

res = {}

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

    te_scores = defaultdict(list)
    for te in t['doms']:
        te_scores[te['dom']] += new_scores[te['span'][0]:te['span'][1]]

    # To get the non-TE parts, I think it's just easiest to set all TE parts to an invalid value
    # and then trim them.
    for te in t['doms']:
        for i in range(te['span'][0], te['span'][1]-1):
            new_scores[i] = -50

    notte_scores = [i for i in new_scores if i > -50]

    # Append score per TE type:
    for te in t['doms']:
        if te['dom'] not in res:
            res[te['dom']] = {'phyloP_tes': [], 'phyloP_nottes': [], 'TPM': []}

            res[te['dom']]['phyloP_tes'].append(sum(te_scores[te['dom']]) / len(te_scores[te['dom']]))
            res[te['dom']]['phyloP_nottes'].append(sum(notte_scores) / len(notte_scores))
            res[te['dom']]['TPM'].append(t['TPM'])

newres = []
for te in res:
    result = {
        'TE': te,
        'phyloP_tes': sum(res[te]['phyloP_tes']) / len(res[te]['phyloP_tes']),
        'phyloP_nottes': sum(res[te]['phyloP_nottes']) / len(res[te]['phyloP_nottes']),
        'TPM': sum(res[te]['TPM']) / len(res[te]['TPM']),
        }
    newres.append(result)

gl = genelist()
gl.load_list(newres)
gl.sort('TE')
gl.saveTSV('phyloP_conservation_table_per_TE_type.tsv')
gl.save('phyloP_conservation_table_per_TE_type.glb')
