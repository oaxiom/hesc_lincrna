import sys, os, itertools
import matplotlib.pyplot as plot
from glbase3 import *

# collect three things:
# 1. The total PhyloP score of the transcript
# 2. score for TE-containing bits
# 3. score for non-TE containign bits;


all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
contains_te = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
contains_not_te = contains_te.map(genelist=all_genes, key='transcript_id', logic='notright')

phylo = flat_track(filename=os.path.expanduser('~/hg38/goldenPath/hg38.phyloP17way.flat'))

res = []

todo = {'contains_te': contains_te, 'contains_not_te': contains_not_te}

for k in todo:
    print(k)
    # just extract all the per-base exon scores; then distribute them to the TE or not-TE categories

    for t in todo[k]:

        # get the tss:
        if t['strand'] == '+':
            tss_loc = t['loc'].pointLeft().expandLeft(1000).expandRight(100)
        elif t['strand'] == '-':
            tss_loc = t['loc'].pointRight().expandRight(1000).expandLeft(100)

        #print(t)
        #print(t['strand'])
        #print(t['loc'], tss_loc)

        all_scores = phylo.get(tss_loc)

        phyloP_all = sum(all_scores) / len(all_scores)

        result = {'name': t['name'],
            'transcript_id': t['transcript_id'],
            'phyloP_all': phyloP_all,
            'TPM': t['TPM'],
            'strand': t['strand'],
            }

        res.append(result)

    gl = genelist()
    gl.load_list(res)
    gl.saveTSV('phyloP_conservation_table-{}.tsv'.format(k))
    gl.save('phyloP_conservation_table-{}.glb'.format(k))
