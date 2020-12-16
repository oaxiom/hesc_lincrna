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
        all_scores = phylo.get(t['loc'])

        new_scores = []

        # get only the exons;
        for es, ee in zip(t['exonStarts'], t['exonEnds']):
            ls = es - t['loc']['left']
            le = ee - t['loc']['left']+1
            new_scores += list(all_scores[ls:le])
        phyloP_all = sum(new_scores) / len(new_scores)

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
