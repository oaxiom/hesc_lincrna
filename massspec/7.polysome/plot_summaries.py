'''

For the masked peptides, see if we can find them in the MS data;

'''

import glob, sys, os, numpy
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../')
import shared
config.draw_mode = 'pdf'

res = {}

all_te_transcripts = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')

has_peptide_hit = glload('../results_gene.glb').removeDuplicates('transcript_id') # ones with a peptide hit.
#has_peptide_hit = has_peptide_hit.map(genelist=all_te_transcripts, key='transcript_id')

no_peptide_hit = glload('../2.blast_searches/super_table.glb').removeDuplicates('transcript_id')
no_peptide_hit = has_peptide_hit.map(genelist=no_peptide_hit, key='transcript_id', logic='notright')

#polysom = glload('../../polysome_rnaseq/rsem_star/polysome_index.glb')
polysom = glload('../../polysome_rnaseq/pack_kallisto/polysome_index.glb')
print(polysom)
#polysom.log(2, .1)
cond_names = polysom.getConditionNames()

res_c = {
    'hit-massspec': {c: [] for c in cond_names},
    'nohit-massspec': {c: [] for c in cond_names},
    }

todo = {'hit': has_peptide_hit,
    'nohit': no_peptide_hit}

for k in todo:
    todo[k] = todo[k].map(genelist=polysom, key='transcript_id')
    todo[k] = todo[k].removeDuplicates('transcript_id')
    for transcript in todo[k]:
        for idx, cond in enumerate(polysom.getConditionNames()):
            res_c['{}-massspec'.format(k)][cond].append(transcript['conditions'][idx])

    shared.boxplots_simple('polysome_ratios-{}.pdf'.format(k), res_c['{}-massspec'.format(k)], None, xlims=[0, 1], col='#FF8A87', vlines=[0.7421478187528863, 0.2960685357901553])

    todo[k] = todo[k].map(genelist=polysom, key='transcript_id')

    todo[k].heatmap(filename='heat-{}.pdf'.format(k), heat_wid=0.1, heat_hei=0.5, bracket=[0, 1])


