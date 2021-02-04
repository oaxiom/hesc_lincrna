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

#localise_ric = glload('../../polysome_rnaseq/rsem_star/polysome_index.glb')
#localise_ric = glload('../../polysome_rnaseq/pack_kallisto/kall_nuc_cyt_ratios_RIC.glb')
#print(localise_ric)

# USe the table from elsewhere in the manuscript;
localise_ric = expression(filename='nucleus_cytosol_0.001_TE.tsv.gz',
    gzip=True,
    format={'force_tsv': True, 'transcript_id': 0, 'skiplines': -1},
    cond_names=['RIC'],
    expn='[column[2],]')

print(localise_ric)

#polysom.log(2, .1)
cond_names = localise_ric.getConditionNames()

res_c = {
    'hit-massspec': {c: [] for c in cond_names},
    'nohit-massspec': {c: [] for c in cond_names},
    }

todo = {'hit': has_peptide_hit,
    'nohit': no_peptide_hit}

for k in todo:
    print(todo[k])
    todo[k] = todo[k].map(genelist=localise_ric, key='transcript_id')
    #todo[k] = todo[k].removeDuplicates('transcript_id')
    for transcript in todo[k]:
        for idx, cond in enumerate(localise_ric.getConditionNames()):
            res_c['{}-massspec'.format(k)][cond].append(transcript['conditions'][idx])

    shared.boxplots_simple('localisation-{}.pdf'.format(k), res_c['{}-massspec'.format(k)], None, xlims=[-4, 4], col='#FF8A87')

    todo[k].sort_sum_expression()

    todo[k].heatmap(filename='heat-{}.pdf'.format(k), heat_wid=0.02, heat_hei=0.5,
        #bracket=[-20, 20],
        row_cluster=False)


