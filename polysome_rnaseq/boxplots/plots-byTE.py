'''

boxplots for polysome, etc.

'''

import sys, os, itertools
import numpy as np
import matplotlib.tri as tri
from glbase3 import *

sys.path.append('../')
import shared_polysome

expn = glload('../pack_kallisto/kall_ratios.glb')
#expn.log(2,.1)
cond_names = expn.getConditionNames()

all_genes = glload('../../transcript_assembly/packed/all_genes.glb')
#dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

contains_te = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
contains_not_te = contains_te.map(genelist=all_genes, key='transcript_id', logic='notright')


'''
dfam_dict = {}
for te in dfam:
    dfam_dict[te['name']] = '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name'])
'''

todo = {'TE': contains_te, 'no TE': contains_not_te}

res_c = {
    'pc-TE': {c: [] for c in cond_names},
    'ncrna-TE': {c: [] for c in cond_names},
    'pc-no TE': {c: [] for c in cond_names},
    'ncrna-no TE': {c: [] for c in cond_names},
    }

for idx, cond in enumerate(cond_names):
    for k in todo:
        mapped = todo[k].map(key='transcript_id', genelist=expn)
        for transcript in mapped:
            if ';C;' in transcript['name']:
                 res_c['pc-{}'.format(k)][cond].append(transcript['conditions'][idx])

            elif ';NC;' in transcript['name']:
                res_c['ncrna-{}'.format(k)][cond].append(transcript['conditions'][idx])

for k in res_c:
    shared_polysome.boxplots('{}.pdf'.format(k), res_c[k], None, xlims=[0, 1], col='#FF8A87')



