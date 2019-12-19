
'''

Build the final annotation tables;

'''


import glob, sys, os, gzip
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6

sys.path.append('../')
import pies
sys.path.append('../../../')
import shared
from glbase3 import glload, utils, expression, genelist

all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
te = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
not_te = te.map(genelist=all_genes, key='transcript_id', logic='notright')

# collect stats:
res = {
    'pc-all-TE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'pc-known-TE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'pc-variant-TE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    #'pc-unknown-TE': {'ES+': 0, 'ES:': 0, 'ES-': 0},

    'pc-all-noTE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'pc-known-noTE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'pc-variant-noTE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    #'pc-unknown-noTE': {'ES+': 0, 'ES:': 0, 'ES-': 0},

    'ncrna-all-TE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'ncrna-known-TE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'ncrna-variant-TE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'ncrna-unknown-TE': {'ES+': 0, 'ES:': 0, 'ES-': 0},

    'ncrna-all-noTE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'ncrna-known-noTE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'ncrna-variant-noTE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    'ncrna-unknown-noTE': {'ES+': 0, 'ES:': 0, 'ES-': 0},
    }

data = {'TE': te, 'noTE': not_te}
expn_key_map = {'enriched': 'ES+', 'unbiased': 'ES:', 'depleted': 'ES-'}


for k in data:
    for g in data[k]:
        expn = expn_key_map[g['expression']]

        if ';C;' in g['name']:
            res['pc-all-{0}'.format(k)][expn] += 1

            if ';~' in g['name']:
                res['pc-variant-{0}'.format(k)][expn] += 1
            elif ';=)' in g['name']:
                res['pc-known-{0}'.format(k)][expn] += 1
            #elif ';!' in g['name']:
            #    res['pc-unknown-{0}'.format(k)][expn] += 1

        elif ';NC;' in g['name']:
            res['ncrna-all-{0}'.format(k)][expn] += 1
            if ';=' in g['name']:
                res['ncrna-known-{0}'.format(k)][expn] += 1
            elif ';~' in g['name']:
                res['ncrna-variant-{0}'.format(k)][expn] += 1
            elif ';!' in g['name']:
                res['ncrna-unknown-{0}'.format(k)][expn] += 1

title_map = {'pc-te': 'protein-coding with TE',
    'pc-noTE': 'protein-coding without TE',
    }

for k in res:
    if k in title_map:
        pies.pie('pies/te_%s.png' % k, [res[k]['ES+'], res[k]['ES:'], res[k]['ES-']], ['ES+', 'ES:', 'ES-'], title_map[k])

pies.split_bar('bar.png'.format(k), res, cols=['#d62728', '#2ca02c', '#ff7f0e', ])

# pickle the results
import pickle
shared.pickle_it('results.pickle', res)
