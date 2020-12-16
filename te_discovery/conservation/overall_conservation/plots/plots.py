'''

Conservation boxplots for all ;

'''

import sys, os, itertools
import numpy as np
import matplotlib.tri as tri
from glbase3 import *

sys.path.append('../../')
import shared_conservation


todo = {'TE': glload('../phyloP_conservation_table-contains_te.glb'),
    'no TE': glload('../phyloP_conservation_table-contains_not_te.glb')}

res_c = {
    'pc':
        {'TE-Variant': [], 'TE-Matching': [], 'no TE-Variant': [], 'no TE-Matching': []},
    'ncrna':
        {'TE-Novel': [], 'TE-Variant': [], 'TE-Matching': [], 'no TE-Novel': [], 'no TE-Variant': [], 'no TE-Matching': []}}

res_es = {
    'pc':
        {'TE-Depleted': [], 'TE-Nonspecific': [],  'TE-Enriched': [], 'no TE-Depleted': [], 'no TE-Nonspecific': [],  'no TE-Enriched': [], },
    'ncrna':
        {'TE-Depleted': [], 'TE-Nonspecific': [],  'TE-Enriched': [], 'no TE-Depleted': [], 'no TE-Nonspecific': [],  'no TE-Enriched': [], }}

for k in todo:
    for transcript in todo[k]:
        if ';C;' in transcript['name']:
            if ';=' in transcript['name']:
                res_c['pc']['{}-Matching'.format(k)].append(transcript['phyloP_all'])
            elif ';~' in transcript['name']:
                res_c['pc']['{}-Variant'.format(k)].append(transcript['phyloP_all'])
            #elif ';!' in transcript['name']:
            #    res_c['pc']['Novel'].append(transcript['phyloP_all'])

            if 'ES+' in transcript['name']:
                res_es['pc']['{}-Enriched'.format(k)].append(transcript['phyloP_all'])
            elif 'ES:' in transcript['name']:
                res_es['pc']['{}-Nonspecific'.format(k)].append(transcript['phyloP_all'])
            elif 'ES-' in transcript['name']:
                res_es['pc']['{}-Depleted'.format(k)].append(transcript['phyloP_all'])

        elif ';NC;' in transcript['name']:
            if ';=' in transcript['name']:
                res_c['ncrna']['{}-Matching'.format(k)].append(transcript['phyloP_all'])
            elif ';~' in transcript['name']:
                res_c['ncrna']['{}-Variant'.format(k)].append(transcript['phyloP_all'])
            elif ';!' in transcript['name']:
                res_c['ncrna']['{}-Novel'.format(k)].append(transcript['phyloP_all'])

            if 'ES+' in transcript['name']:
                res_es['ncrna']['{}-Enriched'.format(k)].append(transcript['phyloP_all'])
            elif 'ES:' in transcript['name']:
                res_es['ncrna']['{}-Nonspecific'.format(k)].append(transcript['phyloP_all'])
            elif 'ES-' in transcript['name']:
                res_es['ncrna']['{}-Depleted'.format(k)].append(transcript['phyloP_all'])


shared_conservation.boxplots('pc.pdf', res_c['pc'], None, xlims=[-0.27, 0.52], col='#FF8A87')
shared_conservation.boxplots('ncrna.pdf', res_c['ncrna'], None, xlims=[-0.27, 0.52], col='#92A7FF')

shared_conservation.boxplots('enrich-pc.pdf', res_es['pc'], None, xlims=[-0.27, 0.52], col='#FF8A87')
shared_conservation.boxplots('enrich-ncrna.pdf', res_es['ncrna'], None, xlims=[-0.27, 0.52], col='#92A7FF')

# scatters;
'''
shared_conservation.hist('hist_all_vs_tes.pdf',
    gl['phyloP_all'], gl['phyloP_tes'],
    'All', 'not-TE',
    ranges=[[-0.2, 0.7], [-0.2, 0.7]],
    hlines = [0, 0.25],
    vlines = [0, 0.25],
    )

shared_conservation.hist('hist_all_vs_nottes.pdf',
    gl['phyloP_all'], gl['phyloP_nottes'],
    'All', 'not-TE',
    ranges=[[-0.2, 0.7], [-0.2, 0.7]],
    hlines = [0, 0.25],
    vlines = [0, 0.25],
    )

for t in ('phyloP_all', 'phyloP_tes', 'phyloP_nottes'):
    shared_conservation.hist(filename='hist_expn_cons_vs_{0}.pdf'.format(t),
        x=gl[t],
        y=np.log2(np.array(gl['TPM'])+0.1),
        xlabel='cons',
        ylabel='expn',
        ranges=[[-0.2, 0.7], [-3, 9]],
        )
'''
