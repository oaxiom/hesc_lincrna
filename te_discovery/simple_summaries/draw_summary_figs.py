
'''

Build the final annotation tables;

'''


import glob, sys, os, gzip
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6

sys.path.append('../../')
import shared
from glbase3 import glload, utils, expression, genelist

gtf = shared.get_pickle('current_gtf/results.pickle')
gencode = shared.get_pickle('gencode/results.pickle')
expn = shared.get_pickle('es_expressed/results.pickle')
expn2 = shared.get_pickle('es_expressed-te-by-expn/results.pickle')

print(gtf)
print(gencode)

data = {
    'hPSC-novel': gtf['novel_ncrna'],
    'hPSC-variant': gtf['ncrna_variant'],
    'hPSC-matching': gtf['ncrna_matching'],
    'GENCODE': gencode['ncrna']
    }
shared.split_bar('bar-ncrna.png', data, key_order=['TE', 'nonTE'], cols=['#ff7f0e', '#1f77b4']) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])


data = {
    'hPSC-variant': gtf['pc_variant'],
    'hPSC-matching': gtf['pc_matching'],
    'GENCODE': gencode['pc']
    }
shared.split_bar('bar-pc.png', data, key_order=['TE', 'nonTE'], cols=['#ff7f0e', '#1f77b4']) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

data = {
    'Variant with TE': expn['pc-variant-TE'],
    'Variant no TE': expn['pc-variant-noTE'],

    'Matching with TE': expn['pc-known-TE'],
    'Matching no TE': expn['pc-known-noTE'],

    #'pc-all-TE': expn['pc-all-TE'],
    #'pc-all-noTE': expn['pc-all-noTE'],
    }

shared.split_bar('bar-expn-pc.png', data, key_order=['ES+', 'ES:', 'ES-'], cols=['#d62728', '#2ca02c', '#ff7f0e', ]) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

data = {
    'hPSC-depleted': expn2['pc-ES-'],
    'hPSC-unbiased': expn2['pc-ES:'],
    'hPSC-enriched': expn2['pc-ES+'],
    }

shared.split_bar('bar-te-by-expn-pc.png', data, key_order=['TE', 'nonTE'], cols=['#ff7f0e', '#1f77b4' ]) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

# And the ncrna versions:
data = {
    'Novel with TE': expn['ncrna-unknown-TE'],
    'Novel no TE': expn['ncrna-unknown-noTE'],

    'Variant with TE': expn['ncrna-variant-TE'],
    'Variant no TE': expn['ncrna-variant-noTE'],

    'Known with TE': expn['ncrna-known-TE'],
    'Known no TE': expn['ncrna-known-noTE'],

    #'ncrna-all-TE': expn['ncrna-all-TE'],
    #'ncrna-all-noTE': expn['ncrna-all-noTE'],
    }

shared.split_bar('bar-expn-ncrna.png', data, key_order=['ES+', 'ES:', 'ES-'], cols=['#d62728', '#2ca02c', '#ff7f0e', ]) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

data = {
    'hPSC-depleted': expn2['ncrna-ES-'],
    'hPSC-unbiased': expn2['ncrna-ES:'],
    'hPSC-enriched': expn2['ncrna-ES+'],
    }

shared.split_bar('bar-te-by-expn-ncrna.png', data, key_order=['TE', 'nonTE'], cols=['#ff7f0e', '#1f77b4' ]) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
