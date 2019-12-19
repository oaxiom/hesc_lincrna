
'''

Build the final annotation tables;

'''


import glob, sys, os, gzip
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6

import pies
sys.path.append('../../')
import shared
from glbase3 import glload, utils, expression, genelist

gtf = shared.get_pickle('current_gtf/results.pickle')
gencode = shared.get_pickle('gencode/results.pickle')
expn = shared.get_pickle('es_expressed/results.pickle')

print(gtf)
print(gencode)

data = {
    'hPSC-novel': gtf['novel_ncrna'],
    'hPSC-variant': gtf['ncrna_variant'],
    'hPSC-known': gtf['ncrna'],
    'GENCODE': gencode['ncrna']
    }
pies.split_bar('bar-ncrna.png', data, key_order=['TE', 'nonTE'], cols=['#ff7f0e', '#1f77b4']) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])


data = {
    'hPSC-variant': gtf['pc_variant'],
    'hPSC-known': gtf['pc'],
    'GENCODE': gencode['pc']
    }
pies.split_bar('bar-pc.png', data, key_order=['TE', 'nonTE'], cols=['#ff7f0e', '#1f77b4']) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

data = {
    'pc-variant-TE': expn['pc-variant-TE'],
    'pc-variant-noTE': expn['pc-variant-noTE'],

    'pc-known-TE': expn['pc-known-TE'],
    'pc-known-noTE': expn['pc-known-noTE'],

    'pc-all-TE': expn['pc-all-TE'],
    'pc-all-noTE': expn['pc-all-noTE'],
    }
pies.split_bar('bar-expn-pc.png', data, key_order=['ES+', 'ES:', 'ES-'], cols=['#d62728', '#2ca02c', '#ff7f0e', ]) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

