
'''

Measure the expression of TEs by family and then by class, plot a violin and a heatmap;

'''

import glob, sys, os, gzip, numpy, math, scipy.stats
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6
from collections import defaultdict
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from glbase3 import glload, utils, expression, genelist, genome_sql, draw, config
config.draw_mode = ['pdf']

sys.path.append('../../../')
import shared

draw_type = 'pdf'

dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
contains_te = glload('../../te_transcripts/transcript_table_merged.mapped.glb')

ages = glload('../../dfam/te_ages.glb')
ages = {i['TE']: i['age'] for i in ages} # lookup;

dfam_dict = {}
for te in dfam:
    dfam_dict[te['name']] = '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name'])

def class_dict():
    return {
        'pc-all': {'TE': [], 'nonTE': []},
        'ncrna-all': {'TE': [], 'nonTE': []},

        'pc-known': {'TE': [], 'nonTE': []},
        'pc-variant': {'TE': [], 'nonTE': []},

        'ncrna-known': {'TE': [], 'nonTE': []},
        'ncrna-variant': {'TE': [], 'nonTE': []},
        'ncrna-unknown': {'TE': [], 'nonTE': []},
        }

def dict_builder():
    return {k: [] for k in class_dict()}

#Below: Split by te_family;
res_type = defaultdict(dict_builder)
gldraw = draw()

for g in contains_te:
    for d in g['doms']:
        tpm = math.log2(g['TPM']+0.1)
        TE = d['dom']
        full_name = dfam_dict[TE]

        if g['coding'] == 'coding':
            res_type[full_name]['pc-all'].append(tpm)

            if ';~' in g['name']:
                res_type[full_name]['pc-variant'].append(tpm)
            elif ';=)' in g['name']:
                res_type[full_name]['pc-known'].append(tpm)

        elif g['coding'] == 'noncoding':
            res_type[full_name]['ncrna-all'].append(tpm)

            if ';=' in g['name']:
                res_type[full_name]['ncrna-known'].append(tpm)
            elif ';~' in g['name']:
                res_type[full_name]['ncrna-variant'].append(tpm)
            elif ';!' in g['name']:
                res_type[full_name]['ncrna-unknown'].append(tpm)

p_scatter = {
    'pc-all': [],
    'ncrna-all': [],

    'pc-known': [],
    'pc-variant': [],

    'ncrna-known': [],
    'ncrna-variant': [],
    'ncrna-unknown': [],
    }

for te in sorted(res_type):
    for t in class_dict().keys(): # pc-all, ncrna-all' etc.
        if True in [typ in te for typ in ['SINE', 'LINE', 'LTR', 'Retroposon']]:
            data = res_type[te][t]

            if len(data) <= 20:
                continue

            if te in ages:
                age = ages[te]
                print(data)
                avg_expn = 2**numpy.mean(data)
                p_scatter[t].append({'name': te, 'age': age, 'exp': avg_expn, 'n': len(data)})

# adjust the p values pls
#print(p_scatter)

for t in class_dict().keys():
    gl = genelist()
    gl.load_list(p_scatter[t])

    spot_cols = []
    for name in gl['name']:
        if 'LINE' in name:
            spot_cols.append('blue')
        elif 'LTR' in name:
            spot_cols.append('red')
        elif 'SINE' in name:
            spot_cols.append('green')
        else:
            spot_cols.append('grey')

    config.draw_mode = 'svg'
    shared.nice_scatter(x=gl['age'], y=gl['exp'], figsize=[2,2], spot_size=12,
        spot_cols=spot_cols, label_t=0.6,
        filename='MA-{0}.png'.format(t), label=gl['name'], hlines=[0.6])
    config.draw_mode = 'pdf'
