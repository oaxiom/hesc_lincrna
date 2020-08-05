
'''

Build the final annotation tables;

'''

import glob, sys, os, gzip
from glbase3 import glload, utils, expression, genelist, draw, config

sys.path.append('../../')
import shared

config.draw_mode = 'pdf'

sam = glload('../hmmer_dfam/data/transcript_table_HSC_SR_PB_merged.transcripts.glb')

dfam = genelist('../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4, 'len': 1})

dfam = {i['name']: {'full_name': '{0}:{1}:{2}'.format(i['type'], i['subtype'], i['name']), 'type': i['type'], 'len': i['len']} for i in dfam}

dom_lens = {}

completion = {'LTR': {},
    'SINE': {},
    'LINE': {},
    'Retroposon': {},
    }

keep_types = set(['LTR', 'SINE', 'LINE', 'Retroposon'])

for match in sam:
    for dom in match['doms']:
        dfam_entry = dfam[dom['dom']]
        te_type = dfam_entry['type']

        if te_type not in keep_types:
            continue

        if dfam_entry['full_name'] not in completion[te_type]:
            completion[te_type][dfam_entry['full_name']] = []

        dom_len = dom['span'][1] - dom['span'][0]

        completion[te_type][dfam_entry['full_name']].append(dom_len / dfam_entry['len']*100.0)

gldraw = draw()

for t in completion:
    labs = []
    vals = []
    for k in sorted(completion[t].keys()):
        if len(completion[t][k]) > 50:
            labs.append(k)
            vals.append(completion[t][k])

    shared.boxplots(filename='box_all_TEs_vert-{0}.pdf'.format(t),
        data=dict(zip(labs, vals)),
        xlims=[-5, 105],
        no_TE_key=None,
        title=None,
        trim_low_samples=False,
        sizer=0.005,
        vert_height=18,
        bot_pad=0.01,
        qs='100%quartiles')

gldraw.boxplot(vals,
    labels=labs,
    size=[24,4],
    ylims=[0,100],
    showfliers=False,
    filename='box_all_TEs.pdf')
