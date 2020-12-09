
'''

Build the final annotation tables;

'''

import glob, sys, os, gzip, numpy
from glbase3 import glload, utils, expression, genelist, draw, config

sys.path.append('../../')
import shared

config.draw_mode = 'pdf'

all_genes = glload('../../transcript_assembly/packed/all_genes.glb')
sam = glload('../hmmer_dfam/data/transcript_table_HSC_SR_PB_merged.transcripts.glb')
sam = sam.map(genelist=all_genes, key='transcript_id')
print(sam)

dfam = genelist('../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4, 'len': 1})

dfam = {i['name']: {'full_name': '{0}:{1}:{2}'.format(i['type'], i['subtype'], i['name']), 'type': i['type'], 'len': i['len']} for i in dfam}

dom_lens = {}

completion ={
    'coding': {'LTR': {}, 'SINE': {}, 'LINE': {}, 'Retroposon': {},},
    'noncoding': {'LTR': {}, 'SINE': {}, 'LINE': {}, 'Retroposon': {},},
    }

keep_types = set(['LTR', 'SINE', 'LINE', 'Retroposon'])

for match in sam:
    for dom in match['doms']:
        dfam_entry = dfam[dom['dom']]
        te_type = dfam_entry['type']

        if te_type not in keep_types:
            continue

        if match['coding'] == 'NA':
            print(match)
            continue

        if dfam_entry['full_name'] not in completion[match['coding']][te_type]:
            completion[match['coding']][te_type][dfam_entry['full_name']] = []

        dom_len = dom['span'][1] - dom['span'][0]

        completion[match['coding']][te_type][dfam_entry['full_name']].append(dom_len / dfam_entry['len']*100.0)

gldraw = draw()

heat_data = {
    'coding': {'LTR': {}, 'SINE': {}, 'LINE': {}, 'Retroposon': {},},
    'noncoding': {'LTR': {}, 'SINE': {}, 'LINE': {}, 'Retroposon': {},},
    }


heat_data = {}
for ti, coding in enumerate(['coding', 'noncoding']):
    for te_type in completion[coding]:
        labs = []
        vals = []
        for TE in sorted(completion[coding][te_type].keys()):
            if len(completion[coding][te_type][TE]) > 10:
                labs.append(TE)
                vals.append(completion[coding][te_type][TE])

            if TE not in heat_data:
                heat_data[TE] = {'name': TE, 'conditions': [0,0]}
            heat_data[TE]['conditions'][ti] = numpy.mean(completion[coding][te_type][TE])

        shared.boxplots(filename='box_{0}_TEs_vert-{1}.pdf'.format(coding, te_type),
            data=dict(zip(labs, vals)),
            xlims=[-5, 105],
            no_TE_key=None,
            title=None,
            trim_low_samples=False,
            sizer=0.005,
            vert_height=18,
            bot_pad=0.01,
            qs='100%quartiles')

    print(heat_data)

e = expression()
e.load_list(list(heat_data.values()), cond_names=['coding', 'noncoding'])

e.heatmap('TEs.pdf',
    heat_wid=0.05,
    digitize=6,
    bracket=[0, 100])
