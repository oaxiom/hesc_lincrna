
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

draw_type = 'png'

all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
contains_te = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
contains_not_te = contains_te.map(genelist=all_genes, key='transcript_id', logic='notright')
contains_not_te.saveTSV('notTE.tsv')
dfam_dict = {}
for te in dfam:
    dfam_dict[te['name']] = '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name'])

def class_dict():
    return {
        'all': {'TE': [], 'nonTE': []},
        'pc-all': {'TE': [], 'nonTE': []},
        'ncrna-all': {'TE': [], 'nonTE': []},

        'pc-known': {'TE': [], 'nonTE': []},
        'pc-variant': {'TE': [], 'nonTE': []},

        'ncrna-known': {'TE': [], 'nonTE': []},
        'ncrna-variant': {'TE': [], 'nonTE': []},
        'ncrna-unknown': {'TE': [], 'nonTE': []},
        }

data = {'TE': contains_te, 'nonTE': contains_not_te}

res = class_dict()

for datatype in data:
    for g in data[datatype]:
        tpm = math.log2(g['TPM']+0.1)
        res['all'][datatype].append(tpm)

        if g['coding'] == 'coding':
            res['pc-all'][datatype].append(tpm)

            if ';~' in g['name']:
                res['pc-variant'][datatype].append(tpm)
            elif ';=)' in g['name']:
                res['pc-known'][datatype].append(tpm)

        elif g['coding'] == 'noncoding':
            res['ncrna-all'][datatype].append(tpm)

            if ';=' in g['name']:
                res['ncrna-known'][datatype].append(tpm)
            elif ';~' in g['name']:
                res['ncrna-variant'][datatype].append(tpm)
            elif ';!' in g['name']:
                res['ncrna-unknown'][datatype].append(tpm)


gldraw = draw()
for typ in res:
    gldraw.beanplot(filename='viol_{0}.png'.format(typ), data=res[typ], figsize=[2,1.8],
        beans=False, ylims=[-2.5, 8])

def dict_builder():
    return {k: [] for k in class_dict()}

#Below: Split by te_type;
res_type = defaultdict(dict_builder)
gldraw = draw()

for g in contains_te:
    for d in g['doms']:
        tpm = math.log2(g['TPM']+0.1)
        TE = d['dom']
        full_name = dfam_dict[TE]
        tetype = full_name.split(':')[0]

        res_type[tetype]['all'].append(tpm)

        if g['coding'] == 'coding':
            res_type[tetype]['pc-all'].append(tpm)
        elif g['coding'] == 'noncoding':
            res_type[tetype]['ncrna-all'].append(tpm)

for te in res_type:
    print(te)
    for t in class_dict().keys(): # pc-all, ncrna-all' etc.
        if te in ['SINE', 'LINE', 'LTR', 'Retroposon']:
            data = {te: res_type[te][t], 'noTE': res[t]['nonTE']}

            p = scipy.stats.mannwhitneyu()

            gldraw.beanplot(filename='by_type/viol_{0}-{1}.png'.format(te, t), data=data,
                figsize=[2,1.8], beans=False, ylims=[-4, 8])
            gldraw.boxplot(filename='by_type/box_{0}-{1}.png'.format(te, t), data=data.values(), labels=data.keys(),
                figsize=[2,1.8], ylims=[-4, 8])


res_type = defaultdict(dict_builder)
gldraw = draw()

for g in contains_te:
    for d in g['doms']:
        tpm = math.log2(g['TPM']+0.1)
        TE = d['dom']
        full_name = dfam_dict[TE]
        tesubtype = ':'.join(full_name.split(':')[0:2])

        res_type[tesubtype]['all'].append(tpm)

        if g['coding'] == 'coding':
            res_type[tesubtype]['pc-all'].append(tpm)

            if ';~' in g['name']:
                res_type[tesubtype]['pc-variant'].append(tpm)
            elif ';=)' in g['name']:
                res_type[tesubtype]['pc-known'].append(tpm)

        elif g['coding'] == 'noncoding':
            res_type[tesubtype]['ncrna-all'].append(tpm)

            if ';=' in g['name']:
                res_type[tesubtype]['ncrna-known'].append(tpm)
            elif ';~' in g['name']:
                res_type[tesubtype]['ncrna-variant'].append(tpm)
            elif ';!' in g['name']:
                res_type[tesubtype]['ncrna-unknown'].append(tpm)

print(res_type.keys())
for te in res_type:
    if True in [i in te for i in ['SINE', 'LINE', 'LTR', 'Retroposon']]:
        print(te)
        for t in class_dict().keys(): # pc-all, ncrna-all' etc.
            data = {te: res_type[te][t], 'noTE': res[t]['nonTE']}
            gldraw.beanplot(filename='by_tesubtype/viol_{0}-{1}.png'.format(te, t), data=data,
                figsize=[2,1.8], beans=False, ylims=[-4, 8])
            gldraw.boxplot(filename='by_tesubtype/box_{0}-{1}.png'.format(te, t), data=data.values(), labels=data.keys(),
                figsize=[2,1.8], ylims=[-4, 8])


