
'''

Measure the expression of TEs by family and then by class, plot a violin and a heatmap;

'''

import glob, sys, os, gzip, numpy, math, scipy.stats
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6

from statsmodels.stats.multitest import multipletests
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

dfam_dict = {}
for te in dfam:
    dfam_dict[te['name']] = '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name'])

def class_dict():
    return {
        'pc-all': {'TE': [], 'nonTE': []},
        'ncrna-all': {'TE': [], 'nonTE': []},
        }

def dict_builder():
    return {k: [] for k in class_dict()}

data = {'TE': contains_te, 'nonTE': contains_not_te}
res = class_dict()

# Fill the data tables
for datatype in data:
    for g in data[datatype]:
        tpm = math.log2(g['TPM']+0.1)

        if g['coding'] == 'coding':
            res['pc-all'][datatype].append(tpm)
        elif g['coding'] == 'noncoding':
            res['ncrna-all'][datatype].append(tpm)

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
        elif g['coding'] == 'noncoding':
            res_type[full_name]['ncrna-all'].append(tpm)

p_scatter = {'pc-all': [], 'ncrna-all': []}

for te in sorted(res_type):
    for t in class_dict().keys(): # pc-all, ncrna-all' etc.
        if True in [typ in te for typ in ['SINE', 'LINE', 'LTR', 'Retroposon']]:
            data = {te: res_type[te][t], 'noTE': res[t]['nonTE']}
            #print(data)

            if len(data[te]) <= 20:
                continue

            p = scipy.stats.mannwhitneyu(data[te], data['noTE'], alternative='two-sided')[1]
            M = utils.fold_change(2**numpy.mean(data['noTE']), 2**numpy.mean(data[te]), pad=0.01)# fold-change
            A = numpy.mean(data[te]) # average expression;

            p_scatter[t].append({'name': te, 'p': p, 'M': M, 'A': A, 'n': len(data[te])})

# adjust the p values pls
#print(p_scatter)

for t in class_dict().keys():
    gl = genelist()
    gl.load_list(p_scatter[t])

    # adjust the P:
    ps = gl['p']
    _, p_adj, _, _ = multipletests(ps, method='fdr_bh')

    for p, item in zip(p_adj, gl):
        item['q'] = p
        item['-log10q'] = -math.log10(p)
    gl._optimiseData()
    gl.saveTSV('tab_{0}.tsv'.format(t), key_order=['name', 'p', 'q'])

    spot_cols = []
    for q, M in zip(gl['-log10q'], gl['M']):
        if q > 1.301: # q=0.05
            if M > 0.58: #~1.5 fold
                spot_cols.append('red')
            elif M <-0.58:
                spot_cols.append('blue')
            else:
                spot_cols.append('grey')
        else:
            spot_cols.append('grey')

    config.draw_mode = 'svg'
    shared.nice_scatter(y=gl['-log10q'], x=gl['M'], figsize=[2,2], spot_size=12,
        spot_cols=spot_cols,
        filename='MA-{0}.png'.format(t), label=gl['name'], hlines=[1.301], vlines=[0])
    config.draw_mode = 'pdf'

    if os.path.exists(t):
        [os.remove(f) for f in glob.glob('{0}/*.pdf'.format(t))]
    else:
        os.mkdir(t)

    for te in gl:
        #print(t, te, res_type[te['name']][t])
        if te['q'] < 0.05:
            data = {te['name']: res_type[te['name']][t], 'noTE': res[t]['nonTE']}
            gldraw.boxplot(filename='{0}/box_{1}-{2}.png'.format(t, te['name'], t), data=data.values(), labels=data.keys(),
                figsize=[2,1.8], ylims=[-4, 8], title='q={0:.2e} n={1}'.format(te['q'], te['n']))
        if True in [i in te['name'] for i in ('HERVH', 'HERVFH21', 'LTR7', 'AluSp', 'AluY', 'HERVK')]:
            data = {te['name']: res_type[te['name']][t], 'noTE': res[t]['nonTE']}
            gldraw.boxplot(filename='{0}-box_{1}-{2}.png'.format(t, te['name'], t), data=data.values(), labels=data.keys(),
                figsize=[2,1.8], ylims=[-4, 8], title='q={0:.2e} n={1}'.format(te['q'], te['n']))
