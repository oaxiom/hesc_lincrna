
'''

Measure the expression of TEs by family and then by class, plot a violin and a heatmap;

'''

import glob, sys, os, gzip, numpy, math, scipy.stats, pickle
import matplotlib
import scipy.stats
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6

#from statsmodels.stats.multitest import multipletests
from collections import defaultdict

import matplotlib.pyplot as plot
import matplotlib.cm as cm
from glbase3 import glload, utils, expression, genelist, genome_sql, draw, config
config.draw_mode = ['pdf']

sys.path.append('../../../')
import shared

draw_type = 'pdf'

all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
contains_te = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
contains_not_te = contains_te.map(genelist=all_genes, key='transcript_id', logic='notright')
contains_not_te.saveTSV('not_TE.tsv')
dfam_dict = {}
for te in dfam:
    dfam_dict[te['name']] = '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name'])

def class_dict():
    return {
        'pc-all': {'TE': [], 'no TE': []},
        'ncrna-all': {'TE': [], 'no TE': []},
        }

def dict_builder():
    return {k: [] for k in class_dict()}

data = {'no TE': contains_not_te, 'TE': contains_te}
res = class_dict()

# Fill the data tables
for datatype in data:
    for g in data[datatype]:
        tpm = math.log2(g['TPM'])
        if g['coding'] == 'coding':
            res['pc-all'][datatype].append(tpm)
        elif g['coding'] == 'noncoding':
            res['ncrna-all'][datatype].append(tpm)

# q
q = {'no TE': 1.0,
    'TE': scipy.stats.ttest_ind(res['pc-all']['no TE'], res['pc-all']['TE'], equal_var=False)[1]
    }#, alternative='two-sided')[1]}
shared.boxplots('pc-all.pdf', res['pc-all'], q)
q = {'no TE': 1.0,
    'TE': scipy.stats.ttest_ind(res['ncrna-all']['no TE'], res['ncrna-all']['TE'], equal_var=False)[1]}
shared.boxplots('ncrna-all.pdf', res['ncrna-all'], q)

q = {'Coding no TE': 1.0,
    'Coding TE': scipy.stats.ttest_ind(res['pc-all']['no TE'], res['pc-all']['TE'], equal_var=False)[1],
    'Non-coding no TE': 1.0,
    'Non-coding TE': scipy.stats.ttest_ind(res['ncrna-all']['no TE'], res['ncrna-all']['TE'], equal_var=False)[1]
    }#, alternative='two-sided')[1]}

'''
res = {
    'Coding no TE': res['pc-all']['no TE'],
    'Coding TE': res['pc-all']['TE'],
    'Non-coding no TE': res['ncrna-all']['no TE'],
    'Non-coding TE': res['ncrna-all']['TE']
    }
shared.boxplots('all.pdf', res, q)
'''
