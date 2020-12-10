
'''

Build the final annotation tables;

'''


import glob, sys, os, gzip
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6

sys.path.append('../../../')
import shared
from glbase3 import glload, utils, expression, genelist

all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
te = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
not_te = te.map(genelist=all_genes, key='transcript_id', logic='notright')

#print(all_genes)
#print()
#print(te)

all_genes_ensg = set([g['ensg'] for g in all_genes if 'ENSG' in g['ensg']])
has_at_least_one_coding_isoform = set([g['ensg'] for g in all_genes if ';C;' in g['name'] and 'ENSG' in g['ensg']])

# coding and non-coding;
shared.pie('coding_genes_with_a_noncoding_transcript.png', [len(all_genes_ensg)-len(has_at_least_one_coding_isoform), len(has_at_least_one_coding_isoform)], ['no non-coding', 'has non-coding transcript'], 'Coding genes with a noncoding transcript')

#print(len(te), len(not_te), len(all_genes), len(te)+len(not_te))

# Piechart, percent of transcripts containing a TE:
shared.pie('pies/te_all.png', [len(not_te), len(te)], ['no-TE', 'TE'], 'All')

# collect stats:
res = {
    'ES+': {'TE': 0, 'nonTE': 0},
    'ES:': {'TE': 0, 'nonTE': 0},
    'ES-': {'TE': 0, 'nonTE': 0},
    'pc-ES+': {'TE': 0, 'nonTE': 0},
    'pc-ES:': {'TE': 0, 'nonTE': 0},
    'pc-ES-': {'TE': 0, 'nonTE': 0},
    'ncrna-ES+': {'TE': 0, 'nonTE': 0},
    'ncrna-ES:': {'TE': 0, 'nonTE': 0},
    'ncrna-ES-': {'TE': 0, 'nonTE': 0},
    }

data = {'TE': te, 'nonTE': not_te}

for k in data:
    for g in data[k]:
        if g['expression'] == 'enriched':
            res['ES+'][k] += 1
        elif g['expression'] == 'unbiased':
            res['ES:'][k] += 1
        elif g['expression'] == 'depleted':
            res['ES-'][k] += 1

        if g['coding'] == 'coding':
            if g['expression'] == 'enriched':
                res['pc-ES+'][k] += 1
            elif g['expression'] == 'unbiased':
                res['pc-ES:'][k] += 1
            elif g['expression'] == 'depleted':
                res['pc-ES-'][k] += 1

        if g['coding'] == 'noncoding':
            if g['expression'] == 'enriched':
                res['ncrna-ES+'][k] += 1
            elif g['expression'] == 'unbiased':
                res['ncrna-ES:'][k] += 1
            elif g['expression'] == 'depleted':
                res['ncrna-ES-'][k] += 1

title_map = {
    'ES+': 'hPSC-enriched',
    'ES:': 'hPSC-neutral',
    'ES-': 'hPSC-depleted',
    'pc-ES+': 'hPSC-enriched (PC)',
    'pc-ES:': 'hPSC-neutral (PC)',
    'pc-ES-': 'hPSC-depleted (PC)',
    'ncrna-ES+': 'hPSC-enriched (ncRNA)',
    'ncrna-ES:': 'hPSC-neutral (ncRNA)',
    'ncrna-ES-': 'hPSC-depleted (ncRNA)',
    }

for k in res:
    shared.pie('pies/te_%s.png' % k, [res[k]['nonTE'], res[k]['TE']], ['no-TE', 'TE'], title_map[k])

shared.split_bar('bar.png'.format(k), res, cols=['#ff7f0e', '#1f77b4'])

# pickle the results
import pickle
shared.pickle_it('results.pickle', res)
