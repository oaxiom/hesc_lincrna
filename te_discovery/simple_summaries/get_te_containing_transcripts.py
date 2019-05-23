
'''

Build the final annotation tables;

'''


import glob, sys, os, gzip
import matplotlib.pyplot as plot
from glbase3 import glload, utils, expression, genelist

def lab_func(pct, allvals):
    absolute = int(pct/100.*sum(allvals))
    return "{:.1f}%\n{:,}".format(pct, absolute)

def pie(filename, data, labels, title=''):
    fig = plot.figure(figsize=[2,2])
    ax = fig.add_subplot(111)

    ax.pie(data, labels=labels, autopct=lambda pct: lab_func(pct, data))

    ax.set_title(title)
    fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.svg'))

all_genes = glload('../../transcript_assembly/packed/all_genes.glb')

te = glload('../te_transcripts/transcript_table_HSC_SR_PB_merged.mapped.glb')
not_te = te.map(genelist=all_genes, key='transcript_id', logic='notright')

has_at_least_one_coding_isoform = set([g['ensg'] for g in all_genes if ';C;' in g['name']])

print(len(te), len(not_te), len(all_genes), len(te)+len(not_te))

print(te)

# Piechart, percent of transcripts containing a TE:
pie('all_te_containing.png', [len(not_te), len(te)], ['no-TE', 'TE'], 'All')

# collect stats:
res = {
    'pc': {'TE': 0, 'nonTE': 0},
    'ncrna': {'TE': 0, 'nonTE': 0},
    'non_coding_version_of_coding_gene': {'TE': 0, 'nonTE': 0},
    }

non_coding_gene_names = ['U1', 'U7', '.', 'LINC', '-AS', 'MIR']

for g in te:
    if ';C;' in g['name']:
        res['pc']['TE'] += 1
    elif ';NC;' in g['name']:
        res['ncrna']['TE'] += 1

    # non_coding_version_of_coding_gene
    if ';NC;' in g['name']:
        if 'ENST' in g['enst']:
            # filter out LINCs, etc:
            if g['ensg'] not in has_at_least_one_coding_isoform:
                continue
            res['non_coding_version_of_coding_gene']['TE'] += 1

for g in not_te:
    if ';C;' in g['name']:
        res['pc']['nonTE'] += 1
    elif ';NC;' in g['name']:
        res['ncrna']['nonTE'] += 1

    # non_coding_version_of_coding_gene
    if ';NC;' in g['name']:
        if 'ENST' in g['enst']:
            # filter out LINCs, etc:
            if g['ensg'] not in has_at_least_one_coding_isoform:
                continue
            res['non_coding_version_of_coding_gene']['nonTE'] += 1

title_map = {'pc': 'protein-coding',
    'ncrna': 'non-coding RNA',
    'non_coding_version_of_coding_gene': 'non-coding version of coding gene'}

for k in res:
    pie('%s_te_containing.png' % k, [res[k]['nonTE'], res[k]['TE']], ['no-TE', 'TE'], title_map[k])

