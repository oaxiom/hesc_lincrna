
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

all_genes = glload('../../../gencode/hg38_gencode_v32.glb')
te = glload('../../te_transcripts/transcript_table_gencode_all.glb')
not_te = te.map(genelist=all_genes, key='transcript_id', logic='notright')

pc = glload('../../../gencode/hg38_gencode_v32.pc.glb')
ncrna = glload('../../../gencode/hg38_gencode_v32.ncrna.glb')

all_genes_ensg = set([g['ensg'] for g in all_genes if 'ENSG' in g['ensg']])
pc_genes_ensg = set(pc['ensg'])
has_at_least_one_coding_isoform = set([g['ensg'] for g in all_genes if g['ensg'] in pc_genes_ensg])

# coding and non-coding;
shared.pie('coding_genes_with_a_noncoding_transcript.png', [len(all_genes_ensg)-len(has_at_least_one_coding_isoform), len(has_at_least_one_coding_isoform)], ['no non-coding', 'has non-coding transcript'], 'Coding genes with a noncoding transcript')

print(len(te), len(not_te), len(all_genes), len(te)+len(not_te))

# Piechart, percent of transcripts containing a TE:
shared.pie('pies/te_all.png', [len(not_te), len(te)], ['no-TE', 'TE'], 'All')

# collect stats:
res = {
    'pc': {'TE': 0, 'nonTE': 0},
    'ncrna': {'TE': 0, 'nonTE': 0},
    'non_coding_version_of_coding_gene': {'TE': 0, 'nonTE': 0},
    }

data = {'TE': te, 'nonTE': not_te}

pc_genes_enst = set(pc['enst'])
ncrna_genes_enst = set(ncrna['enst'])

for k in data:
    for g in data[k]:
        if g['enst'] in pc_genes_enst:
            res['pc'][k] += 1
        elif g['enst'] in ncrna_genes_enst:
            res['ncrna'][k] += 1

        # non_coding_version_of_coding_gene
        if g['enst'] in ncrna_genes_enst:
            if g['ensg'] not in has_at_least_one_coding_isoform:
                continue
            res['non_coding_version_of_coding_gene'][k] += 1


title_map = {'pc': 'protein-coding',
    'ncrna': 'non-coding RNA',
    'all_novel': 'Novel transcripts',
    'pc_novel': 'protein-coding (novel)',
    'ncrna_novel': 'non-coding RNA (novel)',
    'non_coding_version_of_coding_gene': 'non-coding version of coding gene'}

for k in res:
    shared.pie('pies/te_%s.png' % k, [res[k]['nonTE'], res[k]['TE']], ['no-TE', 'TE'], title_map[k])

shared.split_bar('bar.png'.format(k), res, key_order=['TE', 'nonTE'], cols=['#ff7f0e', '#1f77b4']) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

# pickle the results
import pickle
shared.pickle_it('results.pickle', res)
