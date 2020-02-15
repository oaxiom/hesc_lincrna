
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
    'pc': {'TE': 0, 'nonTE': 0},
    'ncrna': {'TE': 0, 'nonTE': 0},
    'non_coding_version_of_coding_gene': {'TE': 0, 'nonTE': 0},
    'novel': {'TE': 0, 'nonTE': 0},
    'novel_pc': {'TE': 0, 'nonTE': 0},
    'novel_ncrna': {'TE': 0, 'nonTE': 0},
    'all_variant': {'TE': 0, 'nonTE': 0},
    'pc_variant': {'TE': 0, 'nonTE': 0},
    'ncrna_variant': {'TE': 0, 'nonTE': 0},
    }

data = {'TE': te, 'nonTE': not_te}

for k in data:
    for g in data[k]:
        if ';!' in g['name']:
            res['novel'][k] += 1
        if ';~' in g['name']:
            res['all_variant'][k] += 1

        if ';C;' in g['name']:
            res['pc'][k] += 1
            if ';~' in g['name']:
                res['pc_variant'][k] += 1
            if ';!' in g['name']:
                res['novel_pc'][k] += 1

        elif ';NC;' in g['name']:
            res['ncrna'][k] += 1
            if ';~' in g['name']:
                res['ncrna_variant'][k] += 1
            if ';!' in g['name']:
                res['novel_ncrna'][k] += 1

        # non_coding_version_of_coding_gene
        if ';NC;' in g['name']:
            if 'ENST' in g['enst']:
                # filter out LINCs, etc:
                if g['ensg'] not in has_at_least_one_coding_isoform:
                    continue
                res['non_coding_version_of_coding_gene'][k] += 1


title_map = {'pc': 'protein-coding',
    'ncrna': 'non-coding RNA',
    'all_variant': 'Variant transcripts',
    'pc_variant': 'protein-coding (variant)',
    'ncrna_variant': 'non-coding RNA (variant)',
    'non_coding_version_of_coding_gene': 'non-coding version of coding gene',
    'novel': 'Novel',
    'novel_pc': 'Novel protein-coding',
    'novel_ncrna': 'Novel non-coding',
    }

for k in res:
    shared.pie('pies/te_%s.png' % k, [res[k]['nonTE'], res[k]['TE']], ['no-TE', 'TE'], title_map[k])

shared.split_bar('bar.png'.format(k), res, cols=['#ff7f0e', '#1f77b4'])

# pickle the results
import pickle
shared.pickle_it('results.pickle', res)
