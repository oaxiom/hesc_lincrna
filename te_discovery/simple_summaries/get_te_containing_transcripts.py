
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
    fig = plot.figure(figsize=[1.3,1.3])
    ax = fig.add_subplot(111)

    wedges, texts, autotexts = ax.pie(data, labels=labels, autopct=lambda pct: lab_func(pct, data))
    plot.setp(autotexts, size=6)
    plot.setp(texts, size=6)

    ax.set_title(title, size=6)
    fig.savefig(filename)
    fig.savefig(filename.replace('.png', '.svg'))

all_genes = glload('../../transcript_assembly/packed/all_genes.glb')

te = glload('../te_transcripts/transcript_table_HSC_SR_PB_merged.mapped.glb')
not_te = te.map(genelist=all_genes, key='transcript_id', logic='notright')

all_genes_ensg = set([g['ensg'] for g in all_genes if 'ENSG' in g['ensg']])
has_at_least_one_coding_isoform = set([g['ensg'] for g in all_genes if ';C;' in g['name'] and 'ENSG' in g['ensg']])
print(has_at_least_one_coding_isoform)

# coding and non-coding;
pie('coding_genes_with_a_noncoding_transcript.png', [len(all_genes_ensg)-len(has_at_least_one_coding_isoform), len(has_at_least_one_coding_isoform)], ['no non-coding', 'has non-coding transcript'], 'Coding genes with a noncoding transcript')

print(len(te), len(not_te), len(all_genes), len(te)+len(not_te))

print(te)

# Piechart, percent of transcripts containing a TE:
pie('pies/te_all.png', [len(not_te), len(te)], ['no-TE', 'TE'], 'All')

# collect stats:
res = {
    'pc': {'TE': 0, 'nonTE': 0},
    'ncrna': {'TE': 0, 'nonTE': 0},
    'non_coding_version_of_coding_gene': {'TE': 0, 'nonTE': 0},
    'all_novel': {'TE': 0, 'nonTE': 0},
    'pc_novel': {'TE': 0, 'nonTE': 0},
    'ncrna_novel': {'TE': 0, 'nonTE': 0},
    }

non_coding_gene_names = ['U1', 'U7', '.', 'LINC', '-AS', 'MIR']

data = {'TE': te, 'nonTE': not_te}

for k in data:
    for g in data[k]:
        if ';~' in g['name']:
            res['all_novel'][k] += 1

        if ';C;' in g['name']:
            res['pc'][k] += 1
            if ';~' in g['name']:
                res['pc_novel'][k] += 1
        elif ';NC;' in g['name']:
            res['ncrna'][k] += 1
            if ';~' in g['name']:
                res['ncrna_novel'][k] += 1

        # non_coding_version_of_coding_gene
        if ';NC;' in g['name']:
            if 'ENST' in g['enst']:
                # filter out LINCs, etc:
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
    pie('pies/te_%s.png' % k, [res[k]['nonTE'], res[k]['TE']], ['no-TE', 'TE'], title_map[k])

