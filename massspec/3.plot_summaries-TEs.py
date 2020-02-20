'''

For the masked peptides, see if we can find them in the MS data;

'''

import glob, sys, os, numpy
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../')
import shared

res = {}

all_matches = glload('results_gene.glb')
all_solo_tes = glload('../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
all_data = all_matches.map(genelist=all_solo_tes, key='transcript_id')
print(all_data)

all_fastas = genelist('2.blast_searches/all_masked_peptides.fa', format=format.fasta) # Just for getting the positions;
all_fastas = {i['name']: i['seq'] for i in all_fastas}
dfam = genelist('../te_discovery/dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
print(all_fastas)

res_per_gene = {}
for solo_te in all_data:
    res_per_gene[solo_te['transcript_id']] = {'te_derived_peptide': [],
        'derived_peptide': []}

delete_chars = set('0123456789.+')
# I need to go through, work out where the peptide is positioned, and then see if it is coming out of a TE:
res_per_TE_type = {}
for m in all_data:
    # Sanitise the peptide for searching
    peptide = ''.join([i for i in m['peptide'] if i not in delete_chars])
    res_per_gene[m['transcript_id']]['derived_peptide'].append(peptide)

    # Get the position in the Peptide_fasta
    name = '|'.join([m['class'], m['name'].replace(' ', ''), m['transcript_id'], m['enst']])
    aa_seq = all_fastas[name]
    left = aa_seq.index(peptide) * 3
    rite = left+len(peptide) * 3

    # see if it's in a TE domain:
    for d in m['doms']:
        if d['span'][1] >= left and d['span'][0] <= rite:
            res_per_gene[m['transcript_id']]['te_derived_peptide'].append(peptide)

            te = dfam.get(key='name', value=d['dom'])[0]
            fullname = '{0}:{1}:{2}'.format(te['type'], te['subtype'], d['dom'])

            if fullname not in res_per_TE_type:
                res_per_TE_type[fullname] = 0
            res_per_TE_type[fullname] += 1

res_unq_pep = {}
for k in res_per_gene:
    # Collapse to number of unique peptides
    res_per_gene[k]['num_unq_derived_peptide'] = len(set(res_per_gene[k]['derived_peptide']))
    res_per_gene[k]['num_unq_te_derived_peptide'] = len(set(res_per_gene[k]['te_derived_peptide']))
    #print(res_per_gene[k]['num_unq_derived_peptide'], res_per_gene[k]['num_unq_te_derived_peptide'])


res = {
    'TE-derived': {'2 or more unique peptides': 0, '1 Unique peptide': 0, 'No peptides': 0},
    'detected': {'2 or more unique peptides': 0, '1 Unique peptide': 0, 'No peptides': 0},
    }
for k in res_per_gene:
    if res_per_gene[k]['num_unq_derived_peptide'] >= 2:
        res['detected']['2 or more unique peptides'] += 1
    elif res_per_gene[k]['num_unq_derived_peptide'] == 1:
        res['detected']['1 Unique peptide'] += 1
    else:
        res['detected']['No peptides'] += 1

    if res_per_gene[k]['num_unq_te_derived_peptide'] >= 2:
        res['TE-derived']['2 or more unique peptides'] += 1
    elif res_per_gene[k]['num_unq_derived_peptide'] == 1:
        res['TE-derived']['1 Unique peptide'] += 1
    else:
        res['TE-derived']['No peptides'] += 1

# Plot 1:
shared.split_bar('peptides.png', res,
    key_order=['2 or more unique peptides', '1 Unique peptide', 'No peptides'],
    cols=['#d62728', '#ff7f0e', '#2ca02c',  ]) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

# Plot 2: Peptide-derived TEs:
fig = plot.figure(figsize=[3,5.8])
fig.subplots_adjust(left=0.6, right=0.98, top=0.998, bottom=0.03)#, hspace=0.2, wspace=0.2)
#ax = plot.subplot2grid((2,2), (0,0), rowspan=2)

ax = fig.add_subplot(111)
res_per_TE_type = {k: res_per_TE_type[k] for k in sorted(res_per_TE_type, key=res_per_TE_type.get, reverse=False)}
print(res_per_TE_type)
ys = numpy.arange(len(res_per_TE_type))
ax.barh(ys, list(res_per_TE_type.values()))
ax.set_yticks(ys)
ax.set_yticklabels(res_per_TE_type.keys())
[t.set_fontsize(6) for t in ax.get_yticklabels()]
[t.set_fontsize(6) for t in ax.get_xticklabels()]
ax.set_ylim([-1, len(ys)])
fig.savefig('detected_TEs.png')
fig.savefig('detected_TEs.pdf')


