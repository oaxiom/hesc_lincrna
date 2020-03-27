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
all_te_transcripts = glload('../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
all_data = all_matches.map(genelist=all_te_transcripts, key='transcript_id')

all_fastas = genelist('2.blast_searches/all_masked_peptides.fa', format=format.fasta) # Just for getting the positions;
all_fastas = {i['name']: i['seq'] for i in all_fastas}
dfam = genelist('../te_discovery/dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

res_per_gene = {}
for solo_te in all_data:
    res_per_gene[solo_te['transcript_id']] = {'te_derived_peptide': [],
        'derived_peptide': []}

# inside TE generated in script 1
res_per_TE_type = {}
for m in all_data:
    if m['insideTE'] == 'No':
        res_per_gene[m['transcript_id']]['derived_peptide'].append(m['peptide_string'])
        continue
    res_per_gene[m['transcript_id']]['te_derived_peptide'].append(m['peptide_string'])

    if m['insideTE'] not in res_per_TE_type:
        res_per_TE_type[m['insideTE']] = set([])
    res_per_TE_type[m['insideTE']].add(m['peptide_string'])

res_unq_pep = {}
for k in res_per_gene:
    # Collapse to number of unique peptides
    res_per_gene[k]['num_unq_derived_peptide'] = len(set(res_per_gene[k]['derived_peptide']))
    res_per_gene[k]['num_unq_te_derived_peptide'] = len(set(res_per_gene[k]['te_derived_peptide']))
    #print(res_per_gene[k]['num_unq_derived_peptide'], res_per_gene[k]['num_unq_te_derived_peptide'])


res = {
    'Inside TE': {'2 or more unique peptides': 0, '1 Unique peptide': 0, 'No peptides': 0},
    'Outside TE': {'2 or more unique peptides': 0, '1 Unique peptide': 0, 'No peptides': 0},
    }

for k in res_per_gene:
    if res_per_gene[k]['num_unq_derived_peptide'] >= 2:
        res['Outside TE']['2 or more unique peptides'] += 1
    elif res_per_gene[k]['num_unq_derived_peptide'] == 1:
        res['Outside TE']['1 Unique peptide'] += 1
    else:
        res['Outside TE']['No peptides'] += 1

    if res_per_gene[k]['num_unq_te_derived_peptide'] >= 2:
        res['Inside TE']['2 or more unique peptides'] += 1
    elif res_per_gene[k]['num_unq_te_derived_peptide'] == 1:
        res['Inside TE']['1 Unique peptide'] += 1
    else:
        res['Inside TE']['No peptides'] += 1

# Plot 1:
shared.bar('peptides.png', res,
    key_order=['2 or more unique peptides', '1 Unique peptide'],
    cols=['#d62728', '#ff7f0e', '#2ca02c',  ],
    figsize=[3,2]) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])

# Save the unqTE table:
gl = genelist()
newl = []
for k in res_per_TE_type:
    for pep in res_per_TE_type[k]:
        newl.append({'name': k, 'num_peps': len(res_per_TE_type[k]), 'peptide': pep})
#newl = [{'name': k, 'num_peps': len(res_per_TE_type[k]), 'peptides':res_per_TE_type[k]} for k in res_per_TE_type]
gl.load_list(newl)
gl.sort('name')
gl.saveTSV('unique_TE-derived_peptides.tsv', key_order=['name', 'num_peps'])

# Plot 2: Peptide-derived TEs:
fig = plot.figure(figsize=[3,5])
fig.subplots_adjust(left=0.6, right=0.90, top=0.998, bottom=0.03)#, hspace=0.2, wspace=0.2)
#ax = plot.subplot2grid((2,2), (0,0), rowspan=2)
ax = fig.add_subplot(111)
sorted_names = sorted(res_per_TE_type.keys())
print(sorted_names)

vals = {k: len(res_per_TE_type[k]) for k in sorted_names}

ys = numpy.arange(len(res_per_TE_type))
labs = vals.keys()
vals = vals.values()

ax.barh(ys, vals)
ax.set_yticks(ys)
ax.set_yticklabels(labs)
xlim = ax.get_xlim()[1]
xlim += (xlim/10)
for y, c in zip(ys, vals):
    ax.text(xlim, y, c, fontsize=6, va='center')
ax.tick_params(right=True)
[t.set_fontsize(6) for t in ax.get_yticklabels()]
[t.set_fontsize(6) for t in ax.get_xticklabels()]
ax.set_ylim([-1, len(ys)])
fig.savefig('detected_TEs.png')
fig.savefig('detected_TEs.pdf')

# Plot 3: As plot 2, but top 10
fig = plot.figure(figsize=[3,2])
fig.subplots_adjust(left=0.6, right=0.90, top=0.998, bottom=0.03)#, hspace=0.2, wspace=0.2)
#ax = plot.subplot2grid((2,2), (0,0), rowspan=2)
ax = fig.add_subplot(111)
res_per_TE_type = {k: len(res_per_TE_type[k]) for k in res_per_TE_type}
res_per_TE_type = {k: res_per_TE_type[k] for k in sorted(res_per_TE_type, key=res_per_TE_type.get, reverse=False)}

ys = numpy.arange(20)
vals = list(res_per_TE_type.values())[-20:]
labs = list(res_per_TE_type.keys())[-20:]

ax.barh(ys, vals)

xlim = ax.get_xlim()[1]
xlim += (xlim/10)
for y, c in zip(ys, vals):
    ax.text(xlim, y, c, fontsize=6, va='center')
ax.tick_params(right=True)
ax.set_yticks(ys)
ax.set_yticklabels(labs)
[t.set_fontsize(6) for t in ax.get_yticklabels()]
[t.set_fontsize(6) for t in ax.get_xticklabels()]
ax.set_ylim([-1, len(ys)])
fig.savefig('detected_TEs_top20.png')
fig.savefig('detected_TEs_top20.pdf')

