'''

For the masked peptides, see if we can find them in the MS data;

'''

import glob, sys, os, numpy
from glbase3 import *
import matplotlib.pyplot as plot

res = {}
'''
for filename in glob.glob('mass_spec_searches/per_transcript_hits*.tsv'):
    if 'prime' in filename:
        continue
    if 'class' in filename:
        continue
    stub = os.path.split(filename)[1].replace('.tsv', '').split('-')[1]
    pep_hits = genelist(filename, format={'name': 0, 'hits': 1, 'force_tsv': True})
    num_hits = sum([i>=1 for i in pep_hits['hits']])

    res[stub] = [num_hits, len(pep_hits)]
'''

all_matches = glload('results_gene.glb')

print(all_matches)
# I need to get tables for all of the peptides put into the search

res_per_gene = {}
for filename in glob.glob('2.blast_searches/masked/*.glb'):
    if 'prime' in filename:
        continue
    if 'class' in filename:
        continue
    stub = os.path.split(filename)[1].replace('.glb', '').split('-')[1]
    pep_hits = glload(filename).getColumns(['name'])

    res_per_gene[stub] = {k: 0 for k in pep_hits['name']}

for transcript in all_matches:
    full_name = '|'.join([transcript['name'], transcript['transcript_id'], transcript['enst']])
    res_per_gene[transcript['class']][full_name] += 1

res = {}
# final result histogram
for stub in res_per_gene:
    pep_hits = res_per_gene[stub].keys()
    num_hits = sum([i>=1 for i in res_per_gene[stub].values()])

    res[stub] = [num_hits, len(pep_hits)]
print(res)

fig = plot.figure(figsize=[3,1.0])
fig.subplots_adjust(left=0.5)
ax = fig.add_subplot(111)

ys = numpy.arange(len(res))

print(res.keys())
order = [
    'table_novel_coding',
    'table_insertion_alternate_cds',
    'table_new_STOP', 'table_new_ATG',
    'table_frameshift_insertion',

    ] # bottom to top

num_hits = numpy.array([res[k][0] for k in order])
len_hits = numpy.array([res[k][1] for k in order]) # starts at 0 so no need to subtract
print(num_hits, len_hits)
percs = (num_hits / len_hits) * 100.0
print(num_hits)
print(len_hits)
print(percs)
ax.barh(ys, len_hits)

ax.barh(ys, num_hits)
print(res.keys())
ax.set_xlabel('Number of proteins with >20 Amino acids')
ax.set_yticks(ys)
ax.set_yticklabels(order)
[t.set_fontsize(6) for t in ax.get_yticklabels()]
[t.set_fontsize(6) for t in ax.get_xticklabels()]
for y, p, x in zip(ys, percs, num_hits):
    ax.text(x+4, y, s='{0} ({1:.1f}%)'.format(x, p), va='center', fontsize=6)

fig.savefig('summary.png')
fig.savefig('summary.pdf')


