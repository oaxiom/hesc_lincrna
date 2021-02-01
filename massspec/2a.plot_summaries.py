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
    if 'table_variant_coding_but_noTE' in filename:
        continue

    stub = os.path.split(filename)[1].replace('.glb', '').split('-')[1]
    pep_hits = glload(filename).getColumns(['name'])

    res_per_gene[stub] = {k: 0 for k in pep_hits['name']}

for transcript in all_matches:
    if 'table_variant_coding_but_noTE' in transcript['class']:
        continue

    full_name = '|'.join([transcript['name'], transcript['transcript_id'], transcript['enst']])
    res_per_gene[transcript['class']][full_name] += 1

res = {}
# final result histogram
for stub in res_per_gene:
    pep_hits = res_per_gene[stub].keys()
    num_hits = sum([i>=1 for i in res_per_gene[stub].values()])
    num_hits2pep = sum([i>=2 for i in res_per_gene[stub].values()])

    res[stub] = [num_hits, num_hits2pep, len(pep_hits)] #these overlap, so don't adjust the bottoms.
print(res)

fig = plot.figure(figsize=[4.0, 2.6])
fig.subplots_adjust(left=0.4)
ax = fig.add_subplot(111)

ys = numpy.arange(len(res))

print(res.keys())
order = [
    'table_frameshift_insertion',
    'table_new_STOP',
    'table_new_ATG',
    'table_insertion_alternate_cds',
    'table_noncoding_to_coding_withTE',
    'table_noncoding_to_coding_noTE',
    'table_novel_coding',
    ] # top to bottom
order.reverse()


num_hits1 = numpy.array([res[k][0] for k in order])
num_hits2 = numpy.array([res[k][1] for k in order])
len_hits = numpy.array([res[k][2] for k in order]) # starts at 0 so no need to subtract
print(num_hits1, num_hits2, len_hits)
percs1 = (num_hits1 / len_hits) * 100.0
percs2 = (num_hits2 / len_hits) * 100.0

ax.barh(ys, len_hits, color='lightgrey')
ax.barh(ys, num_hits1, color='tab:orange')
ax.barh(ys, num_hits2, color='tab:red')

print(res.keys())
ax.set_xlabel('Number of proteins with >20 Amino acids', fontsize=6)
ax.set_yticks(ys)
ax.set_yticklabels(order)
[t.set_fontsize(6) for t in ax.get_yticklabels()]
[t.set_fontsize(6) for t in ax.get_xticklabels()]

for y, p, x in zip(ys, percs1, num_hits1):
    ax.text(x+4, y-0.2, s='1 peptide ({0}; {1:.1f}%)'.format(x, p), va='center', fontsize=6)
for y, p, x in zip(ys, percs2, num_hits2):
    ax.text(x+4, y+0.2, s='2+ peptides ({0}; {1:.1f}%)'.format(x, p), va='center', fontsize=6)

fig.savefig('summary.pdf')


