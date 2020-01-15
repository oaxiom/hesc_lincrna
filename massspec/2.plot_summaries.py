'''

For the masked peptides, see if we can find them in the MS data;

'''

import glob, sys, os, numpy
from glbase3 import *
import matplotlib.pyplot as plot

res = {}
for filename in glob.glob('mass_spec_searches/per_transcript_hits*.tsv'):
    if 'prime' in filename:
        continue
    if 'class' in filename:
        continue
    stub = os.path.split(filename)[1].replace('.tsv', '').split('-')[1]
    pep_hits = genelist(filename, format={'name': 0, 'hits': 1, 'force_tsv': True})
    num_hits = sum([i>=1 for i in pep_hits['hits']])

    res[stub] = [num_hits, len(pep_hits)]

fig = plot.figure(figsize=[2.5,1.4])
ax = fig.add_subplot(111)
ax
ys = numpy.arange(len(res))

print(res.keys())
order = ['table_new_STOP', 'table_new_ATG', 'table_frameshift_insertion',  ] # bottom to top

num_hits = numpy.array([res[k][0] for k in order])
len_hits = numpy.array([res[k][1] for k in order]) # starts at 0 so no need to subtract
percs = (num_hits / len_hits) * 100.0
print(num_hits)
print(len_hits)
print(percs)
ax.barh(ys, len_hits)

ax.barh(ys, num_hits)
print(res.keys())
ax.set_xlabel('Number of proteins with >50 Amino acids')
ax.set_yticks(ys)
ax.set_yticklabels(order)

for y, p, x in zip(ys, percs, num_hits):
    ax.text(x+4, y, s='{:.1f}%'.format(p), va='center')


fig.savefig('summary.png')
fig.savefig('summary.pdf')


