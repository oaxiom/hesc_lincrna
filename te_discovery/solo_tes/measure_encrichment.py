
'''

SOLO TEs are these transcripts that contain intact, or semi intact unspliced transcripts.

As we don't trust the short read data to assemble these, we only consider them from the pacbio data:

'''

import sys, numpy, itertools
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from glbase3 import glload, genelist, expression, config

config.draw_mode = ['png', 'svg']

sys.path.append('../../')
import shared

sys.path.append('../simple_summaries/')
import pies

all_te_transcripts = glload('../te_transcripts/transcript_table_merged.mapped.glb')
dfam = genelist('../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

solotes = glload('solo_tes.glb')

res = {}

for gene in solotes:
    for te_type in gene['te_type'].split('; '):
        if '.' in te_type:
            continue
        if te_type not in res:
            res[te_type] = 0
        res[te_type] += 1

# work out the pair-wise frequency:
print(res.keys())
paired = {k2: {k:0 for k in res.keys()} for k2 in res.keys()}
for gene in solotes:
    tes = [i for i in gene['te_type'].split('; ') if '.' not in i]
    te_types = list(set(tes))

    pairs = list(set(itertools.product(tes, tes)))
    print(pairs)
    for p in pairs:
        if p[0] == p[1]:
            continue
        paired[p[0]][p[1]] += 1

print(paired)

newe = [{'name': k, 'conditions': [paired[k][k2] for k2 in paired.keys()]} for k in paired.keys()]
e = expression()
e.load_list(newe, cond_names=list(paired.keys()))
e = e.filter_low_expressed(5,1)
e = e.sliceConditions(e['name']) # Make it square
e.heatmap('paired.png', heat_wid=0.3, heat_hei=0.2, bracket=[0, 50], grid=True,
    cmap=cm.plasma,
    optimal_ordering=True)


fig = plot.figure(figsize=[3,4])

labs = list(reversed(sorted(res.keys())))
vals = [res[k] for k in labs]
cols = [shared.get_col(k) for k in labs]

ax = fig.add_subplot(1,1,1)
fig.subplots_adjust(left=0.3, top=0.92)

print(vals)

ax.barh(numpy.arange(len(labs)), width=vals, height=0.8, color=cols)
ax.set_yticklabels(labs)
ax.set_yticks(numpy.arange(len(labs)))
[t.set_fontsize(6) for t in ax.get_yticklabels()]
[t.set_fontsize(6) for t in ax.get_xticklabels()]
ax.yaxis.label.set_size(6)

fig.savefig('freq.png')
fig.savefig('freq.svg')
