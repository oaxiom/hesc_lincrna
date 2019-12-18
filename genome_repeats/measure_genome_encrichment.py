
'''

Measure the number of each class; type of TE in the genome, and output a frequency table;

'''

import os, sys, numpy, itertools, pickle
from collections import defaultdict
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from glbase3 import glload, genelist, expression, config, delayedlist

config.draw_mode = ['png', 'svg', 'pdf']

sys.path.append('../')
import shared

genome_repeats = delayedlist(os.path.expanduser('~/hg38/repeats/hg38_rmsk.tsv'), format={'force_tsv': True, 'name': 10, 'class': 11, 'family': 12})

print(genome_repeats)

res = defaultdict(int)

keeps = set(['LTR', 'LINE', 'SINE', 'DNA', 'Retroposon'])

for idx, gene in enumerate(genome_repeats):
    if gene['class'] not in keeps:
        continue

    #if gene['class'] == 'LTR':
    #    te_type = '{g[class]}:{g[family]}:{g[name]}'.format(g=gene)
    #else:
    te_type = '{g[class]}:{g[family]}'.format(g=gene)

    if '?' in te_type:
        continue

    res[te_type] += 1
    if (idx+1) % 1e5 == 0:
        print(idx+1)
        #break


# output the table;

oh = open('te_freqs.txt', 'w')
oh.write('name\tfreq\n')
for k in sorted(res.keys()):
    oh.write('{0}\t{1}\n'.format(k, res[k]))
oh.close()

oh = open('te_freqs.pickle', 'wb')
pickle.dump(res, oh)
oh.close()

print(res)

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
fig.savefig('freq.pdf')

