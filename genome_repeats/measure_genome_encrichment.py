
'''

Measure the number of each class; type of TE in the genome, and output a frequency table;

'''

import os, sys, numpy, itertools, pickle, gzip
from collections import defaultdict
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from glbase3 import glload, genelist, expression, config, delayedlist

config.draw_mode = ['png', 'svg', 'pdf']

sys.path.append('../')
import shared

genome_repeats = delayedlist('hg38_rmsk.txt.gz', format={'force_tsv': True, 'name': 10, 'class': 11, 'family': 12},
    gzip=True)

print(genome_repeats)

res = {} # defaultdicts can be dangerous
res_by_full_name = {}

keeps = set(['LTR', 'LINE', 'SINE', 'DNA', 'Retroposon'])

for idx, gene in enumerate(genome_repeats):
    if gene['class'] not in keeps:
        continue

    te_type = '{g[class]}:{g[family]}'.format(g=gene)
    full_name = '{g[class]}:{g[family]}:{g[name]}'.format(g=gene).replace('-int', '')
    # the -int suffix is deprecated in the dfam

    if '?' in te_type:
        continue
    
    if te_type not in res:
        res[te_type] = 0
    res[te_type] += 1
    if full_name not in res_by_full_name:
        res_by_full_name[full_name] = 0
    res_by_full_name[full_name] += 1
    
    if (idx+1) % 1e5 == 0:
        print('{:,}'.format(idx+1))
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

# And by full name;
oh = open('te_freqs_by_full_name.txt', 'w')
oh.write('name\tfreq\n')
for k in sorted(res_by_full_name.keys()):
    oh.write('{0}\t{1}\n'.format(k, res_by_full_name[k]))
oh.close()

oh = open('te_freqs_by_full_name.pickle', 'wb')
pickle.dump(res_by_full_name, oh)
oh.close()

fig = plot.figure(figsize=[3,4])

labs = list(reversed(sorted(res.keys())))
vals = [res[k] for k in labs]
cols = [shared.get_col(k) for k in labs]

ax = fig.add_subplot(1,1,1)
fig.subplots_adjust(left=0.3, top=0.92)

ax.barh(numpy.arange(len(labs)), width=vals, height=0.8, color=cols)
ax.set_yticklabels(labs)
ax.set_yticks(numpy.arange(len(labs)))
[t.set_fontsize(6) for t in ax.get_yticklabels()]
[t.set_fontsize(6) for t in ax.get_xticklabels()]
ax.yaxis.label.set_size(6)

fig.savefig('freq.pdf')

# And by full name;
fig = plot.figure(figsize=[3,15])
fig.subplots_adjust(top=0.99, left=0.5, right=0.95, bottom=0.01)
labs = list(reversed(sorted(res_by_full_name.keys())))
vals = [res_by_full_name[k] for k in labs]
cols = [shared.get_col(k) for k in labs]

ax = fig.add_subplot(1,1,1)
ax.barh(numpy.arange(len(labs)), width=vals, height=0.8, color=cols)
ax.set_ylim([-1, len(labs)+1])
ax.set_yticklabels(labs)
ax.set_yticks(numpy.arange(len(labs)))
[t.set_fontsize(2) for t in ax.get_yticklabels()]
[t.set_fontsize(2) for t in ax.get_xticklabels()]

fig.savefig('freq_res_by_full_name.pdf')

