
'''

SOLO TEs are these transcripts that contain intact, or semi intact unspliced transcripts.

As we don't trust the short read data to assemble these, we only consider them from the pacbio data:

'''

import sys, numpy, itertools, pickle
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from glbase3 import glload, genelist, expression, config, utils

config.draw_mode = ['png', 'pdf']

sys.path.append('../../')
import shared

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
paired = {k2: {k:0 for k in res.keys()} for k2 in res.keys()}
for gene in solotes:
    tes = [i for i in gene['te_type'].split('; ') if '.' not in i]
    te_types = list(set(tes))

    pairs = list(set(itertools.product(tes, tes)))
    for p in pairs:
        if p[0] == p[1]:
            continue
        paired[p[0]][p[1]] += 1

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
fig.savefig('freq.pdf')

# Measure enrichment versus the genome:

oh = open('../../genome_repeats/te_freqs.pickle', 'rb')
genome_repeat_freqs = pickle.load(oh)
oh.close()

total_tes = sum(genome_repeat_freqs.values())

total_solotes = sum(res.values())

fc = {}
print(res.keys())
for k in sorted(res.keys()):
    if k not in genome_repeat_freqs: # Some of the names are a little different, it's some really rare ones though so just ignore
        print("Warning: couldn't find {0} in genome, number of TEs in transcripts={1}".format(k, res[k]))
        continue
    #print(k, res[k], genome_repeat_freqs[k])

    observed = res[k]
    expected = total_solotes * (genome_repeat_freqs[k] / total_tes)
    res[k] = (observed)/(expected)

    fc[k] = utils.fold_change(expected, observed)

    print('obs={0:.2f}\texp={1:.2f}\tenrich={2:.2f}\t{3}'.format(observed, expected, res[k], k))

    #res[k] = res[k] / genome_repeat_freqs[k] * 1e6

e = expression(loadable_list=[{'name': k, 'conditions': [fc[k]]} for k in fc], cond_names=['enrichment'])
e.sort_sum_expression()
r = e.heatmap(filename='solo_te_enrichment.png', heat_wid=0.03, grid=True, row_cluster=False, heat_hei=0.015*len(e), bracket=[-2,2],
    draw_numbers=True, draw_numbers_threshold=1, draw_numbers_fmt='*')
print('\n'.join(reversed(r['reordered_rows'])))
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
ax.set_xlim([0, 3])

fig.savefig('freq_enrichment.png')
fig.savefig('freq_enrichment.svg')
fig.savefig('freq_enrichment.pdf')

# Key classes to look at:
key_classes = [
    'Retroposon:SVA',
    'DNA:hAT-Tag1',
    'DNA:PiggyBac',
    'LTR:ERVK',
    'LTR:ERV1',
    'LINE:L1',
    'SINE:Alu',
    'DNA:MULE-MuDR']

# Split these down a level for specific classes
for cls in key_classes:
    res = {}
    for gene in solotes:
        if cls in gene['te_type']: # Has one of the enriched classes;
            for TE in gene['te_fullanmes'].split('; '):
                te_family = ':'.join(TE.split(':')[0:2])
                #print("'{0}'".format(te_family))

                if cls in te_family:
                    print(te_family)
                    if TE not in res:
                        res[TE] = 0
                    res[TE] += 1
    fc = {}
    for TE in res:
        observed = res[k]
        expected = total_solotes * (genome_repeat_freqs[k] / total_tes)
        res[k] = (observed)/(expected)
        fc[k] = utils.fold_change(expected, observed)
        print('obs={0:.2f}\texp={1:.2f}\tenrich={2:.2f}\t{3}'.format(observed, expected, res[k], k))

        #res[k] = res[k] / genome_repeat_freqs[k] * 1e6

    e = expression(loadable_list=[{'name': k, 'conditions': [fc[k]]} for k in fc], cond_names=['enrichment'])
    e.sort_sum_expression()
    r = e.heatmap(filename='solo_te_enrichment_by_class.png', heat_wid=0.03, grid=True, row_cluster=False, heat_hei=0.015*len(e), bracket=[-2,2],
        draw_numbers=True, draw_numbers_threshold=1, draw_numbers_fmt='*')
    print(res)
