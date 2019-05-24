
import numpy, sys, os
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared

config.draw_mode = ['png', 'svg']

expn = glload('../data/te_per1Mbp_seq.glb')

expn = expn.sliceConditions(['gencode.ncrna', 'GRCh38.ncrna', 'custom.ncrna'])

expn = expn.norm_multi_fc({'gencode.ncrna':
    ['GRCh38.ncrna', 'custom.ncrna']}, pad=0.1)

print(expn.getConditionNames())

types = []

FC = 1.0

# get the sig ones;
newe = []
for e in expn:
    if e['conditions'][2] > FC:
        merged_type = '%s:%s' % (e['type'], e['subtype'])
        e['merged_type'] = merged_type
        newe.append(e)
        types.append(merged_type)

expn.load_list(newe, cond_names=expn.getConditionNames())

print('FC>%s = %s' % (FC, len(expn)))
print('\n'.join(sorted(set(types))))

row_cols = shared.get_cols(expn['merged_type'])

res = expn.heatmap('fc_gt1.svg', bracket=[-2, 2], figsize=[6,12], row_font_size=6,
    row_colbar=row_cols, heat_wid=0.02*len(expn.getConditionNames()), grid=True)

print('\n'.join(reversed(res['reordered_rows'])))

# Get a heatmap per-type:
all_type = set(expn['type'])
for t in all_type:
    newe = expn.get(key='type', value=t)
    newe.heatmap('type_%s.png' % t, bracket=[-2, 2], figsize=[6,12], row_font_size=6,
        row_colbar=shared.get_cols(newe['merged_type']), heat_wid=0.02*len(expn.getConditionNames()), col_cluster=False,
        heat_hei=0.007*len(newe), grid=True)

# Split out the LTR subtypes
ltrs_only = expn.get(key='type', value='LTR')
all_type = set(ltrs_only['subtype'])
for t in all_type:
    newe = ltrs_only.get(key='subtype', value=t)
    newe.heatmap('type_LTR_%s.png' % t, bracket=[-2, 2], figsize=[6,12], row_font_size=6,
        row_colbar=shared.get_cols(newe['merged_type']), heat_wid=0.02*len(expn.getConditionNames()), col_cluster=False,
        heat_hei=0.007*len(newe), grid=True)

# Plot these for all samples:

# Count the subtypes and get a barchart:
res = {}
for e in expn:
    if e['merged_type'] not in res:
        res[e['merged_type']] = 0
    res[e['merged_type']] += 1

labs = list(reversed(sorted(res.keys())))
vals = [res[k] for k in labs]
cols = [shared.col_keys[k] for k in labs]

fig = plot.figure(figsize=[1.9,3.0])
ax = fig.add_subplot(1,1,1)
fig.subplots_adjust(left=0.5, top=0.92)

print(vals)

ax.barh(numpy.arange(len(labs)), width=vals, height=0.8, color=cols)
ax.set_yticklabels(labs)
ax.set_yticks(numpy.arange(len(labs)))
[t.set_fontsize(6) for t in ax.get_yticklabels()]
[t.set_fontsize(6) for t in ax.get_xticklabels()]
ax.yaxis.label.set_size(6)

fig.savefig('freq.png')
fig.savefig('freq.svg')

gl_freqs = genelist()
gl_freqs.load_list([{'name': k, 'freq': res[k]} for k in res])
gl_freqs.save('gl_freqs.glb')
gl_freqs.saveTSV('gl_freqs.tsv')


