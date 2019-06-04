
import numpy, sys, os
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../')
import shared

config.draw_mode = 'png'

# load in all the tables:
gl_pc = glload('fc_heats_enriched_pc_only/gl_freqs.glb')
gl_ncrna = glload('fc_heats_enriched_ncrna_only/gl_freqs.glb')
gl_all = glload('fc_heats_enriched_all_only/gl_freqs.glb')

# merge all the flies into one genelist:
all_tes = gl_pc['name'] + gl_ncrna['name'] + gl_all['name']
all_tes = set([i for i in all_tes if '.' not in i])

# Count the subtypes and get a barchart:
res = {k: [0,0,0] for k in all_tes}
for n, gl in enumerate([gl_all, gl_pc, gl_ncrna]):
    for e in gl:
        if e['name'] not in res:
            continue
        res[e['name']][n] = e['freq']

title_dict = {0: 'All', 1: 'PC', 2: 'ncRNA'}

fig = plot.figure(figsize=[4,2.4])
for c in [0,1,2]:
    labs = list(reversed(sorted(res.keys())))
    vals = [res[k][c] for k in labs]
    cols = [shared.get_col(k) for k in labs]

    ax = fig.add_subplot(1,3,c+1)
    fig.subplots_adjust(left=0.2, top=0.92)

    print(vals)

    ax.barh(numpy.arange(len(labs)), width=vals, height=0.8, color=cols)
    if c == 0:
        ax.set_yticklabels(labs)
    else:
        ax.set_yticklabels('')
    ax.set_title(title_dict[c])
    ax.set_yticks(numpy.arange(len(labs)))
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    ax.yaxis.label.set_size(6)

fig.savefig('freq.png')
fig.savefig('freq.svg')
fig.savefig('freq.pdf')

print('\n'.join(reversed(labs)))

# save all:
'''
gl_freqs = genelist()
gl_freqs.load_list([{'name': k, 'freq_all': res[k]} for k in res])
gl_freqs.save('gl_freqs.glb')
gl_freqs.saveTSV('gl_freqs.tsv')
'''
