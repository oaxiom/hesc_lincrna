import sys, os, itertools
from collections import Counter
import numpy as np
import matplotlib.pyplot as plot
import matplotlib.tri as tri
from glbase3 import *

import shared_conservation

# collect three things:
# 1. The total PhyloP score of the transcript
# 2. score for TE-containing bits
# 3. score for non-TE containign bits;

dfam = genelist('../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
dfam_dict = {}
for te in dfam:
    dfam_dict[te['name']] = '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name'])

transcripts = glload('../te_transcripts/transcript_table_merged.mapped.glb')
gl = glload('phyloP_conservation_table.glb')
print(gl)

t = 0.25

not_counted = 0
both_conserved = []
te_conserved = []
lncrna_conserved = []

for item in gl:
    if item['phyloP_tes'] > t and item['phyloP_nottes'] > t:
        both_conserved.append(item)
    elif item['phyloP_tes'] > t:
        te_conserved.append(item)
    elif item['phyloP_nottes'] > t:
        lncrna_conserved.append(item)
    else:
        not_counted += 1

print('Not counted     : {0:,}'.format(not_counted))
print('Both conserved  : {0:,}'.format(len(both_conserved)))
print('TE conserved    : {0:,}'.format(len(te_conserved)))
print('lncRNA conserved: {0:,}'.format(len(lncrna_conserved)))
print('Total TE-containing transcripts: {0:,}'.format(len(transcripts)))

gl = genelist()
gl.load_list(both_conserved)
both_conserved = gl
gl = genelist()
gl.load_list(te_conserved)
te_conserved = gl
gl = genelist()
gl.load_list(lncrna_conserved)
lncrna_conserved = gl

all_data = {'both_conserved': both_conserved.map(genelist=transcripts, key='transcript_id'),
    'te_conserved': te_conserved.map(genelist=transcripts, key='transcript_id'),
    'lncrna_conserved': lncrna_conserved.map(genelist=transcripts, key='transcript_id')
    }

for k in all_data:
    # convert to a list of doms:
    doms = []
    for t in all_data[k]:
        doms += [dfam_dict[d['dom']] for d in t['doms']]

    c = Counter(doms)
    c = c.most_common(20)#.items()
    print(c)

    vals = [i[1] for i in c]
    labs = [i[0] for i in c]
    vals.reverse()
    labs.reverse()

    fig = plot.figure(figsize=[2,2])
    fig.subplots_adjust(left=0.5, top=0.97)
    ax = fig.add_subplot(111)

    ys = np.arange(len(vals))

    ax.barh(ys, vals)

    ax.set_xlabel('Number of TE domains')
    ax.set_yticks(ys)
    ax.set_yticklabels(labs)

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    #for y, p, x in zip(ys, percs, num_hits):
    #    ax.text(x+4, y, s='{0} ({1:.1f}%)'.format(x, p), va='center', fontsize=6)

    fig.savefig('class_summary-{0}.pdf'.format(k))

