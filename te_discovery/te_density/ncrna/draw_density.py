
'''

Measure the density of TE insertions, scaled across the length of the transcript.

'''

import glob, sys, os, gzip, numpy, math
import matplotlib.pyplot as plot
import matplotlib.cm as cm
from glbase3 import glload, utils, expression, genelist, genome_sql

sys.path.append('../../../')
import shared
sys.path.append('../')
import shared_draw

draw = 'png'

doms = glload('../../te_transcripts/transcript_table_HSC_SR_PB_merged.mapped.glb')
gencode = glload(os.path.expanduser('~/hg38/hg38_gencode_v29.glb')).getColumns(['enst', 'cds_loc'])
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

doms = doms.map(genelist=gencode, key='enst')

# preprocss the doms list to remove non-coding genes;
newdoms = []
for gene in doms:
    if ';NC;' in gene['name']:
        newdoms.append(gene)
doms = genelist()
doms.load_list(newdoms)

print('List now %s transcripts long' % len(doms))

shared_draw.draw_density('ncrna_all_tes.png', doms, None)

# split them up by TE family
type_subtype = {}
for TE in dfam:
    t_s = '%s:%s' % (TE['type'], TE['subtype'])
    if t_s not in type_subtype:
        type_subtype[t_s] = []
    type_subtype[t_s].append(TE['name'])

res = {}

for ts in type_subtype:
    res[ts] = shared_draw.draw_density('by_te_family/ncrna_%s.png' % (ts,), doms, set(type_subtype[ts]))

# turn it into a megatable:
shared_draw.draw_heatmap('te_heat.png', res)

res = {}
for ts in dfam:
    res[ts['name']] = shared_draw.draw_density('by_te_type/ncrna_%s_%s_%s.png' % (ts['type'], ts['subtype'], ts['name'],), doms, [ts['name'],])

shared_draw.draw_heatmap('te_heat_by_type.png', res)
