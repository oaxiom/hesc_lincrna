
'''

Measure the density of TE insertions, scaled across the length of the transcript.

'''

import glob, sys, os, gzip, numpy, math
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 6

import matplotlib.pyplot as plot
import matplotlib.cm as cm
from glbase3 import glload, utils, expression, genelist, genome_sql


sys.path.append('../../../')
import shared
sys.path.append('../')
import shared_draw

draw = 'png'
type = 'pc'

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
gencode = glload('../../../gencode/hg38_gencode_v32.{}.glb'.format(type)).getColumns(['enst', 'cds_loc'])
gencode_doms = glload('../../te_transcripts/transcript_table_gencode_%s.glb' % type)
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

# preprocss the doms list to remove non-coding genes;
newdoms = []
for gene in doms:
    if ';C;' in gene['name']:
        newdoms.append(gene)
doms = genelist()
doms.load_list(newdoms)

print('List now %s transcripts long' % len(doms))

dataset = {
    'GENCODE': gencode_doms,
    'ES-': doms.get(key='expression', value='depleted'),
    'ES:': doms.get(key='expression', value='unbiased'),
    'ES+': doms.get(key='expression', value='enriched'),
    }

shared_draw.draw_density('%s_all_tes.png' % (type,), dataset, None)

# split them up by TE family
type_subtype = {}
for TE in dfam:
    t_s = '%s:%s' % (TE['type'], TE['subtype'])
    if t_s not in type_subtype:
        type_subtype[t_s] = []
    type_subtype[t_s].append(TE['name'])

res = {}
for ts in type_subtype:
    r = shared_draw.draw_density('by_te_family/%s_%s.png' % (type, ts,), dataset, set(type_subtype[ts]))
    if r:
        res[ts] = r

# split them up by TE subtype
type_subtype = {}
for TE in dfam:
    t_s = '%s:%s' % (TE['type'], TE['subtype'])
    if t_s not in type_subtype:
        type_subtype[t_s] = []
    type_subtype[t_s].append(TE['name'])

res = {}
for ts in type_subtype:
    r = shared_draw.draw_density('by_te_family/%s_%s.png' % (type, ts,), dataset, set(type_subtype[ts]))
    if r:
        res[ts] = r

# turn it into a megatable:
shared_draw.draw_heatmap('te_heat_by_type_%s.png' % (type,), res, dataset)

res = {}
for ts in dfam:
    r = shared_draw.draw_density('by_te_type/%s_%s_%s_%s.png' % (type, ts['type'], ts['subtype'], ts['name'],), dataset, [ts['name'],])
    if r:
        res[ts['name']] = r

shared_draw.draw_heatmap('te_heat_by_subtype_%s.png' % (type, ), res, dataset)
