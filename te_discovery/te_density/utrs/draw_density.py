
'''

Measure the density of TE insertions, scaled across the length of the transcript.

'''

import glob, sys, os, gzip, numpy, math
import matplotlib.pyplot as plot
from glbase3 import glload, utils, expression, genelist, genome_sql

sys.path.append('../../../')
import shared
sys.path.append('../')
import shared_draw

draw = 'png'

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
gencode = glload('../../te_transcripts/transcript_table_gencode_pc.glb')
gencode_sliced = gencode.getColumns(['cds_loc', 'transcript_id', 'loc'])
gencode_sliced = gencode_sliced.renameKey('transcript_id', 'enst')

print(doms)

dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
doms = doms.map(genelist=gencode_sliced, key='enst') # You can't do this. CDS locations are not garunteed to be accurate

# preprocss the doms list to remove non-coding genes;
newdoms = []

type = {'known': [], 'novel': []}

for gene in doms:
    if gene['coding'] == 'noncoding':
        continue
    if ';~' in gene['name']: # CDS locations are not accurate in these; Can I get them from FEELnc?
        type['novel'].append(gene)
        continue
    if len(gene['cds_loc']) == 0:
        continue

    type['known'].append(gene)

gltypes = {'known': genelist(), 'novel': genelist()}
[gltypes[k].load_list(type[k]) for k in gltypes]
print(gltypes)

dataset_all = {
    'GENCODE': gencode_sliced,
    'ES-': doms.get(key='expression', value='depleted'),
    'ES:': doms.get(key='expression', value='unbiased'),
    'ES+': doms.get(key='expression', value='enriched'),
    }

for t in ['known']: # Novel is currently empty;
    dataset = {
        'GENCODE': dataset_all['GENCODE'].map(key='enst', genelist=gltypes[t]),
        'ES-': dataset_all['ES-'].map(key='enst', genelist=gltypes[t]),
        'ES:': dataset_all['ES:'].map(key='enst', genelist=gltypes[t]),
        'ES+': dataset_all['ES+'].map(key='enst', genelist=gltypes[t]),
        }

    shared_draw.draw_density_utrs('all_genes_all_tes-{0}.png'.format(t), dataset, None)

    # preprocss the doms list to remove non-coding genes;

    # split them up by gene-type;

    # split them up by TE family
    type_subtype = {}
    for TE in dfam:
        t_s = '%s:%s' % (TE['type'], TE['subtype'])
        if t_s not in type_subtype:
            type_subtype[t_s] = []
        type_subtype[t_s].append(TE['name'])

    for ts in type_subtype:
        shared_draw.draw_density_utrs('by_te_family/all_genes_{0}-{1}.png'.format(t, ts,), dataset, set(type_subtype[ts]))

