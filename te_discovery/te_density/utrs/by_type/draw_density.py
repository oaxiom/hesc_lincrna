
'''

Measure the density of TE insertions, scaled across the length of the transcript.

'''

import glob, sys, os, gzip, numpy, math
import matplotlib.pyplot as plot
from glbase3 import glload, utils, expression, genelist, genome_sql

sys.path.append('../../../../')
import shared
sys.path.append('../../')
import shared_draw

draw = 'png'

doms = glload('../../../te_transcripts/transcript_table_merged.mapped.glb')
dfam = genelist('../../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
cds = glload('../../../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
cds = {i['transcript_id']: i for i in cds}
newl = []
#for g in doms:
#    g['enst'] = g['enst'].split('.')[0]
#    newl.append(g)
#doms.load_list(newl)
gencode = glload('../../../../transcript_assembly/get_CDS/gencode_cds.glb').map(genelist=doms, key='enst')

# preprocss the doms list to remove non-coding genes;
newdoms = []

type = {'all': [], 'known': [], 'novel': []}

for gene in doms:
    if gene['coding'] == 'noncoding':
        continue

    if '!' in gene['tags']: # I am not considering these, as they look dubious;
        continue

    if gene['transcript_id'] not in cds:
        print('Warning {0} not found'.format(gene['transcript_id']))
        continue

    gene['cds_local_locs'] = cds[gene['transcript_id']]['cds_local_locs']
    gene['tlength'] = cds[gene['transcript_id']]['tlength']

    if gene['cds_local_locs'][0] == gene['cds_local_locs'][1]: # I am unsure about the cds_loc;
        continue

    if '=' in gene['tags']:
        type['known'].append(gene)
    if '~' in gene['tags']:
        type['novel'].append(gene)

    type['all'].append(gene)

gltypes = {'known': genelist(), 'novel': genelist(), 'all': genelist()}
[gltypes[k].load_list(type[k]) for k in gltypes]

dataset_all = {
    'GENCODE': gencode,
    'ES-': doms.get(key='expression', value='depleted'),
    'ES:': doms.get(key='expression', value='unbiased'),
    'ES+': doms.get(key='expression', value='enriched'),
    }

for t in ['known', 'novel', 'all']:
    if os.path.exists(t):
        [os.remove(f) for f in glob.glob('{0}/*.p*'.format(t))]
    else:
        os.mkdir(t)

    dataset = {
        'GENCODE': dataset_all['GENCODE'], #.map(key='enst', genelist=gltypes[t]),
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
        t_s = '%s:%s:%s' % (TE['type'], TE['subtype'], TE['name'])

        shared_draw.draw_density_utrs('{0}/{0}-{1}.png'.format(t, t_s,), dataset, [TE['name']])

