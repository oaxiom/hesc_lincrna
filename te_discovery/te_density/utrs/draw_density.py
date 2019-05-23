
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

doms = glload('../../te_transcripts/transcript_table_HSC_SR_PB_merged.mapped.glb')
gencode = glload(os.path.expanduser('~/hg38/hg38_gencode_v29.glb')).getColumns(['enst', 'cds_loc'])
#print(gencode)
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

print(doms)
doms = doms.map(genelist=gencode, key='enst')

# preprocss the doms list to remove non-coding genes;
newdoms = []
for gene in doms:
    if ';NC;' in gene['name']:
        continue
    if ';~' in gene['name']: # CDS locations are not accurate in these;
        continue
    if len(gene['cds_loc']) == 0:
        continue
    newdoms.append(gene)
doms = genelist()
doms.load_list(newdoms)
print(doms)
print('List now %s transcripts long' % len(doms))

shared_draw.draw_density_utrs('all_genes_all_tes.png', doms, None)

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
    shared_draw.draw_density_utrs('by_te_family/all_genes_%s.png' % (ts,), doms, set(type_subtype[ts]))

