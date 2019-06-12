
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
gencode = glload('../../../gencode/hg38_gencode_v30.glb')
gencode_sliced = gencode.getColumns(['enst', 'cds_loc'])

#print(gencode)
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

print(doms)
doms = doms.map(genelist=gencode_sliced, key='enst')

# preprocss the doms list to remove non-coding genes;
newdoms = []
novel = []
for gene in doms:
    if ';NC;' in gene['name']:
        continue
    if ';~' in gene['name']: # CDS locations are not accurate in these; Can I get them from FEELnc?
        novel.append(gene)
    if len(gene['cds_loc']) == 0:
        continue
    newdoms.append(gene)
known_doms = genelist()
known_doms.load_list(newdoms)
novel_doms = genelist()
novel_doms.load_list(novel)
print(doms)
print(novel)
print('List now %s known transcripts long' % len(doms))
print('List now %s novel transcripts long' % len(novel_doms))
shared_draw.draw_density_utrs('all_genes_all_tes.png', known_doms, novel_doms, gencode, None)

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
    shared_draw.draw_density_utrs('by_te_family/all_genes_%s.png' % (ts,), known_doms, novel_doms, gencode, set(type_subtype[ts]))

