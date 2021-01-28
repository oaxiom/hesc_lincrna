

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload, genelist

draw = 'pdf'

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

#gllocs = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')

gllocs = glload('../../transcript_assembly/packed/all_genes.glb') # needed for noTE

locs = [
    'HPSCSR.276998.297',
    'HPSCSR.72619.2',
    'HPSCSR.21764.4',
    'HPSCSR.327312.6',
    'HPSCSR.304643.1',
    'HPSCLR.4844.3',
    'HPSCSR.162330.3',
    'HPSCSR.216914.100',
    'HPSCSR.165091.37', # HERVK
    'HPSCSR.165091.26', # noTE
    'HPSCSR.329.8',
    'HPSCSR.341.9',
    'HPSCLR.14380.7',
    'HPSCSR.643.6',
    'HPSCSR.88755.6',
    'HPSCSR.30639.18',
    'HPSCSR.126056.12',
    'HPSCSR.17088.41',
    'HPSCSR.163834.5',
    'HPSCSR.60772.41',
    ]
locs = genelist(loadable_list=[{'transcript_id': k} for k in locs])
ll = locs.map(genelist=gllocs, key='transcript_id')

ll.linearData.append({'loc': location(loc='chr20:51,781,374-51,860,697'), 'name': 'SALL4-MLT1J'})
ll.linearData.append({'loc': location(loc='chr6:104,918,050-105,083,558'), 'name': 'LIN28B-AluJb'})

ll._optimiseData()

print(ll)

for gene in ll:
    scale = 1.0
    if draw == 'pdf':
        scale = 0.3

    c.draw.setLocation(loc=gene['loc'].expand(len(gene['loc']) / 10))
    if 'enst' not in gene:
        print(gene['loc'])
        c.draw.exportImage("{0}/{1}.{2}".format(draw, gene['name'], draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason
    else:
        print(gene['name'])
        c.draw.exportImage("{0}/{1}_{2}_{3}.{4}".format(draw, gene['name'], gene['transcript_id'], gene['enst'], draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.


