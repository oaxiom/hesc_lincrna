

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload, genelist

draw = 'pdf'

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

gllocs = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')

locs = [
    'HPSCSR.276998.297',
    'HPSCSR.72619.2',
    'HPSCSR.21764.4',
    'HPSCSR.327312.6',
    'HPSCSR.304643.1',
    'HPSCLR.4844.3',
    'HPSCSR.162330.3',
    'HPSCSR.216914.100',
    #'HPSCSR.
    #'HPSCSR.

    ]
locs = genelist(loadable_list=[{'transcript_id': k} for k in locs])

ll = locs.map(genelist=gllocs, key='transcript_id')

print(ll)

for gene in ll:
    print(gene['name'])
    c.draw.setLocation(loc=gene['loc'].expand(len(gene['loc']) / 10))
    scale = 1.0
    if draw == 'pdf':
        scale = 0.3
    c.draw.exportImage("{0}/{1}_{2}_{3}.{4}".format(draw, gene['name'], gene['transcript_id'], gene['enst'], draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.
