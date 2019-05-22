

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload, genelist

draw = 'png'

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

annot = glload(os.path.expanduser('~/hg38/hg38_ensembl_v95_enst.glb'))
gllocs = glload('../../te_discovery/te_transcripts/transcript_table_HSC_SR_PB_merged.mapped.glb')
#print(gllocs)

locs = ['MFSD13A'
    ]
locs = genelist(loadable_list=[{'name': k} for k in locs])

ll = locs.map(genelist=annot, key='name')
ll = ll.map(genelist=gllocs, key='enst')

print(ll)

for gene in ll:
    print(gene['name'])
    c.draw.setLocation(loc=gene['loc'].expand(len(gene['loc']) / 10))
    scale = 1.0
    if draw == 'svg':
        scale = 0.3
    c.draw.exportImage("%s/%s_%s.%s" % (draw, gene['name'], gene['enst'], draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.
