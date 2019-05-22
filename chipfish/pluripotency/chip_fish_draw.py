

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload, genelist

draw = 'png'

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

annot = glload(os.path.expanduser('~/hg38/hg38_ensembl_v95_enst.glb'))
gllocs = glload('../../te_discovery/transcripts/transcripts_with_te.glb')
#print(gllocs)

locs = ['SOX2', 'NANOG', 'SALL4', 'LIN28A', 'LIN28B', 'SALL1', 'POU5F1A',
    'DPPA2', 'DPPA3', 'DPPA5A', 'PRDM14', 'JARID2', 'SALL2', 'SALL3', 'TCF3',
    'HNRNPU',
    # MA Gang's possibles:
    'HNRNPK', 'DDX1', 'DDX50', 'BRCA2', 'BRCA1', 'TOP1', 'RAP3', 'TRIM25',
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
