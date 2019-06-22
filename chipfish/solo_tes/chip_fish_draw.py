

import sys, os, glob
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload

draw = 'svg'

[os.remove(f) for f in glob.glob('%s/*/*/*.%s' % (draw, draw))]

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

gllocs = glload('../../te_discovery/solo_tes/solo_tes.glb')

for gene in gllocs:
    destination = 'none'

    path = '%s/%s/' % (draw, ';'.join(sorted(gene['te_type'].split('; '))))
    if not os.access('%s' % (path), os.R_OK | os.W_OK):
        os.mkdir('%s' % (path))

    print(gene['name'])
    dd = len(gene['loc']) / 20.0
    c.draw.setLocation(loc=gene['loc'].expand(dd))
    scale = 1.0
    if draw == 'svg':
        scale = 0.3

    c.draw.exportImage("%s/%s_%s.%s" % (path, gene['name'], str(gene['loc']).replace(":", "-"), draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.