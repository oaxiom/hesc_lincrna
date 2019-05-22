

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location

draw = 'png'

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

locs = {
    'AluJb-LIN28B': location(loc='chr6:104,932,879-104,967,549'), # Rough location of hte AluJb-LIN28. A bit hard to identify inconclusivelyFrom the Wang Ting paper;
    'annotation': location(loc='chr1:778782-800871').expand(20000),
    'SVA': location(loc='chr8:86184842-86185042'),
    }

for k in locs:
    print(k, locs[k])
    c.draw.setLocation(loc=locs[k])
    scale = 1.0
    if draw == 'svg':
        scale = 0.3
    c.draw.exportImage("%s/%s_%s.%s" % (draw, k, str(locs[k]).replace(":", "-"), draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.
