'''

For Fig S3

'''

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location

draw = 'pdf'

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs_2.txt"))

locs = {
    'BRCA1': location(loc='chr17:43,033,480-43,144,384'),
    'ATAD3B': location(loc='chr1:1,465,244-1,504,369'),
    'HPSCLR.4844.3': location(loc='chr4:87,906,610-87,949,605').expand(15000),
    'TP53I3': location(loc='chr2:24,075,487-24,087,161'),
    }

for k in locs:
    print(k, locs[k])
    c.draw.setLocation(loc=locs[k])
    scale = 1.0
    if draw == 'pdf':
        scale = 0.3
    c.draw.exportImage("{}/{}.{}".format(draw, k, draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.
