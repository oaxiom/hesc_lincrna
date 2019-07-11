

import sys, os, glob
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload

sys.path.append('../../')
import shared

draw = 'png'

[os.remove(f) for f in glob.glob('%s/*/*/*.%s' % (draw, draw))]

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

gllocs = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')

for k in gllocs:
    destination = 'none'

    destination, alpha = shared.classify_transcript(k['name'])
    if not os.access('%s/%s' % (draw, destination), os.R_OK | os.W_OK):
        os.mkdir('%s/%s' % (draw, destination))

    print(k['name'])
    dd = len(k['loc']) / 20.0
    c.draw.setLocation(loc=k['loc'].expand(dd))
    scale = 1.0
    if draw == 'svg':
        scale = 0.3

    path = '%s/%s/%s' % (draw, destination, alpha)

    if not os.access(path, os.R_OK | os.W_OK):
        os.mkdir(path)
    c.draw.exportImage("%s/%s_%s.%s" % (path, k['name'], str(k['loc']).replace(":", "-"), draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.
