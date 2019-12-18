

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload, genelist

draw = 'png'

sys.path.append('../../')
import shared

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

#annot = glload(os.path.expanduser('~/hg38/hg38_ensembl_v95_enst.glb'))
#annot = annot.renameKey('name', 'gene_symbol')
gllocs = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')

for gene in gllocs:
    if ';ES+;' not in gene['name']:
        continue

    destination, alpha = shared.classify_transcript(gene['name'])
    if not os.access('%s/%s' % (draw, destination), os.R_OK | os.W_OK):
        os.mkdir('%s/%s' % (draw, destination))
    path = '%s/%s/%s' % (draw, destination, alpha)
    if not os.access(path, os.R_OK | os.W_OK):
        os.mkdir(path)


    print(gene['name'])
    c.draw.setLocation(loc=gene['loc'].expand(len(gene['loc']) / 10))
    scale = 1.0
    if draw == 'svg':
        scale = 0.3

    c.draw.exportImage("%s/%s_%s.%s" % (path, gene['name'], str(gene['loc']).replace(":", "-"), draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.
