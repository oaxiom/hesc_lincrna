
import sys, os, glob
import matplotlib.cm as cm
from glbase3 import *
config.draw_mode = "png"
sys.path.append('../../')
import shared

expn = glload("../kallisto/kall_tpm-unmerged.glb")

print(expn)

[os.remove(f) for f in glob.glob('%s/*/*/*.%s' % (draw, draw))]
[os.remove(f) for f in glob.glob('%s/*/*.%s' % (draw, draw))]

for n, gene in enumerate(expn):
    destination, alpha = shared.classify_transcript(gene['name'])
    if not os.access('%s/%s' % (config.draw_mode, destination), os.R_OK | os.W_OK):
        os.mkdir('%s/%s' % (config.draw_mode, destination))

    path = '%s/%s/%s' % (config.draw_mode, destination, alpha)

    if not os.access(path, os.R_OK | os.W_OK):
        os.mkdir(path)

    expn.barh_single_item(value=gene['transcript_id'], key="transcript_id",
        filename='%s/%s.%s.%s.%s' % (path, gene['name'], gene['transcript_id'], gene['enst'], config.draw_mode),
        size=(4, 7), vert_space=0.8, #tree=tree,
        yticklabel_fontsize=6, xticklabel_fontsize=6)
