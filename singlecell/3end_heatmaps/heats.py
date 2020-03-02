
import sys, os
import glob, os, sys, numpy
from operator import itemgetter
from glbase3 import *
import matplotlib.cm as cm
config.draw_mode = 'pdf'

trks = {
    'scRNA-seq (+ strand)': flat_track(filename='../flats/scrnaseq_strp.flat', name='scRNA-seq (+ strand)'),#, norm_factor= 12.22)
    'scRNA-seq (- strand)': flat_track(filename='../flats/scrnaseq_strm.flat', name='scRNA-seq (- strand)'),#, norm_factor= 11.91)
    }

# cid1 = sites affected by low K reprogramming
genes = glload('../custom_3end_gtf/all_3ends.glb')

genes = {
    'scRNA-seq (+ strand)': genes.get(key='strand', value='+').removeDuplicatesByLoc('pointify_expand', 'loc', 500), # remove dupes, to make it look more clear
    'scRNA-seq (- strand)': genes.get(key='strand', value='-').removeDuplicatesByLoc('pointify_expand', 'loc', 500)
    }

big_tab = None
pad_array = None
bracket = [2, 7]

distance = 1000

for k in trks:
    res = trks[k].heatmap(
        filename='{0}.png'.format(k),
        distance=distance,
        genelist=genes[k],
        cmap=cm.inferno,
        norm_by_read_count=True,
        sort_by_intensity=True,
        respect_strand=False,
        size=[3,9],
        log=2,
        log_pad=0.1,
        bracket=[0, 5]
        )

