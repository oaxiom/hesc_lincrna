
import sys, os
import glob, os, sys, numpy
from operator import itemgetter
from glbase3 import *
import matplotlib.cm as cm
config.draw_mode = 'pdf'

trks = {
    'scRNA-seq (+ strand)': flat_track(filename='../flats/scrnaseq_strp.flat', name='scRNA-seq (+ strand)'),
    'scRNA-seq (- strand)': flat_track(filename='../flats/scrnaseq_strm.flat', name='scRNA-seq (- strand)'),
    'scRNA-seq (+ strand) (genes on - strand)': flat_track(filename='../flats/scrnaseq_strp.flat', name='scRNA-seq (+ strand)'),
    'scRNA-seq (- strand) (genes on + strand)': flat_track(filename='../flats/scrnaseq_strm.flat', name='scRNA-seq (- strand)'),
    }

# cid1 = sites affected by low K reprogramming
genes = glload('../custom_3end_gtf/all_3ends.glb')

genes = {
    'scRNA-seq (+ strand)': genes.get(key='strand', value='+').removeDuplicatesByLoc('pointify_expand', 'loc', 500, use_strand=True), # remove dupes, to make it look more clear
    'scRNA-seq (- strand)': genes.get(key='strand', value='-').removeDuplicatesByLoc('pointify_expand', 'loc', 500, use_strand=True),

    'scRNA-seq (+ strand) (genes on - strand)': genes.get(key='strand', value='-').removeDuplicatesByLoc('pointify_expand', 'loc', 500, use_strand=True), # remove dupes, to make it look more clear
    'scRNA-seq (- strand) (genes on + strand)': genes.get(key='strand', value='+').removeDuplicatesByLoc('pointify_expand', 'loc', 500, use_strand=True)
    }

big_tab = None
pad_array = None

distance = 2000

for k in trks:
    res = trks[k].heatmap(
        filename='heat-{0}.png'.format(k),
        distance=distance,
        genelist=genes[k],
        cmap=cm.magma,
        norm_by_read_count=True,
        sort_by_intensity=True,
        respect_strand=False, # Don't need to respect strand, otherwise it gest confusing as both heatmaps appear the same
        size=[6,21],
        log=2,
        log_pad=0.01,
        bracket=[-2, 5]
        )

