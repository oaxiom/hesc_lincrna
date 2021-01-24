'''

For the masked peptides, see if we can find them in the MS data;

'''

import glob, sys, os, numpy
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../')
import shared

res = {}

all_matches = glload('results_gene.glb')
all_te_transcripts = glload('../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
cds = glload('../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
all_data = all_matches.map(genelist=all_te_transcripts, key='transcript_id').map(genelist=cds, key='transcript_id')

res = {
    'TE bps': {'n': 0},
    'Not TE bps': {'n': 0},
    }

for m in all_data:
    print(m)


    cds_left = m['cds_local_locs'][0]
    cds_right = m['cds_local_locs'][1]

    total_bps = cds_right - cds_left
    num_te_bps = 0

    for d in m['doms']:
        if d['span'][1] >= cds_left and d['span'][0] <= cds_right:
            # work out the num bps:
            max_left = max([d['span'][0], cds_left])
            max_rite = min([d['span'][1], cds_right])
            num_te_bps += max_rite - max_left

    res['TE bps']['n'] += num_te_bps
    res['Not TE bps']['n'] += total_bps - num_te_bps

shared.bar('percent_bps.png', res,
    key_order=['n'],
    cols=['#d62728'],
    figsize=[3,2])
