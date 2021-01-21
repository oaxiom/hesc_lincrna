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

all_data = all_matches.map(genelist=all_te_transcripts, key='transcript_id')
all_data = all_data.removeDuplicates(['transcript_id', 'peptide_string'])

res = {'Inside TE': {'m': 0}, 'Outside TE': {'m': 0}}

for pep in all_data:
    print(pep['insideTE'])
    if pep['insideTE'] == 'No':
        res['Outside TE']['m'] += 1
    else:
        res['Inside TE']['m'] += 1

# Plot 1:
shared.bar('peptides-peps.png', res,
    key_order=['m'],
    cols=['#d62728', '#ff7f0e'],
    figsize=[3,2]) #, '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'])
