import sys, os, itertools
import numpy as np
import matplotlib.tri as tri
from glbase3 import *

import shared_conservation

# collect three things:
# 1. The total PhyloP score of the transcript
# 2. score for TE-containing bits
# 3. score for non-TE containign bits;

gl = glload('phyloP_conservation_table.glb')
print(gl)

t = 0.25

not_counted = 0
both_conserved = []
te_conserved = []
lncrna_conserved = []

for item in gl:
    if item['phyloP_tes'] > t and item['phyloP_nottes'] > t:
        both_conserved.append(item)
    elif item['phyloP_tes'] > t:
        te_conserved.append(item)
    elif item['phyloP_nottes'] > t:
        lncrna_conserved.append(item)
    else:
        not_counted += 1

print('Not counted     : {0:,}'.format(not_counted))
print('Both conserved  : {0:,}'.format(len(both_conserved)))
print('TE conserved    : {0:,}'.format(len(te_conserved)))
print('lncRNA conserved: {0:,}'.format(len(lncrna_conserved)))
