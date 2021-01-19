'''

boxplots for polysome, etc.

'''

import sys, os, itertools
import numpy as np
import matplotlib.tri as tri
from glbase3 import *

sys.path.append('../../')
import shared

#expn = glload('../rsem_star/norm_input.glb')
expn = glload('../pack_kallisto/polysome_index.glb')
#expn.log(2,.1)
cond_names = expn.getConditionNames()

res_c = {
    'pc':
        {c: [] for c in cond_names},
    'ncrna':
        {c: [] for c in cond_names},
    }

for idx, cond in enumerate(cond_names):
    for transcript in expn:
        if ';C;' in transcript['name']:
             res_c['pc'][cond].append(transcript['conditions'][idx])

        elif ';NC;' in transcript['name']:
            res_c['ncrna'][cond].append(transcript['conditions'][idx])

cm = np.median(res_c['pc']['Polysome High'])
ncm = np.median(res_c['ncrna']['Polysome High'])

print(cm, ncm)

shared.boxplots_simple('poly-pc.pdf'.format(cond), res_c['pc'], None, xlims=[0, 1], col='#FF8A87', vlines=[0.6195055272311432 0.3118535407277264])
shared.boxplots_simple('poly-ncrna.pdf'.format(cond), res_c['ncrna'], None, xlims=[0, 1], col='#92A7FF', vlines=[0.6195055272311432 0.3118535407277264])



