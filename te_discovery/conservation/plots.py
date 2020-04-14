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

shared_conservation.scat('scat_te_vs_not_tes.pdf',
    gl['phyloP_tes'], gl['phyloP_nottes'],
    'TE', 'not-TE',
    xlims=[-0.6, 0.7],
    ylims=[-0.6, 0.7],
    )

shared_conservation.hist('hist_te_vs_not_tes.pdf',
    gl['phyloP_tes'], gl['phyloP_nottes'],
    'TE', 'not-TE',
    ranges=[[-0.6, 0.7], [-0.6, 0.7]],
    )

shared_conservation.contour('cont_te_vs_not_tes.pdf',
    gl['phyloP_tes'], gl['phyloP_nottes'],
    'TE', 'not-TE',
    ranges=[[-0.6, 0.7], [-0.6, 0.7]],
    )

for t in ('phyloP_tes', 'phyloP_nottes'):
    shared_conservation.scat(filename='scat_expn_cons_vs_{0}.pdf'.format(t),
        x=gl[t],
        y=np.log2(np.array(gl['TPM'])+0.1),
        xlabel='cons',
        ylabel='expn',
        xlims=[-0.6, 0.7],
        ylims=[-3, 9],
        )

    shared_conservation.hist(filename='hist_expn_cons_vs_{0}.pdf'.format(t),
        x=gl[t],
        y=np.log2(np.array(gl['TPM'])+0.1),
        xlabel='cons',
        ylabel='expn',
        ranges=[[-0.6, 0.7], [-3, 9]],
        )

    shared_conservation.contour(filename='cont_expn_cons_vs_{0}.pdf'.format(t),
        x=gl[t],
        y=np.log2(np.array(gl['TPM'])+0.1),
        xlabel='cons',
        ylabel='expn',
        ranges=[[-0.6, 0.7], [-3, 9]],
        vmax=100,
        )
