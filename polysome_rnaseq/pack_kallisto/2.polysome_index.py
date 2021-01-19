"""

post EDASeq clean-up and annotation.


"""

import sys, os, glob
from glbase3 import *
config.draw_mode = ['pdf']

expn = glload("kall_tpm-unmerged.glb")

print(expn.getConditionNames())

expn = expn.sliceConditions(
    ['hesc_cyto_rp1', 'hesc_cyto_rp2',
    'hesc_monosome_rp1', 'hesc_monosome_rp2',
    'hesc_poly_low_rp1', 'hesc_poly_low_rp2',
    'hesc_poly_high_rp1', 'hesc_poly_high_rp2']
    )

newe = []
__rejected = 0
for transcript in expn:
    # assume cyto = 100%
    td = transcript['conditions']
    total_cyt = (td[0] + td[1]) / 2
    if total_cyt < 0.1:
        __rejected += 1
        continue
    monosome = (td[2] + td[3]) / 2
    poly_low = (td[4] + td[5]) / 2
    poly_hi = (td[6] + td[7]) / 2

    monosome = monosome / total_cyt
    poly_low = poly_low / total_cyt
    poly_hi = poly_hi / total_cyt

    transcript['conditions'] = [monosome, poly_low, poly_hi]

    newe.append(transcript)

print('Rejected = {} Kept = {}'.format(__rejected, len(newe)))
expn = expression(loadable_list=newe, cond_names=['Monosome', 'Polysome Low', 'Polysome High'])

expn.save('polysome_index.glb')
print(expn)
