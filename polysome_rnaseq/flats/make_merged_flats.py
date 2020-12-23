
from glbase3 import *


wig_to_flat(['../star/hesc_cyto.rp1.Signal.UniqueMultiple.str1.out.wig.gz',
    '../star/hesc_cyto.rp2.Signal.UniqueMultiple.str1.out.wig.gz'], 'hesc_cyto.flat', 'hESC Cytoplasm', skip_non_standard_chroms=True, gzip=True)

wig_to_flat(['../star/hesc_monosome.rp1.Signal.UniqueMultiple.str1.out.wig.gz',
    '../star/hesc_monosome.rp2.Signal.UniqueMultiple.str1.out.wig.gz'], 'hesc_monosome.flat', 'hESC Monosome', skip_non_standard_chroms=True, gzip=True)

wig_to_flat(['../star/hesc_nuc.rp1.Signal.UniqueMultiple.str1.out.wig.gz',
    '../star/hesc_nuc.rp2.Signal.UniqueMultiple.str1.out.wig.gz'], 'hesc_nuc.flat', 'hESC Nucleus', skip_non_standard_chroms=True, gzip=True)

wig_to_flat(['../star/hesc_poly_high.rp1.Signal.UniqueMultiple.str1.out.wig.gz',
    '../star/hesc_poly_high.rp2.Signal.UniqueMultiple.str1.out.wig.gz'], 'hesc_polyhigh.flat', 'hESC Polysome high', skip_non_standard_chroms=True, gzip=True)

wig_to_flat(['../star/hesc_poly_low.rp1.Signal.UniqueMultiple.str1.out.wig.gz',
    '../star/hesc_poly_low.rp2.Signal.UniqueMultiple.str1.out.wig.gz'], 'hesc_polylow.flat', 'hESC Polysome low', skip_non_standard_chroms=True, gzip=True)

