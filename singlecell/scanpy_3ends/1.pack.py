"""

Pack the scRNA-seq data using scanpy, prep for scran normalisation

"""

import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from anndata import AnnData
from matplotlib import rcParams
from matplotlib import colors
plt.rcParams['figure.figsize'] = (8,8)
sc.settings.verbosity = 3
sc.set_figure_params(dpi=200, dpi_save=200)
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['font.size'] = 10
sc.settings.autoshow = False

from sc_utils import sparsify

sam1 = sparsify("../te_count_ends/ss.Hs_psc_c11.rp1.tsv.gz", csv=False)           ; sam1.obs['cell_type'] = "iPSC" ; sam1.obs['replicate'] = "c11#1"
sam2 = sparsify("../te_count_ends/ss.Hs_psc_wibr3.rp1.tsv.gz", csv=False)      ; sam2.obs['cell_type'] = "hESC" ; sam2.obs['replicate'] = "WIBR3#1"
#sam3 = sparsify("../te_count_ends/ss.Hs_psc_wibr3nai.rp1.tsv.gz", csv=False)      ; sam3.obs['cell_type'] = "hESC"  ; sam3.obs['replicate'] = "WIBR3-naive#1"
sam4 = sparsify("../te_count_ends/ss.hIPSC_scRNA_Sample1.tsv.gz", csv=False)   ; sam4.obs['cell_type'] = "iPSC" ; sam4.obs['replicate'] = "WTC#1"
#sam5 = sparsify("../te_count_ends/ss.hIPSC_scRNA_Sample2.tsv.gz", csv=False)   ; sam5.obs['cell_type'] = "iPSC-primed" ; sam5.obs['replicate'] = "WTC#2" # This is the one they sequenced a few cells very deep
sam6 = sparsify("../te_count_ends/ss.hIPSC_scRNA_Sample3.tsv.gz", csv=False)   ; sam6.obs['cell_type'] = "iPSC" ; sam6.obs['replicate'] = "WTC#3"
sam7 = sparsify("../te_count_ends/ss.hIPSC_scRNA_Sample4.tsv.gz", csv=False)   ; sam7.obs['cell_type'] = "iPSC" ; sam7.obs['replicate'] = "WTC#4"
sam8 = sparsify("../te_count_ends/ss.hIPSC_scRNA_Sample5.tsv.gz", csv=False)   ; sam8.obs['cell_type'] = "iPSC" ; sam8.obs['replicate'] = "WTC#5"

print('Loaded Samples...')

# Do very simple prefiltering:
samples = [sam1, sam2,
    #sam3,
    sam4,
    #sam5,
    sam6, sam7, sam8]

# Quick pre-filtering, these should be low, otherwise it can mess up downstream analysis, but also can get rid of trivial uninteresting things
[sc.pp.filter_cells(sam, min_genes=500) for sam in samples]
[sc.pp.filter_cells(sam, max_counts=200000) for sam in samples]
[sc.pp.filter_cells(sam, min_counts=2000) for sam in samples]
# Do not filter gene here; concatenate joins on the union, so if a gene fails in a single sample, it will also be deleted from all other samples;

print('Concatenating')
adata = sam1.concatenate(samples[1:])

adata.X = adata.X.astype('float32')

sc.pp.filter_genes(adata, min_cells=1) # Keep this very low, just to remove trivial failures

print('Total number of cells: {:d}'.format(adata.n_obs))
print('Total number of genes: {:d}'.format(adata.n_vars))

adata.write('./raw_data.h5ad')

oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()

