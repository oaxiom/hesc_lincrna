import logging, matplotlib, os, sys, glob, numpy, gzip
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import pandas as pd
from glbase3 import genelist

adata = sc.read('../../learned.h5ad')
print(adata)

# Getting farcical, just dump the table so you can figure out what the hell is going on.
oh = gzip.open('gene_names.filtered.tsv.gz', 'wt')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()


# rows = genes; cols=cells;
print('Save dense matrix')
data_mat = adata.raw.X.T.tocsr()
# Save out a dense array of the sparse array:
oh = gzip.open('dense_array.tsv.gz', 'wt')
for i in range(data_mat.shape[0]):
    if i % 1000 == 0:
        print('{0}/{1}'.format(i, data_mat.shape[0]))

    oh.write('\t'.join([str(i) for i in data_mat.getrow(i).toarray()[0,:] if i > 0.0]))
    oh.write('\n')
oh.close()
