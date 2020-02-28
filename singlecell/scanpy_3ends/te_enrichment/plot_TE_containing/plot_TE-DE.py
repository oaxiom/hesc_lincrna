import os, sys, glob, numpy
import scanpy as sc
from glbase3 import genelist, glload



for filename in glob.glob('../*.glb'):
    print(filename)
    # Just do a simple annotation of yh the TE-contianig DE genes from the scRNA-seq
    stub = os.path.split(filename)[1].split('.')[0].replace('de_genes-', '')
    de_list = glload(filename)

