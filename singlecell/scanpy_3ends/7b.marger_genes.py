import logging, matplotlib, os, sys
import scanpy as sc
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import colors
import seaborn as sb
from rpy2.robjects.packages import importr
from glbase3 import glload
#from gprofiler import gprofiler
plt.rcParams['figure.figsize']=(8,8) #rescale figures
sc.settings.verbosity = 1
sc.set_figure_params(dpi=200, dpi_save=200)

sc.settings.figdir = 'markers'

transcript_id = glload('../../transcript_assembly/packed/all_genes.glb')

transcript_id_lookup = {}
for ti in transcript_id:
    gene_symbol = ti['name'].split(' ')[0]
    if gene_symbol not in transcript_id_lookup:
        transcript_id_lookup[gene_symbol] = []
    transcript_id_lookup[gene_symbol].append({'name': ti['name'], 'transcript_id': ti['transcript_id']})

adata = sc.read('learned.h5ad') # You can skip the script 3 if using te 2b.

print(adata.var_names)

oh = open('gene_names.all.tsv', 'w')
for g in adata.var_names:
    oh.write('%s\n' % g)
oh.close()
'''
marker_genes_dict = {
    'stemness': [transcript_id['POU5F1'], transcript_id['SOX2'], transcript_id['UTF1'], transcript_id['NANOG'], transcript_id['DNTM3B'], transcript_id['DPPA2']],
    'lineage_markers': [transcript_id['LEFTY2'],
    }

sc.pl.stacked_violin(adata, marker_genes_dict, groupby='leiden_r0.60', rotation=90, dendrogram=True, show=False, save='markers.pdf')
sc.pl.dotplot(adata, marker_genes_dict, groupby='leiden_r0.60', dot_max=0.5, dendrogram=True, standard_scale='var', show=False, save='markers.pdf')
sc.pl.heatmap(adata, marker_genes_dict, groupby='leiden_r0.60', show=False, save='markers.pdf')
'''

genes_to_do = ['POU5F1', 'SOX2', 'UTF1', 'NANOG', 'DPPA2', 'LEFTY2', 'KLF4', 'LIN28A',
    'HNRNPU', 'TFCP2L1', # Naive;
    'EOMES', 'OTX2', 'MIXL1', 'DKK1','SOX17', 'SOX7',  # DE markers
    'SOX1', # EC makers;
    'GATA3', 'GATA2', # ME markers
    # PGC markers;
    'CDH1', 'EPCAM', 'PECAM1', 'OCLN', 'PKP1', 'TCF4', 'CDH11',# E markers
    'CDH2', 'VIM', 'SNAI1', 'SNAI2', 'ZEB1', 'ZEB2', 'TGFBR1', 'TWIST2', 'TWIST1', 'KLF8'# M markers;
    ]

for gene_name in genes_to_do:
    if gene_name not in transcript_id_lookup:
        print('{0} gene name not found in lookup!'.format(gene_name))
        continue

    for transcript in transcript_id_lookup[gene_name]:
        try:
            sc.pl.umap(adata, color=transcript['transcript_id'], size=10, legend_loc='on data',
            vmax=3, show=False, save='markers-{0}-{1}.pdf'.format(transcript['transcript_id'], transcript['name']))
        except KeyError: # this specific transcript_id is missing;
            print('{0} {1} not found'.format(transcript['transcript_id'], transcript['name']))

# These cytokines have very low expression, so move the scale down a bit;
genes_to_do = [
    'WNT4', 'BMP4', 'TGFB1'
    ]

for gene_name in genes_to_do:
    if gene_name not in transcript_id_lookup:
        print('{0} gene name not found in lookup!'.format(gene_name))
        continue

    for transcript in transcript_id_lookup[gene_name]:
        try:
            sc.pl.umap(adata, color=transcript['transcript_id'], size=10, legend_loc='on data',
            vmax=2, show=False, save='markers-{0}-{1}.pdf'.format(transcript['transcript_id'], transcript['name']))
        except KeyError: # this specific transcript_id is missing;
            print('{0} {1} not found'.format(transcript['transcript_id'], transcript['name']))

