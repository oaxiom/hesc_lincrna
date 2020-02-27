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

transcript_id = glload('../../transcript_assembly/packed/all_genes.glb')

transcript_id_lookup = {}
for ti in transcript_id:
    gene_symbol = ti['name'].split(' ')[0]
    if gene_symbol not in transcript_id_lookup:
        transcript_id_lookup[gene_symbol] = []
    transcript_id_lookup[gene_symbol].append({'name': ti['name'], 'transcript_id': ti['transcript_id']})

adata = sc.read('learned.h5ad')

genes_to_do_g1s= [
    # G1/S
    'MCM5',
    'PCNA',
    'TYMS',
    'FEN1',
    'MCM2',
    'MCM4',
    'RRM1',
    'UNG',
    'GINS2',
    'MCM6',
    'CDCA7',
    'DTL',
    'PRIM1',
    'UHRF1',
    'MLF1IP',
    'HELLS',
    'RFC2',
    'RPA2',
    'NASP',
    'RAD51AP1',
    'GMNN',
    'WDR76',
    'SLBP',
    'CCNE2',
    'UBR7',
    'POLD3',
    'MSH2',
    'ATAD2',
    'RAD51',
    'RRM2',
    'CDC45',
    'CDC6',
    'EXO1',
    'TIPIN',
    'DSCC1',
    'BLM',
    'CASP8AP2',
    'USP1',
    'CLSPN',
    'POLA1',
    'CHAF1B',
    'BRIP1',
    'E2F8',
    ]

sc.settings.figdir = 'cell_cycle_g1s'

for gene_name in genes_to_do_g1s:
    if gene_name not in transcript_id_lookup:
        print('{0} gene name not found in lookup!'.format(gene_name))
        continue

    for transcript in transcript_id_lookup[gene_name]:
        try:
            sc.pl.umap(adata, color=transcript['transcript_id'], size=30, legend_loc='on data', vmax=3, show=False, save='markers-{0}-{1}.pdf'.format(transcript['transcript_id'], transcript['name']))
        except KeyError: # this specific transcript_id is missing;
            print('{0} {1} not found'.format(transcript['transcript_id'], transcript['name']))
        #sc.pl.umap(adata, color=marker_genes_dict[k], color_map='plasma', size=10, vmax=3, legend_loc='on data', show=False, save='markers-{0}.pdf'.format(k))


genes_to_do_g2m = [
    'HMGB2',
    'CDK1',
    'NUSAP1',
    'UBE2C',
    'BIRC5',
    'TPX2',
    'TOP2A',
    'NDC80',
    'CKS2',
    'NUF2',
    'CKS1B',
    'MKI67',
    'TMPO',
    'CENPF',
    'TACC3',
    'FAM64A',
    'SMC4',
    'CCNB2',
    'CKAP2L',
    'CKAP2',
    'AURKB',
    'BUB1',
    'KIF11',
    'ANP32E',
    'TUBB4B',
    'GTSE1',
    'KIF20B',
    'HJURP',
    'HJURP',
    'CDCA3',
    'HN1',
    'CDC20',
    'TTK',
    'CDC25C',
    'KIF2C',
    'RANGAP1',
    'NCAPD2',
    'DLGAP5',
    'CDCA2',
    'CDCA8',
    'ECT2',
    'KIF23',
    'HMMR',
    'AURKA',
    'PSRC1',
    'ANLN',
    'LBR',
    'CKAP5',
    'CENPE',
    'CTCF',
    'NEK2',
    'G2E3',
    'GAS2L3',
    'CBX5',
    'CENPA',
    ]

sc.settings.figdir = 'cell_cycle_g2m'

for gene_name in genes_to_do_g2m:
    if gene_name not in transcript_id_lookup:
        print('{0} gene name not found in lookup!'.format(gene_name))
        continue

    for transcript in transcript_id_lookup[gene_name]:
        try:
            sc.pl.umap(adata, color=transcript['transcript_id'], size=30, legend_loc='on data', vmax=3, show=False, save='markers-{0}-{1}.pdf'.format(transcript['transcript_id'], transcript['name']))
        except KeyError: # this specific transcript_id is missing;
            print('{0} {1} not found'.format(transcript['transcript_id'], transcript['name']))
        #sc.pl.umap(adata, color=marker_genes_dict[k], color_map='plasma', size=10, vmax=3, legend_loc='on data', show=False, save='markers-{0}.pdf'.format(k))
