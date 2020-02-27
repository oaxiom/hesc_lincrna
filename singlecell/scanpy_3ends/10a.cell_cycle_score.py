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

# Cell cycle genes are taken from: Tirosh et al, 2015.

sc.settings.figdir = 'cell_cycle_scoring'

genes_to_do_g1s= [
    # G1/S
    'MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1',
    'UNG','GINS2','MCM6','CDCA7','DTL','PRIM1','UHRF1','MLF1IP',
    'HELLS','RFC2','RPA2','NASP','RAD51AP1','GMNN','WDR76','SLBP',
    'CCNE2','UBR7','POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45',
    'CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2','USP1','CLSPN',
    'POLA1','CHAF1B','BRIP1','E2F8',
    ]

genes_to_do_g2m = [
    'HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A',
    'NDC80','CKS2','NUF2','CKS1B','MKI67','TMPO','CENPF',
    'TACC3','FAM64A','SMC4','CCNB2','CKAP2L','CKAP2','AURKB',
    'BUB1','KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP',
    'HJURP','CDCA3','HN1','CDC20','TTK','CDC25C','KIF2C','RANGAP1',
    'NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA',
    'PSRC1','ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3',
    'GAS2L3','CBX5','CENPA',
    ]

# I need to collect all the actual IDs in use:
scrnaseq_data_valid_tids = set(adata.var_names)

tids_g1s = []
for gene_name in genes_to_do_g1s:
    toadd = []
    for tid in transcript_id_lookup[gene_name]:
        if tid in scrnaseq_data_valid_tids:
            toadd.append(tid)
    if toadd:
        tids_g1 += toadd

tids_g2m = []
for gene_name in genes_to_do_g2m:
    toadd = []
    for tid in transcript_id_lookup[gene_name]:
        if tid in scrnaseq_data_valid_tids:
            toadd.append(tid)
    if toadd:
        tids_g2m += toadd

sc.tl.score_genes_cell_cycle(adata, s_genes=tids_g1s, g2m_genes=tids_g2m)

adata_cc_genes = adata[:, tids_g1s+tids_g2m]
sc.tl.pca(adata_cc_genes)
sc.pl.pca_scatter(adata_cc_genes, color='phase')
