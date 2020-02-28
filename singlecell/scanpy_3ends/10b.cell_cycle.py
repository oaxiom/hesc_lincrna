import logging, matplotlib, os, sys, glob
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

[os.remove(f) for f in glob.glob('cell_cycle_g1s/*.pdf')]
[os.remove(f) for f in glob.glob('cell_cycle_g2m/*.pdf')]

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
    'ABCC5', 'ABHD10', 'ANKRD18A', 'ASF1B', 'ATAD2', 'BBS2',
    'BIVM', 'BLM', 'BMI1', 'BRCA1', 'BRIP1', 'C5orf42', 'C11orf82',
    'CALD1', 'CALM2', 'CASP2', 'CCDC14', 'CCDC84', 'CCDC150',
    'CDC7', 'CDC45', 'CDCA5', 'CDKN2AIP', 'CENPM', 'CENPQ',
    'CERS6', 'CHML', 'COQ9', 'CPNE8', 'CREBZF', 'CRLS1',
    'DCAF16', 'DEPDC7', 'DHFR', 'DNA2', 'DNAJB4', 'DONSON',
    'DSCC1', 'DYNC1LI2', 'E2F8', 'EIF4EBP2', 'ENOSF1', 'ESCO2',
    'EXO1', 'EZH2', 'FAM178A', 'FANCA', 'FANCI', 'FEN1', 'GCLM',
    'GOLGA8A', 'GOLGA8B', 'H1F0', 'HELLS', 'HIST1H2AC', 'HIST1H4C',
    'INTS7', 'KAT2A', 'KAT2B', 'KDELC1', 'KIAA1598', 'LMO4', 'LYRM7',
    'MAN1A2', 'MAP3K2', 'MASTL', 'MBD4', 'MCM8', 'MLF1IP', 'MYCBP2',
    'NAB1', 'NEAT1', 'NFE2L2', 'NRD1', 'NSUN3', 'NT5DC1', 'NUP160',
    'OGT', 'ORC3', 'OSGIN2', 'PHIP', 'PHTF1', 'PHTF2', 'PKMYT1',
    'POLA1', 'PRIM1', 'PTAR1', 'RAD18', 'RAD51', 'RAD51AP1',
    'RBBP8', 'REEP1', 'RFC2', 'RHOBTB3', 'RMI1', 'RPA2', 'RRM1',
    'RRM2', 'RSRC2', 'SAP30BP', 'SLC38A2', 'SP1', 'SRSF5', 'SVIP',
    'TOP2A', 'TTC31', 'TTLL7', 'TYMS', 'UBE2T', 'UBL3', 'USP1',
    'ZBED5', 'ZWINT',
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
   'ANLN', 'AP3D1', 'ARHGAP19', 'ARL4A', 'ARMC1', 'ASXL1', 'ATL2', 'AURKB',
   'BCLAF1', 'BORA', 'BRD8', 'BUB3', 'C2orf69', 'C14orf80', 'CASP3', 'CBX5',
   'CCDC107', 'CCNA2', 'CCNF', 'CDC16', 'CDC25C', 'CDCA2', 'CDCA3', 'CDCA8',
   'CDK1', 'CDKN1B', 'CDKN2C', 'CDR2', 'CENPL', 'CEP350', 'CFD', 'CFLAR',
   'CHEK2', 'CKAP2', 'CKAP2L', 'CYTH2', 'DCAF7', 'DHX8', 'DNAJB1', 'ENTPD5',
   'ESPL1', 'FADD', 'FAM83D', 'FAN1', 'FANCD2', 'G2E3', 'GABPB1', 'GAS1',
   'GAS2L3', 'H2AFX', 'HAUS8', 'HINT3', 'HIPK2', 'HJURP', 'HMGB2', 'HN1',
   'HP1BP3', 'HRSP12', 'IFNAR1', 'IQGAP3', 'KATNA1', 'KCTD9', 'KDM4A',
   'KIAA1524', 'KIF5B', 'KIF11', 'KIF20B', 'KIF22', 'KIF23', 'KIFC1',
   'KLF6', 'KPNA2', 'LBR', 'LIX1L', 'LMNB1', 'MAD2L1', 'MALAT1', 'MELK',
   'MGAT2', 'MID1', 'MIS18BP1', 'MND1', 'NCAPD3', 'NCAPH', 'NCOA5',
   'NDC80', 'NEIL3', 'NFIC', 'NIPBL', 'NMB', 'NR3C1', 'NUCKS1',
   'NUMA1', 'NUSAP1', 'PIF1', 'PKNOX1', 'POLQ', 'PPP1R2', 'PSMD11',
   'PSRC1', 'RANGAP1', 'RCCD1', 'RDH11', 'RNF141', 'SAP30', 'SKA3',
   'SMC4', 'STAT1', 'STIL', 'STK17B', 'SUCLG2', 'TFAP2A', 'TIMP1', 'TMEM99',
   'TMPO', 'TNPO2', 'TOP2A', 'TRAIP', 'TRIM59', 'TRMT2A', 'TTF2', 'TUBA1A',
   'TUBB', 'TUBB2A', 'TUBB4B', 'TUBD1', 'UACA', 'UBE2C', 'VPS25', 'VTA1',
   'WSB1', 'ZNF587', 'ZNHIT2',
    ]

sc.settings.figdir = 'cell_cycle_g2m'

for gene_name in genes_to_do_g2m:
    if gene_name not in transcript_id_lookup:
        print('{0} gene name not found in lookup!'.format(gene_name))
        continue

    for transcript in transcript_id_lookup[gene_name]:
        try:
            sc.pl.umap(adata, color=transcript['transcript_id'], size=10, legend_loc='on data', vmax=3, show=False, save='markers-{0}-{1}.pdf'.format(transcript['transcript_id'], transcript['name']))
        except KeyError: # this specific transcript_id is missing;
            print('{0} {1} not found'.format(transcript['transcript_id'], transcript['name']))
        #sc.pl.umap(adata, color=marker_genes_dict[k], color_map='plasma', size=10, vmax=3, legend_loc='on data', show=False, save='markers-{0}.pdf'.format(k))
