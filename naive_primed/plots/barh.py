
import sys, os
import matplotlib.cm as cm
from glbase3 import *
config.draw_mode = "png"

annot = glload(os.path.expanduser('~/hg38/hg38_ensembl_v95_enst.glb'))
gllocs = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
expn = glload("../kallisto/kall_tpm-unmerged.glb")

#tree = expn.tree(filename="tree.png", row_names=expn["name"], color_threshold=0.0)

gene_list = ['SOX2', 'NANOG', 'SALL4', 'LIN28A', 'LIN28B', 'SALL1', 'POU5F1A',
    'DPPA2', 'DPPA3', 'DPPA5A', 'PRDM14', 'JARID2', 'SALL2', 'SALL3', 'TCF3',
    'ZFP42', 'C9ORF135', 'ST6GAL1', 'LRP4', 'MSTO1', 'PRODH',# From Pontis et al., 2019 CSC
    'ESRRB', 'LIN28A', 'LIN28B', 'PRDM14',
    'POU5F1', 'SOX2', 'NANOG', 'NCOR1', 'NCOR2', 'SALL1', 'KLF4', 'SALL1', 'NR5A1', 'NR5A2', 'NR5A3',
    'KLF2', 'KLF5', 'LEFTY1', 'LEFTY2', 'FGF4', 'NODAL',
    # Naive-specific genes;
    'ESRRB',  'TFCP2L1', 'ZFP42', 'MT1H', 'DPPA3', 'DPPA4', 'DPPA5', 'ZNF486', 'CR1L', 'DNMT3L', 'ZNF534',
    # Diffenretiation genes;
    'GATA2', 'GATA3', 'GATA4', 'SOX17', 'CER1',
    # 2C genes
    'NR0B1', 'CDX2', 'DUXF3',
    # Down in naive:
    'SFRP1', 'ZIC2', 'KDR', 'OTX2', 'DUSP6', 'SPRY4', 'THY1', 'ID3', 'ZIC5',
    # MA Gang's possibles:
    'HNRNPK', 'DDX1', 'DDX50', 'BRCA2', 'BRCA1', 'TOP1', 'RAP3', 'TRIM25', 'HNRNPU',
    ]

locs = genelist(loadable_list=[{'name': k} for k in gene_list])
ll = locs.map(genelist=annot, key='name')
ll = ll.map(genelist=gllocs, key='enst')

for g in ll:
    expn.barh_single_item(value=g['transcript_id'], key="transcript_id", filename="barh/%s-%s.png" % (g['name'], g['transcript_id']),
        size=(4, 7), vert_space=0.8, #tree=tree,
        yticklabel_fontsize=6, xticklabel_fontsize=6)

gene_list = ['HSCSR.31106.2',

    ]

for g in gene_list:
    expn.barh_single_item(value=g, key="transcript_id", filename="barh/%s.png" % g,
        size=(4, 7), vert_space=0.8, #tree=tree,
        yticklabel_fontsize=6, xticklabel_fontsize=6)

