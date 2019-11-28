
'''

Draw protein domain-style plots, but containing the locations of the TEs

The plots should be the mRNA, with a protein coding exon (if present), and the location of the TE;

This script collectes the data

'''

import glob, sys, os, gzip
from glbase3 import glload, utils, expression, genelist, genome_sql

sys.path.append('../')
import draw_domains_share
sys.path.append('../../../')
import shared

draw = 'png'

#[os.remove(f) for f in glob.glob('%s/*.%s' % (draw, draw))]

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
gencode_db = genome_sql(filename=os.path.expanduser('~/hg38/hg38_gencode_v29.sql'))
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

genes = set(['SOX2', 'NANOG', 'SALL4', 'LIN28A', 'LIN28B', 'SALL1', 'POU5F1A',
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

    # Wang Jichang paper, Fig 3a. These ones have HERV spliced into their message
    'SCGB3A2', 'NCR1', 'KLKB1', 'IL34', 'PLP1', 'ESRG', 'RPL39L',
    ])

for n, gene in enumerate(doms):
    #print(gene['name'].split(' ')[0])
    if gene['name'].split(' ')[0] not in genes:
        continue

    draw_domains_share.draw_domain(gene, '%s/%s.%s.%s.%s' % (draw, gene['name'], gene['transcript_id'], gene['enst'], draw), gencode_db, dfam)

    if (n+1) % 1000 == 0:
        print('Processed: {:,} domains'.format(n+1))
        #break
