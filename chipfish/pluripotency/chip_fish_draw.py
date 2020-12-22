

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload, genelist

draw = 'pdf'

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

#annot = glload(os.path.expanduser('~/hg38/hg38_ensembl_v95_enst.glb'))
#annot = annot.renameKey('name', 'gene_symbol')
gllocs = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')

locs = ['SOX2', 'NANOG', 'SALL4', 'LIN28A', 'LIN28B', 'SALL1', 'POU5F1A',
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
    # Headline genes from Ihry et al., Cell Rep.
    # Significantly down-regualte POU5F1 P<0.05
    'TGFBR2',
    'GGCT',
    'MAML2',
    'POU5F1',
    'TGFBR1',
    'TMEM107',
    'ZNF469',
    'SMARCA4',
    'PROK2',
    'PAQR7',
    'MINDY4',
    # Odd stuff:
    'LIN28B-AS1',
    # Wang Jichang paper, Fig 3a. These ones have HERV spliced into their message
    'SCGB3A2', 'NCR1', 'KLKB1', 'IL34', 'PLP1', 'ESRG', 'RPL39L',

    ]
locs = genelist(loadable_list=[{'gene_symbol': k} for k in locs])
ll = locs.map(genelist=gllocs, key='gene_symbol')

print(ll)

for gene in ll:
    print(gene['name'])
    c.draw.setLocation(loc=gene['loc'].expand(len(gene['loc']) / 10))
    scale = 1.0
    if draw == 'svg':
        scale = 0.3
    c.draw.exportImage("%s/%s_%s.%s" % (draw, gene['name'], gene['transcript_id'], draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.
