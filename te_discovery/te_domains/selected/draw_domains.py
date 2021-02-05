
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

draw = 'pdf'

#[os.remove(f) for f in glob.glob('%s/*.%s' % (draw, draw))]

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
CDSs = glload('../../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

cds_doms = doms.map(genelist=CDSs, key='transcript_id')
print(cds_doms)

genes = set([
    'AC005099.1', 'UQCRHL', 'AC005099.1', 'AC005099.1', # GRP 4;
    'SOX2', 'NANOG', 'SALL4', 'LIN28A', 'LIN28B', 'SALL1', 'POU5F1', 'BRCA1', 'BRCA2',
    'DPPA2', 'DPPA3', 'DPPA5', 'PRDM14', 'JARID2', 'SALL2', 'SALL3', 'TCF3',
    'DNMT3L', 'LEFTY2', 'FGF4', 'NODAL', 'CER1', 'NLRP7', 'SFRP1', 'ZIC2', 'KDR',
    'ZFP42', 'C9ORF135', 'ST6GAL1', 'LRP4', 'MSTO1', 'PRODH', # From Pontis et al., 2019 CSC
    'SFRP2', 'OTX2', 'KLF4', 'KLF5',
    'HNRNPK', 'HNRNPU', 'PMPCB',
    'ABL1',
    'ADGRB1',
    'ADGRG5',
    'AL022318.4',
    'AL390728.4',
    'AL390728.4',
    'AL390728.4',
    'ANKRD36',
    'ARHGAP40',
    'ARHGEF6',
    'ATP6V1FNB',
    'ATP10A',
    'C1orf109',
    'C8orf37',
    'C10orf105',
    'CAND1',
    'CCDC144CP',
    'CDC37L1',
    'CFAP52',
    'CFAP73',
    'CNBD2',
    'COG2',
    'CRAMP1',
    'CREB3L2',
    'CTSLP3',
    'CUL9',
    'CYLD',
    'DNAH5',
    'DOCK3',
    'DRC3',
    'EVC',
    'FAM149B1',
    'FBN3',
    'FBXL4',
    'FHIT',
    'FHOD3',
    'GALNS',
    'GARNL3',
    'GDI2',
    'GOLGA6L4',
    'H6PD',
    'HERC2P2',
    'HERC2P2',
    'HERC2P3',
    'HERC2P9',
    'HMCN2',
    'IGF2R',
    'INF2',
    'KCNK17',
    'LBP',
    'MCC',
    'MCF2L2',
    'METTL2A',
    'MGAT3',
    'NBPF9',
    'NBPF19',
    'NDUFV3',
    'NME4',
    'NOTUM',
    'NPR1',
    'OPRD1',
    'PDK1',
    'PIDD1',
    'PIGG',
    'PPIL3',
    'PTK6',
    'PTPRN2',
    'RAB11FIP5',
    'RAP1GAP2',
    'RDH13',
    'RPL5',
    'SAMD12',
    'SCAF4',
    'SCAND2P',
    'SCN5A',
    'SELENON',
    'SHQ1',
    'SMG1P3',
    'ST8SIA6',
    'TAF4',
    'TEX22',
    'THAP3',
    'TNIP1',
    'TRIOBP',
    'TXLNGY',
    'YPEL1',
    'ZBTB24',
    'ZDHHC24',
    'ZNF44',
    'ZNF415',
    'ZNF788P',
    'ZZEF1',

    ])

for n, gene in enumerate(cds_doms):
    #print(gene['name'].split(' ')[0])
    if gene['name'].split(' ')[0] not in genes:
        continue

    draw_domains_share.draw_domain(gene, '%s/%s.%s.%s.%s' % (draw, gene['name'], gene['transcript_id'], gene['enst'], draw), dfam)

    if (n+1) % 1000 == 0:
        print('Processed: {:,} domains'.format(n+1))
        #break

transcript_ids = [
    # GRP 4:
    'HPSCSR.157359.37','HPSCSR.1172.31', 'HPSCSR.157359.37', 'HPSCSR.157359.37',
    ]

gl = genelist()
gl.load_list([{'transcript_id': i} for i in transcript_ids])
tids = doms.map(genelist=gl, key='transcript_id')

for gene in tids:
    draw_domains_share.draw_domain(gene, '%s/%s.%s.%s.%s' % (draw, gene['name'], gene['transcript_id'], gene['enst'], draw), dfam)

    if (n+1) % 1000 == 0:
        print('Processed: {:,} domains'.format(n+1))

