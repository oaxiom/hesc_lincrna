
'''

Draw protein domain-style plots, but containing the locations of the TEs

The plots should be the mRNA, with a protein coding exon (if present), and the location of the TE;

This version will draw all canonical domains, even if they have no TE. 
From GENCODE;

'''

import glob, sys, os, gzip
from glbase3 import glload, utils, expression, genelist, genome_sql

sys.path.append('../')
import draw_domains_share
sys.path.append('../../../')
import shared

draw = 'pdf'

if not os.path.isdir(draw):
    os.mkdir(draw)
    
[os.remove(f) for f in glob.glob('%s/*.%s' % (draw, draw))]

GENCODE = glload('../../../gencode/hg38_gencode_v32.glb')
CDSs = glload('../../../transcript_assembly/get_CDS/gencode_cds.glb')
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

cds_doms = GENCODE.map(genelist=CDSs, key='enst')
print(cds_doms)

genes = set([
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

