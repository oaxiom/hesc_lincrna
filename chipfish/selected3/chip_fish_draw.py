

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload, genelist

draw = 'pdf'

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs_3.txt"))

#gllocs = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
gllocs = glload('../../transcript_assembly/packed/all_genes.glb') # needed for noTE

locs = [
    'ERVH48-1',
    'PMPCB',
    'PCAT14',
    'TFCP2L1', # LBP9
    # HERVK-containing
    'LINC02018', 'PDCL3P4',
    'AC108519.1',
    'AF228730.5',
    'AC068587.4',
    'WRAP73',
    'CEP104',
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
    ]
locs = genelist(loadable_list=[{'gene_symbol': k} for k in locs])
ll = locs.map(genelist=gllocs, key='gene_symbol')

print(ll)

for gene in ll:
    scale = 1.0
    if draw == 'pdf':
        scale = 0.3

    c.draw.setLocation(loc=gene['loc'].expand(len(gene['loc']) / 10))
    if 'enst' not in gene:
        print(gene['loc'])
        c.draw.exportImage("{0}/{1}.{2}".format(draw, gene['name'], draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason
    else:
        print(gene['name'])
        c.draw.exportImage("{0}/{1}_{2}_{3}.{4}".format(draw, gene['name'], gene['transcript_id'], gene['enst'], draw), scale=scale, type=draw) # Cannot draw png and svg interleaved for some reason.


