

import sys, os
sys.path.append(os.path.join(os.path.expanduser("~"), "chipfish"))
import app as chipfish
from glbase_wrapper import location, glload, genelist

draw = 'pdf'

c = chipfish.app()
c.startup(os.path.expanduser("../trk_TEs.txt"))

#gllocs = glload('../../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')
gllocs = glload('../../transcript_assembly/packed/all_genes.glb') # needed for noTE

locs = [
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


