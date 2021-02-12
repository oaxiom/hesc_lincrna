'''

Work out the average CDS sizes;

'''

import glob, sys, os, numpy
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../')
import shared
config.draw_mode = 'pdf'

res = {}

all_transcripts = glload('../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
all_exons = glload('../../transcript_assembly/packed/all_genes.glb')

all_transcripts = all_transcripts.map(key='transcript_id', genelist=all_exons)
all_expressed = all_transcripts.getColumns(['enst', 'transcript_id'])
all_coding = glload('../../gencode/hg38_gencode_v32.pc.glb').map(genelist=all_transcripts, key='enst')
print(len(all_coding))

has_peptide_hit = glload('../results_gene.glb').removeDuplicates('transcript_id') # ones with a peptide hit.
has_peptide_hit = has_peptide_hit.map(genelist=all_transcripts, key='transcript_id')

no_peptide_hit = glload('../2.blast_searches/super_table.glb').removeDuplicates('transcript_id')
no_peptide_hit = has_peptide_hit.map(genelist=no_peptide_hit, key='transcript_id', logic='notright')
no_peptide_hit = no_peptide_hit.map(genelist=all_transcripts, key='transcript_id')

res = {
    'Canonical': [],
    'MS Hit': [],
    'No MS hit': [],

    }
    
data = {
    'Canonical': all_coding,
    'MS Hit': has_peptide_hit,
    'No MS hit': no_peptide_hit,

    }

__nostop = 0
__notfound = 0

for k in data:
    for transcript in data[k]:
        orf_size = transcript['cds_local_locs'][1] - transcript['cds_local_locs'][0]
            
        res[k].append(orf_size)

print('__nostop = {}'.format(__nostop))
print('__notfound = {}'.format(__notfound))

shared.boxplots_simple('cds_sizes-box.pdf', res, None, 
    #xlims=[0, 2500], 
    col='#FF8A87', 
    vlines=[0, 55, 1000],
    )

shared.violin_simple('cds_sizes-viol.pdf', res, None, 
    xlims=[0, 2000], 
    col='#FF8A87', 
    vlines=[0, 55])
