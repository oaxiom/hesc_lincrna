'''

Work out the distance of the STOP from the 3' exon-exon junction

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

has_peptide_hit = glload('../results_gene.glb').removeDuplicates('transcript_id') # ones with a peptide hit.
has_peptide_hit = has_peptide_hit.map(genelist=all_transcripts, key='transcript_id')

no_peptide_hit = glload('../2.blast_searches/super_table.glb').removeDuplicates('transcript_id')
no_peptide_hit = has_peptide_hit.map(genelist=no_peptide_hit, key='transcript_id', logic='notright')
no_peptide_hit = no_peptide_hit.map(genelist=all_transcripts, key='transcript_id')

res = {
    'hit': [],
    'no hit': []
    }
    
data = {
    'hit': has_peptide_hit,
    'no hit': no_peptide_hit,
    }

__nostop = 0

for k in data:
    for transcript in data[k]:
        strand = transcript['strand']
        
        cds_key_to_use = None
        if 'cds_genome_loc' in transcript and transcript['cds_genome_loc']:
            cds_key_to_use = 'cds_genome_loc'
        else:
            cds_key_to_use = 'cds_local_to_genome'
        
        if strand == '+':
            STOP = transcript[cds_key_to_use]['right']
            min_dist_to_ejc = transcript['loc']['right'] - STOP
        else:
            STOP = transcript[cds_key_to_use]['left']
            min_dist_to_ejc = STOP - transcript['loc']['left']
        #print(STOP, min_dist_to_ejc, transcript)
        
        if STOP == 0:
            __nostop += 1
            continue

        for exons in zip(transcript['exonStarts'], transcript['exonEnds']):
            if STOP >= exons[0] and STOP <= exons[1]:
                if strand == '+':
                    d = exons[1] - STOP
                else:
                    d = STOP - exons[0]
                print(d)
                min_dist_to_ejc = min([min_dist_to_ejc, d])
        res[k].append(min_dist_to_ejc)

print('__nostop = {}'.format(__nostop))

# percent within 55 bp;
perc_55 = {    
    'hit': [0, 0], # <55, total
    'no hit': [0, 0],
    }
qs = {}
for k in res:
    for bp in res[k]:
        if bp <= 55:
            perc_55[k][0] += 1
        perc_55[k][1] += 1

    qs[k] = (perc_55[k][0] / perc_55[k][1]) * 100.0
    print(k, '{:.2f}'.format((perc_55[k][0] / perc_55[k][1]) * 100.0), perc_55[k])    


shared.boxplots_simple('distances_from_ejc.pdf', res, qs, 
    #xlims=[0, 2500], 
    col='#FF8A87', 
    vlines=[0, 55, 1000],
    )

shared.violin_simple('distances_from_ejc-viol.pdf', res, None, 
    xlims=[0, 2000], 
    col='#FF8A87', 
    vlines=[0, 55])
