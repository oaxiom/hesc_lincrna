import sys, os
from glbase3 import *
sys.path.append('../../')
import shared

all_genes = glload('../../transcript_assembly/packed/all_genes.glb')
cds = glload('../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
cds = {gene['transcript_id']: gene for gene in cds}
#print(cds)
tes = glload('../te_transcripts/transcript_table_merged.mapped.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}

newl = []
for transcript in all_genes:
    if transcript['coding'] != 'coding':
        continue
    if '!' not in transcript['tags']:
        continue
    if 'LR;' not in transcript['tags']:
        continue

    # Only novel coding, long read-derived transcripts

    print(transcript)

    # Add the doms key to all transcripts;
    if transcript['transcript_id'] in tes:
        te = tes[transcript['transcript_id']]
        transcript['doms'] = te['doms']

    transcript['cds_local_locs'] = cds[transcript['transcript_id']]['cds_local_locs']
    transcript['cds_info'] = 'TRUE'

    newl.append(transcript)

gl = genelist()
gl.load_list(newl)
gl.saveTSV('table_novel_coding.tsv')
gl.save('table_novel_coding.glb')

print('Novel coding passing:', len(gl))

