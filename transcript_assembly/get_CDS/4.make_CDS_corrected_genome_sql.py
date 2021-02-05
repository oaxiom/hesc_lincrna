
import gzip
from glbase3 import genelist, location, genome, glload

# Combined GTF:
print('LR+SR GTF...')

all_transcripts = glload('../../transcript_assembly/packed/all_genes.glb')
all_transcripts = {gene['transcript_id']:gene for gene in all_transcripts}
data = glload('coding_genes_with_local_CDS-corrected.glb')

newgl = []
for trans in data:
    tid = trans['transcript_id']
    trans_data = all_transcripts[tid]

    #print(trans)
    #print(trans_data)

    if 'cds_gencode_loc' in trans and trans['cds_gencode_loc']:
        cds = trans['cds_gencode_loc']
    elif 'cds_local_to_genome' in trans and trans['cds_local_to_genome']:
        cds = trans['cds_local_to_genome']
    else: # noncoding;
        cds = location(chr=trans_data['loc']['chr'], left=-1, right=-1) # Kind of hacky;

    newt = {
        'loc': trans_data['loc'],
        'exonCounts': trans_data['exonCounts'],
        'exonStarts': trans_data['exonStarts'],
        'exonEnds': trans_data['exonEnds'],
        'name': trans_data['name'],
        'strand': trans_data['strand'],
        'cds_loc': cds,
        }
    print(newt)
    newgl.append(newt)

gl = genome()
gl.load_list(newgl)
gl.save('cds_genome.glb')
gl.saveTSV('cds_genome.tsv')

