
import sys
from glbase3 import *
sys.path.append('../../')
import shared

# Replace the predictions with confirmed

newd = glload('coding_genes_with_local_CDS-predicted.glb')

# Now go back and put the gencode CDS into the ;= transcripts:
gencode_map = glload('gencode_cds.glb')
gencode_map = {gene['enst']:gene for gene in gencode_map}

gencode_ann = glload('../../gencode/hg38_gencode_v32.pc.glb')
gencode_ann = {gene['enst']:gene for gene in gencode_ann}

all_transcripts = glload('../../transcript_assembly/packed/all_genes.glb')
all_transcripts = {gene['transcript_id']:gene for gene in all_transcripts}

newl = []
for gene in newd:
    enst = gene['enst']
    print(gene)
    # Overwrite if it's available in gencode_ann
    if enst in gencode_ann and ';=)' in gene['name']:
        gene['cds_local_locs'] = gencode_map[enst]['cds_local_locs']
        gene['cds_gencode_loc'] = gencode_ann[enst]['cds_loc'] # Must be by definition!

        #print(gencode_ann[enst])

        gloc = shared.convert_local_to_genome(
            gencode_map[enst]['cds_local_locs'][0], gencode_map[enst]['cds_local_locs'][1],
            gencode_ann[enst]['loc'],
            gencode_ann[enst]['exonStarts'],
            gencode_ann[enst]['exonEnds'],
            gencode_ann[enst]['strand'])

        gene['cds_local_to_genome'] = location(chr=gene['cds_gencode_loc']['chr'], left=gloc[0], right=gloc[1])
        gene['strand'] = gencode_ann[enst]['strand']

    #elif enst in gencode_map and ';=)' in gene['name']:
    #    gene['cds_local_locs'] = gencode_map[enst]['cds_local_locs']
    #    gene['cds_genome_loc'] = gencode_map[enst]['cds_genome']

    elif gene['coding']: # coding, but no cds_genome, you need to figure it out;
        # convert the cds_local_locs to genome_coordinates based on
        # exonStarts and exonEnds;
        tid = gene['transcript_id']
        gene['cds_gencode_loc'] = None # Unknown
        gloc = shared.convert_local_to_genome(
            gene['cds_local_locs'][0], gene['cds_local_locs'][1],
            all_transcripts[tid]['loc'],
            all_transcripts[tid]['exonStarts'],
            all_transcripts[tid]['exonEnds'],
            all_transcripts[tid]['strand'])

        gene['cds_local_to_genome'] = location(chr=all_transcripts[tid]['loc']['chr'], left=gloc[0], right=gloc[1])
        gene['strand'] = all_transcripts[tid]['strand']
    else:
        print(gene)
        1/0

    newl.append(gene) # add the prediction

newd = genelist()
newd.load_list(newl)
newd.save('coding_genes_with_local_CDS-corrected.glb') # Now with actual CDS from GENCODE;
newd.saveTSV('coding_genes_with_local_CDS-corrected.tsv')
