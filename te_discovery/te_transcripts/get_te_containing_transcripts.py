
'''

Build the final annotation tables;

'''

import glob, sys, os, gzip
from glbase3 import glload, utils, expression, genelist

all_genes = glload('../../transcript_assembly/packed/all_genes.glb')
sam = '../hmmer_dfam/data/transcript_table_HSC_SR_PB_merged.transcripts.glb'
print('... %s' % sam)
doms = glload(sam)

doms = doms.map(all_genes, key='transcript_id')
print(doms.keys())
doms.save('transcript_table_merged.mapped.glb') # fixed name version
doms.saveTSV('transcript_table_merged.mapped.tsv') # fixed name version

# And do the same for the gencode data:
gencode = glload('../../gencode/hg38_gencode_v32.glb')
doms = glload('../hmmer_dfam/data/transcript_table_gencode.v32.transcripts.glb') # to fix;
# I need to go in and fix the transcript_id key, which is just the FASTA name;
for g in doms:
    g['transcript_id'] = g['transcript_id'].split('|')[0].split('.')[0]
    g['enst'] = g['transcript_id']

doms._optimiseData()
doms = doms.map(gencode, key='enst')
print(doms.keys())
doms.save('transcript_table_gencode_all.glb')
doms.saveTSV('transcript_table_gencode_all.tsv')

# And the PC and ncrna lists:
doms = glload('../hmmer_dfam/data/transcript_table_Homo_sapiens.GRCh38.ncrna.glb') # to fix;
# I need to go in and fix the transcript_id key, which is just the FASTA name;
for g in doms:
    g['transcript_id'] = g['transcript_id'].split('|')[0].split('.')[0]
    g['enst'] = g['transcript_id']
doms._optimiseData()
doms = doms.map(gencode, key='enst')
print(doms.keys())
doms.save('transcript_table_gencode_ncrna.glb')
doms.saveTSV('transcript_table_gencode_ncrna.tsv')

doms = glload('../hmmer_dfam/data/transcript_table_gencode.v32.pc_transcripts.glb') # to fix;
# I need to go in and fix the transcript_id key, which is just the FASTA name;
for g in doms:
    g['transcript_id'] = g['transcript_id'].split('|')[0].split('.')[0]
    g['enst'] = g['transcript_id']
doms._optimiseData()
doms = doms.map(gencode, key='enst')
print(doms.keys())
doms.save('transcript_table_gencode_pc.glb')
doms.saveTSV('transcript_table_gencode_pc.tsv')
