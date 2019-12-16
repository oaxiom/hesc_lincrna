'''

Add new tables that split the custom gtf into coding and non-coding

'''

import glob, sys, os, gzip
from glbase3 import utils, expression, genelist, glload

gl = glload('transcript_table_HSC_SR_PB_merged.transcripts.glb')
coding_noncoding = genelist(filename='../../../transcript_assembly/coding_noncoding/coding_table.txt.gz', format={'force_tsv': True, 'transcript_id': 0, 'coding': 11}, gzip=True)
all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')

new_c = []
new_nc = []

num_fasta = {'c': 0, 'nc': 0}
num_bps = {'c': 0, 'nc': 0}

for item in gl:
    c_nc = coding_noncoding.get(key='transcript_id', value=item['transcript_id'])
    #print(item)
    if not c_nc:    # 10 are missing for some reason
        continue

    c_nc = c_nc['coding'][0]

    if c_nc == 'coding':
        new_c.append(item)
    elif c_nc == 'noncoding':
        new_nc.append(item)

new_gl = genelist()
new_gl.load_list(new_c)
new_gl.save('transcript_table_HSC_SR_PB_merged.pc.glb')

new_gl = genelist()
new_gl.load_list(new_nc)
new_gl.save('transcript_table_HSC_SR_PB_merged.ncrna.glb')

# Update the summary.tsv with the new sizes:

for c_nc in coding_noncoding:
    #print(item)
    ann = all_genes.get(key='transcript_id', value=c_nc['transcript_id'])[0]
    if not ann:    # 10 are missing for some reason
        continue

    c_nc = c_nc['coding']

    transcript_length = 1
    for e in zip(ann['exonStarts'], ann['exonEnds']):
        transcript_length += (e[1] - e[0]) + 1 # Actual sizes are open;

    #print(ann)
    #print(transcript_length)

    if c_nc == 'Coding':
        num_fasta['c'] += 1
        num_bps['c'] += transcript_length
    elif c_nc == 'Noncoding':
        num_fasta['nc'] += 1
        num_bps['nc'] += transcript_length
    else:
        print('No annot %s' % c_nc)


ih = open('raw_data/summary.tsv', 'r')
oh = open('lib_sizes.txt', 'w')

for line in ih:
    oh.write(line)

ih.close()
oh.write('HSC_SR_PB_merged.ncrna\t%s\t%s\n' % (num_fasta['nc'], num_bps['nc']))
oh.write('HSC_SR_PB_merged.pc\t%s\t%s\n' % (num_fasta['c'], num_bps['c']))



