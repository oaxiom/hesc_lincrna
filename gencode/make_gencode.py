
'''

Pack the GTF into a gneome_sql for chipFish

'''

import sys, os, gzip
from glbase3 import *

def add_entry(done ):
    if '.' in trans['loc']['chr']: # strange contigs
        skipped += 1
        return

    #print(trans)
    #1/0

    if trans['cds_loc']:
        trans['cds_loc'] = location(chr=trans['loc']['chr'], left=trans['cds_loc'][0], right=trans['cds_loc'][1])
    else: # A non-coding, so fill in with the TSS instead:
        if trans['strand'] == '+':
            trans['cds_loc'] = trans['loc'].pointLeft()
        else:
            trans['cds_loc'] = trans['loc'].pointRight()

    entry = {'ensg': trans['ensg'],
        'enst': trans['enst'],
        'name': trans['name'],
        'loc': trans['loc'],
        'cds_loc': trans['cds_loc'],
        'transcript_id': trans['transcript_id'],
        'strand': trans['strand'],
        'exonStarts': trans['exonStarts'],
        'exonEnds': trans['exonEnds'],
        'exonCounts': trans['exonCounts'],
        'transcript_type': trans['transcript_type'],
        }

    if trans['transcript_type'] not in transcript_types_dont_keep:
        transcript_types_seen.add(trans['transcript_type'])

        newgl.append(entry)
        gsql.add_feature(trans['loc'], trans['cds_loc'], trans['exonCounts'], trans['exonStarts'], trans['exonEnds'], '%s (%s)' % (trans['name'], trans['enst'].split('.')[0]), trans['strand'], 'gene')

        if trans['transcript_type'] == 'protein_coding':
            pc.append(entry)
        elif trans['transcript_type'] == 'lncRNA':
            ncrna.append(entry)

    done += 1
    if done % 10000 == 0:
        print('Processed: {:,}'.format(done))

    return done

version = 32

# Combined GTF:
print('GENCODE GTF...')
gsql = genome_sql(new=True, filename='hg38_gencode_v%s.sql' % version)
oh = gzip.open('gencode.v%s.annotation.gtf.gz' % version, 'rt') # you need to download this one

newgl = []

pc = []
ncrna = []

done = 0
skipped = 0
trans = None
transcript_types_seen = set([])
transcript_types_dont_keep = set(['miRNA', 'snRNA', 'snoRNA', 'Mt_tRNA', 'Mt_rRNA', 'TEC', 'scRNA', 'scaRNA','vaultRNA',
    'ribozyme',
    'IG_pseudogene', 'IG_J_gene', 'TR_J_pseudogene', 'IG_J_pseudogene', 'TR_D_gene',
    'IG_V_gene', 'IG_C_pseudogene',  'TR_J_gene','TR_V_pseudogene','sRNA', 'TR_V_gene', 'IG_D_gene', 'IG_C_gene', 'IG_V_pseudogene', 'TR_C_gene',])

for idx, line in enumerate(oh):
    if '#' in line[0]:
        continue
    line = line.strip().split('\t')

    # I need to assemble the full transcript data first
    if line[2] == 'transcript':
        if trans and trans['exonCounts'] > 0: # Write the previous transcript out
            done = add_entry(done)

        # Start a new transcript;
        gtf_dec = {}
        for item in line[8].split(';'):
            if item:
                item = item.strip(' ').replace('"', '').split(' ')
                gtf_dec[item[0]] = item[1]

        trans = {'strand': line[6],
            'loc': location(chr=line[0], left=line[3], right=line[4]),
            'cds_loc': None,
            'exonCounts': 0,
            'exonStarts': [],
            'exonEnds': [],
            'name': gtf_dec['gene_name'],
            'ensg': gtf_dec['gene_id'],
            'enst': gtf_dec['transcript_id'],
            'transcript_id': gtf_dec['transcript_id'],
            'transcript_type': gtf_dec['transcript_type'],
            }

    if line[2] == 'exon':
        # Isn't this different, depending upon the strand?
        trans['exonStarts'].append(int(line[3]))
        trans['exonEnds'].append(int(line[4]))
        trans['exonCounts'] += 1

    if line[2] == 'CDS':
        if not trans['cds_loc']:
            trans['cds_loc'] = [1e20, -1]
        trans['cds_loc'] = [min(int(line[3]), trans['cds_loc'][0]), max(int(line[4]), trans['cds_loc'][1])]
        #print(trans['cds_loc'], trans)

done = add_entry(done) # add the last one

print(transcript_types_seen)

print('Processed: %s transcripts' % done)
print('Skipped  : %s transcripts' % skipped)
oh.close()
gsql.finalise()

gl = genelist()
gl.load_list(newgl)
gl.saveTSV('hg38_gencode_v%s.tsv' % version, key_order=['ensg', 'transcript_id', 'name', 'enst'])
gl.save('hg38_gencode_v%s.glb' % version)

gl = genelist()
gl.load_list(pc)
gl.saveTSV('hg38_gencode_v%s.pc.tsv' % version, key_order=['ensg', 'transcript_id', 'name', 'enst'])
gl.save('hg38_gencode_v%s.pc.glb' % version)

gl = genelist()
gl.load_list(ncrna)
gl.saveTSV('hg38_gencode_v%s.ncrna.tsv' % version, key_order=['ensg', 'transcript_id', 'name', 'enst'])
gl.save('hg38_gencode_v%s.ncrna.glb' % version)
