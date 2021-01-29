
'''

Pack the GTF into a glb and a flat table;

chr1	StringTie	exon	817373	817712	1000	-	.

gene_id "HSCSR.153";
transcript_id "HSCSR.153.40";
exon_number "5";
reference_id "ENST00000447500";
ref_gene_id "ENSG00000230092";
ref_gene_name "AL669831.4";
evidence "SR";
percent_exon_match "100.0";
percent_splice_match "100.0";
decision "same";
gene_name "AL669831.4";
transcript_name "AL669831.4-201";
ens_gene_id "ENSG00000230092";
ens_transcript_id "ENST00000447500";
ens_exon_id "ENSE00001746491";
ens_exon_overlap_perc "1.0";
ens_transcript_biotype "processed_transcript";

Be careful with this script, it's a bit precarious

'''

import gzip
from glbase3 import genelist, glload, genome_sql

packed_assembly = glload('../../transcript_assembly/packed/all_genes.glb')

corrected_cds = glload('coding_genes_with_local_CDS-corrected.glb')
corrected_cds_lookup = {i['transcript_id']: i for i in corrected_cds.linearData}

# Combined GTF:
print('LR+SR GTF...')
gsql = genome_sql(new=True, filename='hg38_hpscs_srlr_guessed_cds.sql')

newgl = []

for idx, transcript in enumerate(packed_assembly):
    print(transcript)

    gsql.add_feature(transcript['loc'], trans['cds_loc'], trans['exonCounts'], trans['exonStarts'], trans['exonEnds'], trans['name'], trans['strand'], 'gene')

    1/0
    # Start a new transcript;
    gtf_dec = {}
    for item in line[8].split(';'):
        if item:
            item = item.strip(' ').replace('"', '').split(' ')
            if len(item) == 2: # A few bad decorators;
                gtf_dec[item[0]] = item[1]
    #print(gtf_dec)
    if 'GENCODE_gene_id' not in gtf_dec:
        gtf_dec['GENCODE_gene_id'] = gtf_dec['gene_id']
        gtf_dec['gene_name'] = gtf_dec['transcript_id'] # If one is missing, the other is also missing;
        gtf_dec['GENCODE_transcript_id'] = gtf_dec['transcript_id'] # enst

    trans = {'strand': line[6], 'loc': location(chr=line[0], left=line[3], right=line[4]),
        'exonCounts': 0,
        'exonStarts': [],
        'exonEnds': [],
        'transcript_class': gtf_dec['transcript_class'],
        'evidence': gtf_dec['evidence'],
        'name': gtf_dec['gene_name'],
        'ensg': gtf_dec['GENCODE_gene_id'],
        'enst': gtf_dec['GENCODE_transcript_id'],
        'gene_id': gtf_dec['gene_id'],
        'transcript_id': gtf_dec['transcript_id'],
        'evidence': gtf_dec['evidence'],
        'tags': '',
        'coding': '',
        'expression': '',
        'TPM': float(gtf_dec['TPM'])}

    if line[2] == 'exon':
        # Isn't this different, depending upon the strand?
        trans['exonStarts'].append(int(line[3]))
        trans['exonEnds'].append(int(line[4]))
        trans['exonCounts'] += 1

# Add the last entry;
add_entry(trans, gsql, newgl, done)

print('Processed: {:,} transcripts'.format(done))
print('Skipped  : {:,} transcripts'.format(skipped))
oh.close()
gsql.finalise()

gl = genelist()
gl.load_list(newgl)
gl.saveTSV('all_genes.tsv', key_order=['ensg', 'transcript_id', 'name', 'enst'])
gl.save('all_genes.glb')

