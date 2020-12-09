
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
from glbase3 import genelist

def add_entry(trans, gsql, newgl, done):
    transcript_id = trans['transcript_id']

    exon_state = '--'
    if trans['exonCounts'] > 1:
        exon_state = 'ME'
    elif trans['exonCounts'] == 1:
        exon_state = 'SE'

    # Get hte extra features;
    if transcript_id in data_features_lookup:
        e = data_features_lookup[transcript_id]['expn']
        c_nc = data_features_lookup[transcript_id]['coding']
    else:
        print('WARNING: {0} not in data_features_lookup table, skipping'.format(transcript_id))
        e = 'U'
        return False

    new_name = '%s (%s;%s;%s;%s;%s)' % (trans['name'], exon_state, coding_noncoding_map[c_nc], expn_map[e], trans['evidence'], trans['decision'])

    toadd = {'ensg': trans['ensg'], 'enst': trans['enst'], 'name': new_name,
        'gene_symbol': trans['name'], 'loc': trans['loc'],
        'transcript_id': transcript_id,
        'exonCounts': trans['exonCounts'],
        'exonStarts': trans['exonStarts'], 'exonEnds': trans['exonEnds'],
        'strand': trans['strand'],
        'tags': '%s; %s; %s; %s; %s' % (exon_state, c_nc, e, trans['evidence'], trans['decision']),
        'coding': c_nc,
        'expression': e,
        'TPM': trans['TPM'],
        }

    #print(toadd)

    newgl.append(toadd)

    gsql.add_feature(trans['loc'], trans['loc'].pointLeft(), trans['exonCounts'], trans['exonStarts'], trans['exonEnds'], new_name, trans['strand'], 'gene')

    done += 1
    if done % 10000 == 0:
        print('Processed: {:,}'.format(done))

    return done

data_features = genelist(filename='../feature_table/assembly_hPSC_detailed.tsv.gz', format={'force_tsv': True, 'transcript_id': 0, 'coding': 3, 'expn': 2, 'decision': 4}, gzip=True)

# convert ot fast lookups:
data_features_lookup = {i['transcript_id']: i for i in data_features.linearData}

#print(data_expression_data)
#print(data_coding_noncoding)

decision = {'same': '=',
    'different': '~'}

coding_noncoding_map = {'coding': 'C',
    'noncoding': 'NC',
    'NA': 'U'}

expn_map = {'enriched': 'ES+', 'unbiased': 'ES:', 'depleted': 'ES-'}

from glbase3 import *

# Combined GTF:
print('LR+SR GTF...')
gsql = genome_sql(new=True, filename='hg38_hpscs_srlr.sql')
oh = gzip.open('../gtf/current_gtf.gtf.gz', 'rt')

newgl = []

done = 0
skipped = 0
trans = None

for idx, line in enumerate(oh):
    if '#' in line[0]:
        continue
    line = line.strip().split('\t')

    #print(line)

    # I need to assemble the full transcript data first
    if line[2] == 'transcript':
        if trans and trans['exonCounts'] > 0: # Write the previous transcript out
            done = add_entry(trans, gsql, newgl, done)

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

        if 'decision' not in gtf_dec:
            gtf_dec['decision'] = '!'
        else:
            gtf_dec['decision'] = decision[gtf_dec['decision']]

        trans = {'strand': line[6], 'loc': location(chr=line[0], left=line[3], right=line[4]),
            'exonCounts': 0,
            'exonStarts': [],
            'exonEnds': [],
            'decision': gtf_dec['decision'],
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

