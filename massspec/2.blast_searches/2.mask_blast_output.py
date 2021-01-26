'''
blastp -outfmt 6

 0.	 qseqid	 query (e.g., gene) sequence id
 1.	 sseqid	 subject (e.g., reference genome) sequence id
 2.	 pident	 percentage of identical matches
 3.	 length	 alignment length
 4.	 mismatch	 number of mismatches
 5.	 gapopen	 number of gap openings
 6.	 qstart	 start of alignment in query
 7.	 qend	 end of alignment in query
 8.	 sstart	 start of alignment in subject
 9.	 send	 end of alignment in subject
 10.	 evalue	 expect value
 11.	 bitscore	 bit score

   'qaccver saccver pident length mismatch gapopen qstart qend sstart send
   evalue bitscore', which is equivalent to the keyword 'std'
'''

# This script takes the blastp output, extracts the >95% ~ matches, and then deletes those hits in the
# FASTA peptide sequences,
# leaving a list of unique peptide sequences

import glob, sys, os
from glbase3 import *

[os.remove(f) for f in glob.glob('masked/*.tsv')]
[os.remove(f) for f in glob.glob('masked/*.glb')]

form = {'query_name': 0,
    'hit_name': 1,
    'pident': 2,
    'match_len': 3,
    'qstart': 6,
    'qend': 7,
    'sstart': 8,
    'send': 9,
    'skiplines': -1,
    'force_tsv': True,
    }

min_length = 20

super_table = None

oh_all_seqs = open('all_masked_peptides.fa', 'w')

for filename in glob.glob('blaster/table_*.tsv'):
    if 'table_inframe_insertion.tsv' in filename:
        continue

    stub = os.path.split(filename)[1].replace('.tsv', '')
    blasta = genelist(filename, format=form)

    # Quick lookups;
    blasta_lookup = {}
    for b in blasta:
        if b['query_name'] not in blasta_lookup:
            blasta_lookup[b['query_name']] = []
        blasta_lookup[b['query_name']].append(b)

    # Parse the FASTA
    fasta = genelist('../1.sequences/fasta/{}.fasta'.format(stub), format=format.fasta)
    fasta_lookup = {}

    # A few squeek through somehow. I think these are a few crazy lucky ones that are the same length
    # as the exisiting CDS. Probably a very subtle frameshift.
    if 'class_not_found' in stub:
        continue
    if '_no_disruption_' in stub:
        continue

    res = []
    for f in fasta:
        query_len = len(f['seq'])
        if f['name'] not in blasta_lookup:
            f['blastp_status'] = 'Not Found'
            if len(f['seq']) < min_length:
                print('Warning: {0} <{1} Amino acids, skipping'.format(hit['query_name'], min_length))
            else:
                res.append(f) # Full sequence is uniquq
                oh_all_seqs.write('>{}|{}\n'.format(stub, f['name']))
                oh_all_seqs.write(f['seq'])
                oh_all_seqs.write('\n')

        else: # Was found by BLASTP
            blastp_status = 'Found, but no good hits'
            remaining_sequence = f['seq']

            for hit in blasta_lookup[f['name']]:
                if query_len == hit['match_len'] and hit['pident'] >= 90.0: # basically a 100 % match;
                    print('Warning: {} full length and >90% blast match, skipping'.format(hit['query_name']))
                    remaining_sequence = None
                    break

                if hit['pident'] >= 90.0:
                    # I just keep masking it from all the hits
                    hit_len = ''.join(['n'] * (hit['match_len']-1))
                    remaining_sequence = remaining_sequence[0:hit['qstart']] + hit_len + remaining_sequence[hit['qend']:]
                    #print(hit)
                    #print(f['seq'])
                    #print('{0}{1}{2}'.format(remaining_sequence[0:hit['qstart']], hit_len, remaining_sequence[hit['qend']:]))
                    # Sanity check that it's the correct length:
                    if len(remaining_sequence) != len(f['seq']):
                        1/0
                    blastp_status = 'Found, masked'

            if remaining_sequence:
                # check it could still be found in MS data:
                num_notN = len(remaining_sequence) - remaining_sequence.count('n')
                if num_notN < min_length:
                    print('Warning: {0} masked to <{1} Amino acids, skipping'.format(hit['query_name'], min_length))
                elif len(remaining_sequence) < min_length:
                    print('Warning: {0} <{1} Amino acids, skipping'.format(hit['query_name'], min_length))
                elif True not in [aa in remaining_sequence for aa in ('K', 'R')]: # Lys-C/Trypsin mix cutters
                    print('Warning: {} no Lys-C cleavage sites'.format(hit['query_name'], min_length))
                else:
                    # Sometimes the starting M fails to get masked;
                    if remaining_sequence[0:2] == 'Mn':
                        remaining_sequence = remaining_sequence[1:]
                    res.append({'name': f['name'],
                        'seq': remaining_sequence,
                        'blastp_status': blastp_status})

                    # save a version for DB searching
                    if 'no_disruption' not in stub and 'class_not_found' not in stub:
                        oh_all_seqs.write('>{0}|{1}\n'.format(stub, f['name']))

                        # tr -s 'n', and replace with a '-' to signify a gap, and stop the ability to search back across peptides that have a matching segment(s) in the middle.
                        # Don't do this, and just replace the 'n' with a '-'. This allows me to keep the actual position, and so

                        oh_all_seqs.write(''.join([a for a in remaining_sequence.replace('n', '-') if a not in ('_')]))
                        oh_all_seqs.write('\n')

    if res:
        resgl = genelist()
        resgl.load_list(res)
        print('Number of surviving peptides: {0}'.format(len(resgl)))
        resgl.sort('name')
        resgl.saveTSV('masked/masked_results-{0}.tsv'.format(stub), key_order=['name', 'blastp_status'])
        resgl.save('masked/masked_results-{0}.glb'.format(stub))
        # resgl.saveFASTA

        if super_table:
            super_table += resgl
        else:
            super_table = resgl

newgl = []
for item in super_table:
    n = item['name'].split('|')
    del item['name']
    del item['seq']
    item['name'] = n[0]
    item['transcript_id'] = n[1]
    item['enst'] = n[2]
    newgl.append(item)

super_table.load_list(newgl)
super_table.save('super_table.glb')
super_table.saveTSV('super_table.tsv')

oh_all_seqs.close()
