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

oh_all_seqs = open('all_masked_peptides.fa', 'w')

for filename in glob.glob('blaster/table_*.tsv'):
    stub = os.path.split(filename)[1].replace('.tsv', '')
    blasta = genelist(filename, format=form)

    # Quick lookups;
    blasta_lookup = {}
    for b in blasta:
        if b['query_name'] not in blasta_lookup:
            blasta_lookup[b['query_name']] = []
        blasta_lookup[b['query_name']].append(b)

    # Parse the FASTA
    fasta = genelist('../1.fasta/fasta/{0}.fasta'.format(stub), format=format.fasta)
    fasta_lookup = {}

    res = []
    for f in fasta:
        query_len = len(f['seq'])
        if f['name'] not in blasta_lookup:
            f['blastp_status'] = 'Not Found'
            if len(f['seq']) < min_length:
                print('Warning: {0} <{1} Amino acids, skipping'.format(hit['query_name'], min_length))
            else:
                res.append(f) # Full sequence is uniquq

        else: # Was found by BLASTP
            blastp_status = 'Found, but no good hits'
            remaining_sequence = f['seq']

            for hit in blasta_lookup[f['name']]:
                if query_len == hit['match_len']: # 100 % match;
                    print('Warning: {0} 100% blast match, skipping'.format(hit['query_name']))
                    remaining_sequence = None
                    break
                if hit['pident'] > 95.0:
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
                # check it's not too many NNN's
                num_notN = len(remaining_sequence) - remaining_sequence.count('n')
                if num_notN < min_length:
                    print('Warning: {0} masked to <{1} Amino acids, skipping'.format(hit['query_name'], min_length))
                elif len(remaining_sequence) < min_length:
                    print('Warning: {0} <{1} Amino acids, skipping'.format(hit['query_name'], min_length))
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
                        squeeze_ns = [] # probably a 1-liner could do this...
                        last = None
                        for c in remaining_sequence:
                            if last == 'n':
                                last = c
                                continue
                            last = c
                            to_add = c
                            if to_add == 'n':
                                to_add = '-'
                            squeeze_ns.append(to_add)

                        print(squeeze_ns)
                        if squeeze_ns[0] == '-': squeeze_ns =squeeze_ns[1:]
                        oh_all_seqs.write(''.join([a for a in squeeze_ns if a not in ('_')]))
                        oh_all_seqs.write('\n')

    if res:
        resgl = genelist()
        resgl.load_list(res)
        print('Number of surviving peptides: {0}'.format(len(resgl)))
        resgl.sort('name')
        resgl.saveTSV('masked/masked_results-{0}.tsv'.format(stub), key_order=['name', 'blastp_status'])
        resgl.save('masked/masked_results-{0}.glb'.format(stub))
        # resgl.saveFASTA

oh_all_seqs.close()
