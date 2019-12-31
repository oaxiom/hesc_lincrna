'''
blastp -outfmt 6

 1.	 qseqid	 query (e.g., gene) sequence id
 2.	 sseqid	 subject (e.g., reference genome) sequence id
 3.	 pident	 percentage of identical matches
 4.	 length	 alignment length
 5.	 mismatch	 number of mismatches
 6.	 gapopen	 number of gap openings
 7.	 qstart	 start of alignment in query
 8.	 qend	 end of alignment in query
 9.	 sstart	 start of alignment in subject
 10.	 send	 end of alignment in subject
 11.	 evalue	 expect value
 12.	 bitscore	 bit score

   'qaccver saccver pident length mismatch gapopen qstart qend sstart send
   evalue bitscore', which is equivalent to the keyword 'std'
'''

# This script takes the blastp output, extracts the 100% ~ matches, and then deletes those hits in the
# FASTA peptide sequences,
# leaving a list of unique peptide sequences

import glob, sys, os
from glbase3 import *

form = {'query_name': 0,
    'hit_name': 1,
    'pident': 2,
    'match_len': 3,
    'qstart': 4,
    'qend': 5,
    'sstart': 6,
    'send': 7,
    'skiplines': -1,
    'force_tsv': True,
    }

for filename in glob.glob('blaster/table_*.tsv'):
    stub = os.path.split(filename)[1].replace('.tsv', '')
    blasta = genelist(filename, format=form)

    blasta_lookup = {}
    for b in blasta:
        if b['query_name'] not in blasta_lookup:
            blasta_lookup[b['query_name']] = []
        blasta_lookup[b['query_name']].append(b['query_name'])

    fasta = genelist('../fasta/fasta/{0}.fasta'.format(stub), format=format.fasta)
    fasta_lookup = {}
    for f in fasta:

        for b in blasta_lookup[f['name']]:
            query = fasta_lookup[hit['query_name']]
            query_len = len(query)

            if query_len == hit['match_len']: # 100 % match;
                print(hit['pident'])
                print('Warning: {0} 100% blast match, skipping'.format(hit['query_name']))
                print(hit)
                break

