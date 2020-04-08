'''

For the masked peptides, see if we can find them in the MS data;

'''

import glob, sys, os, numpy, regex
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../')
import shared

res = {}

all_matches = glload('../results_gene.glb')

pep_hit = []

for pep in all_matches:
    if pep['insideTE'] != 'None':
        if 'LTR:ERVK:HERVK' in pep['insideTE']:
            pep_hit.append(pep['peptide_string'])

print('Found {0} HERVK peptides'.format(len(pep_hit)))
pep_hit = list(set(pep_hit))
print('Found {0} unique HERVK peptides'.format(len(pep_hit)))

print(pep_hit)

phoenix = genelist('phoenix.fa', format=format.fasta)

# See if simple matches will do it first.
finds = {seq['name']: [] for seq in phoenix}
finds['Not Found'] = []
for p in pep_hit:
    found = False
    for seq in phoenix:
        if p in seq['seq']:
            finds[seq['name']].append({'pep_seq': p, 'num_mismatch': 0})
            found = True
            continue
        # try n bp mismatch until pass:
        for num_mismatch in (1,2,3): # There are no useful matches after 3
            m = regex.findall("(%s){e<=%s}" % (p, num_mismatch), seq['seq'])
            if m: # In testing, only ever 1 hit;
                finds[seq['name']].append({'pep_seq': p, 'num_mismatch': num_mismatch})
                found = True
                break
    if not found: # Do another round looking for way wilder ones:
        for seq in phoenix:
            # try n bp mismatch until pass:
            for num_mismatch in (4,5,6,7, 8, 9, 10, 11, 12, 13, 14, 15, 16):
                m = regex.findall("(%s){s<=%s}" % (p, num_mismatch), seq['seq'])
                if m: # In testing, only ever 1 hit;
                    if found: # We have a find already:
                        if num_mismatch < found['num_mismatch']: # Only replace if a better match
                            found = {'typ': seq['name'], 'pep_seq': p, 'num_mismatch': num_mismatch}
                    else:
                        found = {'typ': seq['name'], 'pep_seq': p, 'num_mismatch': num_mismatch}
                    break

            #if found: # break again
            #    break

        # Add the best match:
        if found:
            typ = found['typ']
            del found['typ']
            finds[typ].append(found)
        else:
            finds['Not Found'].append({'pep_seq': p, 'num_mismatch': -1})

# Save a summary table:
oh = open('pep_matches.tsv', 'wt')
oh.write('{0}\n'.format('\t'.join(['HERVK', 'peptide', 'num_mismatch'])))
for typ in finds:
    for hit in finds[typ]:
        oh.write('{0}\n'.format('\t'.join([typ, hit['pep_seq'], str(hit['num_mismatch'])])))
oh.close()

print('Found {0}/{1} ({2:.1f}%) peptides in the HERVK seq'.format(len(sum(finds.values(), [])), len(pep_hit), len(sum(finds.values(), []))/len(pep_hit) * 100))
print(finds)
