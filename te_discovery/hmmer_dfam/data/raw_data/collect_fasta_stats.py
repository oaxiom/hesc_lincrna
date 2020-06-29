'''
Summarise some FASSTA stats just to get things like number of base pairs, number of entries.
'''

import sys, os, glob

res = {}

for f in glob.glob('*/*.fa'):
    name = os.path.split(f)[1].replace('.fa', '') 
    print('Doing', name)
    oh = open(f, 'r')
    res[name] = {'num_fasta': 0, 'num_bps': 0}
    for line in oh:
        if line[0] == '>':
            res[name]['num_fasta'] += 1
        else:
            res[name]['num_bps'] += len(line.strip()) # Assume no other charachters;

    oh.close()

oh = open('summary.tsv', 'w')

oh.write('%s\n' % ('\t'.join(['name', 'num_fasta', 'num_bps'])))
for sam in sorted(res):
    oh.write('%s\n' % ('\t'.join([sam, str(res[sam]['num_fasta']), str(res[sam]['num_bps'])])))

oh.close()
        

