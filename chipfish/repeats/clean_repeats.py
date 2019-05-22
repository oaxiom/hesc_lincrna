

'''

clean the rmsk data to only include categories I am interested in:

'''

from glbase3 import *

rmsk_track_form = {"force_tsv": True, 'loc': 'location(chr=column[5], left=column[6], right=column[7])',
    'strand': 9, 'repName': 10, 'repClass': 11, 'repFamily': 12}

ll = delayedlist(filename='hg38_rmsk.tsv', format=rmsk_track_form)
ll_len = len(ll)

'''
These are the repClasses:

cat hg38_rmsk.tsv | cut -f 12 | sort | uniq

DNA
DNA?
LINE
LTR
LTR?
Low_complexity
RC
RC?
RNA
Retroposon
SINE
SINE?
Satellite
Simple_repeat
Unknown
rRNA
repClass
scRNA
snRNA
srpRNA
tRNA
'''

keep_classes = frozenset(['LINE', 'LTR', 'SINE', 'DNA', 'Retroposon'])

gsql = genome_sql(new=True, filename='hg38_rmsk_EREs.sql')
oh = open('hg38_rmsk_EREs.tsv', 'w')
keep = 0
p = progressbar(ll_len)
for idx, entry in enumerate(ll):
    if entry['repClass'] in keep_classes:
        oh.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (entry['loc']['chr'], entry['loc']['left'], entry['loc']['right'], entry['strand'], entry['repName'], entry['repClass'], entry['repFamily']))
        gsql.add_feature(entry['loc'], entry['loc'], 1, [1], [1], entry['repName'], entry['strand'], entry['repClass'])
        keep += 1
    p.update(idx)

print('Kept:', keep)
oh.close()
gsql.finalise()
