
'''

SOLO TEs are these transcripts that contain intact, or semi intact unspliced transcripts.

As we don't trust the short read data to assemble these, we only consider them from the pacbio data:

'''

import sys
from glbase3 import glload, genelist

sys.path.append('../simple_summaries/')
import pies

all_te_transcripts = glload('../te_transcripts/transcript_table_merged.mapped.glb')
dfam = genelist('../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

solotes = []

type_subtype_counts = {}

stats = {'known_unkonwn': {'known': 0, 'novel': 0, 'unknown': 0},
    'coding_noncoding': {'coding': 0, 'noncoding': 0}
    }

for trans in all_te_transcripts:
    if 'LR' not in trans['name']:
        continue
    if trans['exonCounts'] > 1:
        continue
    # add a typ_subtype key:
    ts = set([])
    for d in trans['doms']:
        te = dfam.get(key='name', value=d['dom'])[0]
        type_subtyp = '%s:%s' % (te['type'], te['subtype'])
        ts.add(type_subtyp)
        if type_subtyp not in type_subtype_counts:
            type_subtype_counts[type_subtyp] = 0
        type_subtype_counts[type_subtyp] += 1

    trans['te_type'] = '; '.join(ts)
    solotes.append(trans)

    # collect stats;
    if ';!' in trans['name']:
        stats['known_unkonwn']['unknown'] += 1
    elif ';=' in trans['name']:
        stats['known_unkonwn']['known'] += 1
    elif ';~' in trans['name']:
        stats['known_unkonwn']['novel'] += 1

    if ';C;' in trans['name']:
        stats['coding_noncoding']['coding'] += 1
    elif ';NC;' in trans['name']:
        stats['coding_noncoding']['noncoding'] += 1

newgl = genelist()
newgl.load_list(solotes)
newgl.save('solo_tes.glb')
newgl.saveTSV('solo_tes.tsv', key_order=['ensg', 'enst', 'te_type',])

# collect some stats and pies;
for k in stats:
    pies.pie('pie_%s.png' % k, stats[k].values(), stats[k].keys(), k)
