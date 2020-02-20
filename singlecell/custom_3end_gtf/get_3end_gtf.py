
from glbase3 import *

gtf = glload(filename='../../transcript_assembly/packed/all_genes.glb')

ends = set([])

# Make a list of all the 'ensgs' that can be unambiguosly identified from their 3'ends:

lines = []
non_unique = 0
for idx, item in enumerate(gtf):
    if (idx+1) % 10000 == 0:
        print('Processed {:,}'.format(idx+1))
        #break

    if item['strand'] == '+':
        end = (item['loc']['chr'], item['loc']['right']-1000, item['loc']['right']+500, item['strand'])
    elif item['strand'] == '-':
        end = (item['loc']['chr'], item['loc']['left']-500, item['loc']['left']+1000, item['strand'])
    else:
        print('WARNING: no strand {0}'.format(item['name']))

    if end not in ends:
        # never seen this end, add and output it
        ends.add(end)
        if ';C;' in item['name']:
            typ = 'protein_coding'
        elif ';NC;' in item['name']:
            typ = 'lncRNA'

        tags = [
            'gene_id "{0}"'.format(item['transcript_id']),
            'transcript_id "{0}"'.format(item['transcript_id']),
            'transcript_name "{0}"'.format(item['name']),
            'gene_type "{0}"'.format(typ),
            'ensg "{0}"'.format(item['ensg']),
            'gene_name "{0}"'.format(item['name']),
            ]
        tags = "{0};".format('; '.join(tags))
        line = [item['loc']['chr'], '3ends', 'exon', str(end[1]), str(end[2]), '1000', item['strand'], '.', tags]
        lines.append({'loc': location(chr=end[0], left=end[1], right=end[2]), 'strand': item['strand'], 'line': '{0}\n'.format('\t'.join(line))})
    else:
        non_unique += 1

print('Exact end duplicates:')
print('Unique ends: {0}'.format(len(ends)))
print('Non-unique ends: {0}'.format(non_unique))

# Remove multiple ends:
all_ends = genelist()
all_ends.load_list(lines)
all_ends = all_ends.removeDuplicatesByLoc(mode='overlap', delta=0, use_strand=True)# keep the first entry you come across, delete_any_matches=True) # I then need to remove the ends that have multiple hits

oh = open('ends.gtf', 'w')
for item in all_ends:
    oh.write(item['line'])
oh.close()


#gtf = delayedlist(filename='../../transcript_assembly/gtf/current_gtf.gtf.gz', format=format.gtf, gzip=True)
