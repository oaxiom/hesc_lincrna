import logging, matplotlib, os, sys, glob, numpy, gzip, pickle, math
from glbase3 import genelist, glload

transcript_data = glload('../../../../transcript_assembly/packed/all_genes.glb')
transcript_data = {i['transcript_id']: i for i in transcript_data}

transcript_names = []
oh = gzip.open('gene_names.filtered.tsv.gz', 'rt')
for line in oh:
    transcript_names.append(line.strip())
oh.close()

oh = open('../../dense_array.tsv', 'rt')

res = {
    'm': {'coding': [], 'noncoding': []},
    's': {'coding': [], 'noncoding': []},
    'cv': {'coding': [], 'noncoding': []},
    'num_cells': {'coding': [], 'noncoding': []},
    }


min_counts = 2

# Save out a dense array of the sparse array:
__skipped = 0
for idx, line in enumerate(oh):
    line = line.strip().split()
    #print(line)
    # Value here is ln(counts+1)
    e = [float(i) for i in line if float(i) >= min_counts] # this is in counts ...

    if not e:
        __skipped += 1
        continue
    tdata = transcript_data[transcript_names[idx]]

    #print(tdata['name'], len(e), numpy.min(e), numpy.max(e), numpy.mean(e), numpy.std(e))

    if tdata['coding'] != 'NA': # There are 9 weird transcripts that FELLnc didn't place;
        res['m'][tdata['coding']].append(numpy.mean(e))
        res['s'][tdata['coding']].append(numpy.std(e))
        res['cv'][tdata['coding']].append(numpy.mean(e) / numpy.std(e))
        res['num_cells'][tdata['coding']].append((len(e) / 30001) * 100)

    if (idx+1) % 1000 == 0:
        print('{0}'.format(idx+1))

print('Skipped {:,} transcripts with counts < {}'.format(__skipped, min_counts-1))

oh.close()

oh = open('res.pickle', 'wb')
pickle.dump(res, oh)
oh.close()

