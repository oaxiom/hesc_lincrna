

from glbase3 import *
import gzip

# Some kind of problem in the middle of the file:
#form = {'force_tsv': True, 'seq': 0, 'score': 35}# , 'debug': 100000}
#ms = genelist('raw_data/hipsci.proteomics.maxquant.uniprot.MQ1_MQ2.20151023.peptides.txt.gz', format=form, gzip=True)

newl = []
oh = gzip.open('raw_data/human_peptides.txt.gz', 'rt')
for line in oh:

    if 'Sequence' not in line:
        t = line.strip().split('\t')
        newl.append({'seq': t[0]})#, 'score': float(t[35])})

ms = genelist()
ms.load_list(newl)

print(ms)

ms.save('mass_spec_pride.glb')
