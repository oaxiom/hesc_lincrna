
import glob
from glbase3 import *
import gzip

# Some kind of problem in the middle of the file:
#form = {'force_tsv': True, 'seq': 0, 'score': 35}# , 'debug': 100000}
#ms = genelist('raw_data/hipsci.proteomics.maxquant.uniprot.MQ1_MQ2.20151023.peptides.txt.gz', format=form, gzip=True)

newl = []

peptides = {}

for filename in glob.glob('raw_data/PXD010557-hipsci/andromeda/*.res'):
    print(filename)
    oh = gzip.open('raw_data/human_peptides.txt.gz', 'rt')
    for line in oh:
        if '>' not in line:
            t = line.strip().split('\t')
            if t[0] not in peptides:
                peptides[t[0]] = 0
            peptides[t[0]] += 0
    oh.close()
    print('Found {0} peptides'.format(len(peptides)))

ms = genelist()
ms.load_list([{'seq': k, 'score': v} for k, v in peptides.items()])
print(ms)
print('Found {0} peptides'.format(len(ms)))
ms.save('mass_spec_hipsci.glb')
