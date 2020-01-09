

from glbase3 import *
import gzip

# Some kind of problem in the middle of the file:
#form = {'force_tsv': True, 'seq': 0, 'score': 35}# , 'debug': 100000}
#ms = genelist('raw_data/hipsci.proteomics.maxquant.uniprot.MQ1_MQ2.20151023.peptides.txt.gz', format=form, gzip=True)

fa = genelist('raw_data/APD_Hs_all.fasta', format=format.fasta)
fa = fa.removeDuplicates('seq')

fa.save('mass_spec_peptideatlas.glb')

print(fa)
