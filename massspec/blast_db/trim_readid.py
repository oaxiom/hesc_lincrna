
from glbase3 import *

fasta = genelist('../../transcript_assembly/get_CDS/gencode.v32.pc_translations.fa.gz', format=format.fasta, gzip=True)

for f in fasta:
    t = f['name'].split('|')
    f['name'] = '|'.join([t[0], t[6]])

fasta._optimiseData()
print(fasta)

fasta = fasta.removeDuplicates('seq')

fasta.saveFASTA('gencode.v32.pc_translations.fa', name='name')
