
'''

Pack the FASTA into a glb for faster access

'''

import gzip
from glbase3 import genelist, format, glload

annotations = glload('../packed/all_genes.glb').getColumns(['transcript_id', 'enst', 'name', 'coding', 'expression', 'tags'])
print(annotations)
fasta = genelist('transcripts.fa.gz', format=format.fasta, gzip=True)
fasta = fasta.renameKey('name', 'transcript_id')

fasta = annotations.map(genelist=fasta, key='transcript_id')

fasta.save('transcripts.glb')
