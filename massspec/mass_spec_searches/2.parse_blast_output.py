'''

For the masked peptides, see if we can find them in the MS data;

'''

import glob, sys, os
from glbase3 import *

ms_data = glload('../mass_spec.glb')

for filename in glob.glob('../blast_searches/masked/*.glb'):
    stub = os.path.split(filename)[1].replace('.glb', '').split('-')[1]
    blasta = glload(filename)

    # Parse the FASTA
    fasta = genelist('../fasta/fasta/{0}.fasta'.format(stub), format=format.fasta)
    fasta_lookup = {}

    res_peps = []
    res_fastas  = {} # for each fasta, a
    for peptide_fragment in ms_data:
        #print(peptide_fragment)
        for f in fasta:
            if peptide_fragment['seq'] in f['seq']:
                posl = peptide_fragment['seq'].find(f['seq'])
                posr = peptide_fragment['seq'].find(f['seq']) + len(peptide_fragment['seq'])

                res_peps.append({'peptide_fragment': peptide_fragment,
                    'context': f['seq'][max(posl-5, 0):posl].lower() + f['seq'][posl:posr] + f['seq'][posr:posr+5].lower(),
                    'pos': (posl, posr),
                    'peptide_score': peptide_fragment['score']
                    # and match, contect, etc.
                    })

                if f['name'] not in res_fastas:
                    res_fastas[f['name']] = 0
                res_fastas[f['name']] += 1

    resgl = genelist()
    resgl.load_list(res_peps)

    resgl.saveTSV('peptide_hits-{0}.tsv'.format(stub), key_order=['name', 'blastp_status'])
    resgl.save('peptide_hits-{0}.glb'.format(stub) )
    # resgl.saveFASTA

    print()
    print(stub)
    print('Number of matching peptides: {0}'.format(len(resgl)))
    for k in sorted(res_fastas):
        print('{0}:\t {1}'.format(k, res_fastas[k]))
    print()

