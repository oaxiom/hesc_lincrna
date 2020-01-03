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
    p = progressbar(len(ms_data))
    for idx, peptide_fragment in enumerate(ms_data):
        #print(peptide_fragment)
        for f in fasta:
            if peptide_fragment['seq'] in f['seq']:
                posl = f['seq'].find(peptide_fragment['seq'])
                posr = posl + len(peptide_fragment['seq'])

                res_peps.append({
                    'name': f['name'],
                    'seq_length': len(f['seq']),
                    'seq': f['seq'],
                    'peptide_fragment': peptide_fragment['seq'],
                    'context': f['seq'][max(posl-5, 0):posl].lower() + f['seq'][posl:posr] + f['seq'][posr:posr+5].lower(),
                    'pos': (posl, posr),
                    'peptide_score': peptide_fragment['score']
                    # and match, contect, etc.
                    })

                if f['name'] not in res_fastas:
                    res_fastas[f['name']] = 0
                res_fastas[f['name']] += 1
        p.update(idx)

    if res_peps:
        resgl = genelist()
        resgl.load_list(res_peps)
        resgl.sort('name')
        resgl.saveTSV('peptide_hits-{0}.tsv'.format(stub), key_order=['name', 'seq_length', 'pos', 'peptide_fragment', 'context', 'peptide_score', 'seq' ])
        resgl.save('peptide_hits-{0}.glb'.format(stub) )
        # resgl.saveFASTA

        # per-transcript number of hits:
        resgl = genelist()
        resgl.load_list([{'name': k, 'num_peptides': res_fastas[k]} for k in res_fastas])
        resgl.sort('name')
        resgl.saveTSV('per_transcript_hits-{0}.tsv'.format(stub), key_order=['name', 'num_peptides'])
        #resgl.save('per_transcript_hits-{0}.glb'.format(stub) )
        print()
        print(stub)
        print('Number of matching peptides: {0}'.format(len(resgl)))
        for k in sorted(res_fastas):
            print('{0}:\t {1}'.format(k, res_fastas[k]))
        print()


