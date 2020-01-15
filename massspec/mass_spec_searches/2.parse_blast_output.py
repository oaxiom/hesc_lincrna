'''

For the masked peptides, see if we can find them in the MS data;

'''

import glob, sys, os
from glbase3 import *

#ms_data = glload('../mass_spec.glb')
#ms_data = glload('../mass_spec_pride.glb')
#ms_data = glload('../mass_spec_hipsci.glb')
ms_data = glload('../mass_spec_peptideatlas.glb')
print(ms_data)

for filename in glob.glob('../blast_searches/masked/*.glb'):

    stub = os.path.split(filename)[1].replace('.glb', '').split('-')[1]
    blasta = glload(filename)
    blasta = blasta.removeDuplicates('seq')

    res_fastas = {} #f['name']: 0 for f in blasta}
    res_peps = []
    p = progressbar(len(ms_data))
    for idx, peptide_fragment in enumerate(ms_data):
        #print(peptide_fragment)
        for f in blasta:
            if 'R' in f['seq'] or 'K' in f['seq'] or 'L' in f['seq']: # at least one possible pepetide cant be cut
                if f['name'] not in res_fastas:
                    res_fastas[f['name']] = 0
            else:
                continue # skip it; no possible peptide can be detected;

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
                    #'peptide_score': peptide_fragment['score']
                    # and match, contect, etc.
                    })

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
        num_matches = sum([res_fastas[i]>=1 for i in res_fastas])
        print('Number of matching genes with >=1 peptides: {0} ({1:.1f}%)'.format(num_matches, num_matches/len(resgl)*100))

        #for k in sorted(res_fastas):
        #    print('{0}:\t {1}'.format(k, res_fastas[k]))
        print()


