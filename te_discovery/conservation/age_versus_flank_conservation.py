
'''

Measure the conservation of the lncRNAs versus the Family age of the TE and the flanks

'''

from glbase3 import glload

dfam = genelist('../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

gl = glload('phyloP_conservation_table_per_TE_type.glb')
print(gl)

# fix the dfam annotations to get htem to match better
dfam_dict = {}
for te in dfam:
    dfam_dict[te['name']] = '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name'])

ages = glload('../dfam/te_ages.glb')
print(ages)

