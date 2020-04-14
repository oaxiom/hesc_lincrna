
'''

Measure the conservation of the lncRNAs versus the Family age of the TE and the flanks

'''

from glbase3 import glload, genelist
import shared

dfam = genelist('../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

gl = glload('phyloP_conservation_table_per_TE_type.glb')
print(gl)

# fix the dfam annotations to get htem to match better
dfam_dict = {}
for te in dfam:
    dfam_dict[te['name']] = '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name'])

ages = glload('../dfam/te_ages.glb')
ages = {i['TE']: i['age'] for i in ages} # lookup;

x = []
y = []
for te in gl:
    print(te['TE'])
    if te['TE'] not in dfam_dict:
        continue
    if dfam_dict[te['TE']] in ages:
        y.append(te['phyloP_tes'])
        x.append(ages[dfam_dict[te['TE']]])

shared.scat('scat_age_vs_te.pdf',
    x, y,
    'Age', 'TE phyloP',
    None, None,
    )

shared.contour('cont_age_vs_te.pdf',
    x, y,
    'Age', 'TE phyloP',
    [[0, 190], [-0.6, 0.6]],
    )

x = []
y = []
for te in gl:
    print(te['TE'])
    if te['TE'] not in dfam_dict:
        continue
    if dfam_dict[te['TE']] in ages:
        y.append(te['phyloP_nottes'])
        x.append(ages[dfam_dict[te['TE']]])

shared.scat('scat_age_vs_notte.pdf',
    x, y,
    'Age', 'not-TE phyloP',
    None, None,
    )

shared.contour('cont_age_vs_te.pdf',
    x, y,
    'Age', 'TE phyloP',
    [[0, 190], [-0.6, 0.6]],
    )
