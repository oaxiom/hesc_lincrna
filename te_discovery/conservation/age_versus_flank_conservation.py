
'''

Measure the conservation of the lncRNAs versus the Family age of the TE and the flanks

'''
import sys
from glbase3 import glload, genelist, config
sys.path.append('../../')
import shared
config.draw_mode = 'pdf'

dfam = genelist('../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})

gl = glload('phyloP_conservation_table_per_TE_type.glb')
print(gl)

# fix the dfam annotations to get htem to match better
dfam_dict = {}
for te in dfam:
    dfam_dict[te['name']] = '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name'])

ages = glload('../dfam/te_ages.glb')
ages = {i['TE']: i['age'] for i in ages} # lookup;

def get_cols(TE_names):
    spot_cols = []
    for name in TE_names:
        if 'LINE' in name:
            spot_cols.append('blue')
        elif 'LTR' in name:
            spot_cols.append('red')
        elif 'SINE' in name:
            spot_cols.append('green')
        else:
            spot_cols.append('grey')
    return spot_cols

x = []
y = []
c = []
for te in gl:
    if te['TE'] not in dfam_dict:
        continue
    if dfam_dict[te['TE']] in ages:
        y.append(te['phyloP_tes'])
        x.append(ages[dfam_dict[te['TE']]])
        c.append(dfam_dict[te['TE']])

shared.nice_scatter(x=x, y=y,
    figsize=[2,2], spot_size=12,
    spot_cols=get_cols(c),
    label_t=0.25,
    filename='scat_age_vs_te.pdf',
    label=c,
    hlines=[0])

x = []
y = []
c = []
for te in gl:
    if te['TE'] not in dfam_dict:
        continue
    if dfam_dict[te['TE']] in ages:
        y.append(te['phyloP_nottes'])
        x.append(ages[dfam_dict[te['TE']]])
        c.append(dfam_dict[te['TE']])

shared.nice_scatter(x=x, y=y,
    figsize=[2,2], spot_size=12,
    spot_cols=get_cols(c),
    label_t=0.25,
    filename='scat_age_vs_notte.pdf',
    label=c, hlines=[0])
