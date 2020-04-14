
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

ages_gl = glload('../dfam/te_ages.glb')
ages = {}
for te in ages_gl:
    if '-int' in te['TE']:
        ages[te['TE'].replace('-int', '')] = te['age'] # Harmonise the names
    ages[te['TE']] = te['age'] # lookup;

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

keep_labels = set([
    #'DNA:TcMar-Tigger:Tigger1': res_type['DNA:TcMar-Tigger:Tigger1']['ncrna-all'],
    'SINE:Alu:AluSg',
    'SINE:Alu:AluSp',
    'SINE:Alu:AluY',
    'SINE:Alu:FRAM',

    'LINE:L2:L2',

    'LINE:L1:L1M2_orf2',
    'LINE:L1:L1M5_orf2',
    'LINE:L1:L1ME1_3end',
    'LINE:L1:L1HS_5end',
    'LINE:L1:L1M1_5end',
    'LINE:L1:L1M3_orf2',
    'LINE:L1:L1MC4_5end',
    'LINE:L1:L1P1_orf2',
    'LINE:L1:L1P4_orf2',
    'LINE:L1:L1PA3_3end',
    'LINE:L1:L1PA4_3end',

    'LTR:ERVL-MaLR:MLT1B',
    'LTR:ERV1:LTR7',
    'LTR:ERV1:LTR7B',
    'LTR:ERV1:LTR7Y',
    'LTR:ERV1:LTR8',
    'LTR:ERV1:LTR12C',
    'LTR:ERV1:LTR12E',
    'LTR:ERV1:HERVH',
    'LTR:ERV1:HERVFH21',
    'LTR:ERV1:HERVH48',
    'LTR:ERV1:HERV-Fc2',
    'LTR:ERV1:HERVS71',
    'LTR:ERV1:MER50-int',
    #'LTR:ERV1:HUERS-P3b', # No valid bundles
    #'LTR:ERV1:MER110', # Just 1 valid gene;

    'LTR:ERVK:HERVK',

    'LTR:ERVL-MaLR:MST-int',
    'LTR:ERVL-MaLR:THE1-int',
    #'', res_type['']['pc-all'],
    ])

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

spot_cols = get_cols(c)
labs = []
for c in c:
    if c in keep_labels:
        labs.append(c)
    else:
        labs.append('')

print(labs)

shared.nice_scatter(x=x, y=y,
    figsize=[2,2], spot_size=12,
    spot_cols=spot_cols,
    label_t=-20,
    filename='scat_age_vs_te.pdf',
    label=labs,
    doR=True,
    ylims=[-0.25, 0.6],
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

spot_cols = get_cols(c)
labs = []
for c in c:
    if c in keep_labels:
        labs.append(c)
    else:
        labs.append('')

print(labs)

shared.nice_scatter(x=x, y=y,
    figsize=[2,2], spot_size=12,
    spot_cols=spot_cols,
    label_t=-20,
    filename='scat_age_vs_notte.pdf',
    ylims=[-0.25, 0.6],
    doR=True,
    label=labs, hlines=[0])
