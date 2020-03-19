
import pickle, numpy, sys
import matplotlib
import matplotlib.pyplot as plot
from glbase3 import *
from collections import defaultdict
matplotlib.rc('font', serif='Helvetica')

sys.path.append('../../')
import shared

def class_dict():
    return {
        'pc-all': {'TE': [], 'nonTE': []},
        'ncrna-all': {'TE': [], 'nonTE': []},

        'pc-known': {'TE': [], 'nonTE': []},
        'pc-variant': {'TE': [], 'nonTE': []},

        'ncrna-known': {'TE': [], 'nonTE': []},
        'ncrna-variant': {'TE': [], 'nonTE': []},
        'ncrna-unknown': {'TE': [], 'nonTE': []},
        }

def dict_builder():
    return {k: [] for k in class_dict()}

# Load pickles:
oh = open('by_transcript_class/res.pickle', 'rb')
res = pickle.load(oh)
oh.close()
oh = open('by_transcript_class/res_type.pickle', 'rb')
res_type = pickle.load(oh)
oh.close()

# selected classes;

# The PC TEs:
qvals = genelist('by_transcript_class/tab_pc-all.tsv', format={'force_tsv': True, 'name': 0, 'q': 2})

# figure is bottom to top:
data = [
    'DNA:TcMar-Tigger:Tigger1',
    'SINE:MIR:MIR',
    'SINE:Alu:AluSp',
    'SINE:Alu:AluJb',
    'SINE:Alu:FRAM',
    'LINE:L2:L2',
    'LINE:L1:L1M2_orf2',
    'LINE:L1:L1M5_orf2',
    'LINE:L1:L1HS_5end',
    'LTR:ERVL-MaLR:MLT1B',
    'LTR:ERV1:LTR7',
    'LTR:ERV1:LTR7Y',
    'LTR:ERV1:HERVH',
    'LTR:ERV1:HERV-Fc2',
    'LTR:ERV1:HERVE',
    #'', res_type['']['pc-all'],
    ]

data = {te: res_type[te]['pc-all'] for te in data}
data['no TE'] = res['pc-all']['nonTE']

qs = {}
for k in data:
    if k == 'no TE':
        qs['no TE'] = 0.0
        continue
    qs[k] = qvals.get(key='name', value=k, mode='lazy')[0]['q']
shared.boxplots('pc.pdf', data, qs)

# The ncRNA TEs:
qvals = genelist('by_transcript_class/tab_ncrna-all.tsv', format={'force_tsv': True, 'name': 0, 'q': 2})

# figure is bottom to top:
data = [
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
    'LTR:ERV1:HUERS-P3b',
    'LTR:ERV1:MER110',

    'LTR:ERVK:HERVK',

    'LTR:ERVL-MaLR:MST-int',
    'LTR:ERVL-MaLR:THE1-int',
    #'', res_type['']['pc-all'],

    ]

data = {te: res_type[te]['ncrna-all'] for te in data}
data['no TE'] = res['ncrna-all']['nonTE']

qs = {}
for k in data:
    if k == 'no TE':
        qs['no TE'] = 0.0
        continue
    qs[k] = qvals.get(key='name', value=k, mode='lazy')[0]['q']
shared.boxplots('ncrna.pdf', data, qs)
