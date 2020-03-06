
import pickle, numpy
import matplotlib
import matplotlib.pyplot as plot
from glbase3 import *
from collections import defaultdict
matplotlib.rc('font', serif='Helvetica')

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

def boxplots(filename, data, qs):
    fig = plot.figure(figsize=[2.8,4])
    mmheat_hei = 0.1+(0.022*len(data))
    fig.subplots_adjust(left=0.4, right=0.8, top=mmheat_hei, bottom=0.1)
    ax = fig.add_subplot(111)
    ax.tick_params(right=True)
    m = numpy.median(data['no TE'])

    ax.axvline(m, ls=":", lw=0.5, color="grey") # add a grey line at zero for better orientation
    ax.axvline(m-1, ls=":", lw=0.5, color="grey") # add a grey line at zero for better orientation
    ax.axvline(m+1, ls=":", lw=0.5, color="grey") # add a grey line at zero for better orientation

    dats = numpy.array(list(data.values()))
    r = ax.boxplot(dats,
        showfliers=False,
        whis=True,
        patch_artist=True,
        widths=0.5, vert=False)

    plot.setp(r['medians'], color='black', lw=2) # set nicer colours
    plot.setp(r['boxes'], color='black', lw=0.5)
    plot.setp(r['caps'], color="grey", lw=0.5)
    plot.setp(r['whiskers'], color="grey", lw=0.5)

    ax.set_yticklabels(data.keys())

    gtm = '#FF8A87'
    ltm = '#92A7FF'

    '''
    cols = []
    for k in data:
        if numpy.median(data[k]) > m:
            cols.append(gtm)
        else:
            cols.append(ltm)
        print(k, cols[-1])
    cols[-1] = 'grey'
    #cols.reverse()
    print(len(r['boxes']), len(cols))
    print(r['boxes'])
    print(cols)
    for patch, color in zip(r['boxes'], cols):
        patch.set_facecolor(color)
    '''

    for i, q, k, p in zip(range(0, len(data)), qs, data, r['boxes']):
        #ax.text(6.3, i+0.5, q, ha='left', va='center', fontsize=6,)
        if qs[q] < 0.05:
            ax.text(6.5, i+1, '*{0:.1e}'.format(qs[q]), ha='left', va='center', fontsize=6,)
            if numpy.median(data[k]) > m:
                p.set_facecolor(gtm)
            else:
                p.set_facecolor(ltm)
        else:
            ax.text(6.5, i+1, '{0:.1e}'.format(qs[q]), ha='left', va='center', fontsize=6,)
            p.set_facecolor('lightgrey')

        if k == 'no TE':
            p.set_facecolor('grey')

    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]

    fig.savefig(filename)

# The PC TEs:
qvals = genelist('by_transcript_class/tab_pc-all.tsv', format={'force_tsv': True, 'name': 0, 'q': 2})

# figure is bottom to top:
data = {
    'DNA:TcMar-Tigger:Tigger1': res_type['DNA:TcMar-Tigger:Tigger1']['pc-all'],
    'SINE:MIR:MIR': res_type['SINE:MIR:MIR']['pc-all'],
    'SINE:Alu:AluSp': res_type['SINE:Alu:AluSp']['pc-all'],
    'LINE:L2:L2': res_type['LINE:L2:L2']['pc-all'],
    'LINE:L1:L1M2_orf2': res_type['LINE:L1:L1M2_orf2']['pc-all'],
    'LINE:L1:L1M5_orf2': res_type['LINE:L1:L1M5_orf2']['pc-all'],
    'LINE:L1:L1HS_5end': res_type['LINE:L1:L1HS_5end']['pc-all'],
    'LTR:ERVL-MaLR:MLT1B': res_type['LTR:ERVL-MaLR:MLT1B']['pc-all'],
    'LTR:ERV1:LTR7': res_type['LTR:ERV1:LTR7']['pc-all'],
    'LTR:ERV1:LTR7Y': res_type['LTR:ERV1:LTR7Y']['pc-all'],
    'LTR:ERV1:HERVH': res_type['LTR:ERV1:HERVH']['pc-all'],
    #'', res_type['']['pc-all'],
    'no TE': res['pc-all']['nonTE'],
    }
qs = {}
for k in data:
    if k == 'no TE':
        qs['no TE'] = 0.0
        continue
    qs[k] = qvals.get(key='name', value=k, mode='lazy')[0]['q']
boxplots('pc.pdf', data, qs)

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
boxplots('ncrna.pdf', data, qs)
