import os, sys
from glbase3 import genelist
import matplotlib.pyplot as plt

sys.path.append('../../')
import shared

gl = genelist('cell_cycle.csv', format={'grp': 17, 'phase': 20, 'S_score': 18, 'G2M_score': 19})

cluster = {}

t = 0.05

for item in gl:
    if item['grp'] not in cluster:
        cluster[item['grp']] = {'G1': 0, 'S': 0, 'G2M': 0}

    # redo the scoring; Scanpy is waaaay too generous
    phase = None
    if item['S_score'] < t and item['G2M_score'] < t:
        phase = 'G1'
    elif item['S_score'] > t and item['G2M_score'] >t:
        # In this case you just choose the stage chosen by Scanpy
        phase = item['phase']
    elif item['S_score'] > t and item['G2M_score'] < t:
        phase = 'S'
    elif item['S_score'] < t and item['G2M_score'] > t:
        phase = 'G2M'

    cluster[item['grp']][phase] += 1

for k in sorted(cluster):
    #print(k, cluster[k])
    print('\nCluster:', k)
    tot = sum(cluster[k].values())
    print('   G1  : {0:.1f}%'.format(cluster[k]['G1']/tot*100))
    print('   S   : {0:.1f}%'.format(cluster[k]['S']/tot*100))
    print('   G2/M: {0:.1f}%'.format(cluster[k]['G2M']/tot*100))

cluster = {c: cluster[c] for c in [4, 3, 2, 1, 0]}

shared.split_bar('cell_cycle.png', cluster, )
