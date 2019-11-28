
'''

Measure the density of TE insertions, and measure the frequencies at the TSS, TTS or in the middle of the transcript.

Output a heatmap with 3 columns;

'''

import glob, sys, os, gzip, numpy, math
from collections import defaultdict
import matplotlib.pyplot as plot
from glbase3 import glload, utils, expression, genelist, genome_sql

sys.path.append('../../../')
import shared
sys.path.append('../')
import shared_draw

draw = 'png'

type = 'ncrna'

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # Add a new key here, PSC and notPSC
gencode = glload('../../te_transcripts/transcript_table_gencode_%s.glb' % type)

dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
# rebuild to a lookup table:
dfam_lookup = {}
for item in dfam:
    dfam_lookup[item['name']] = '{0}:{1}:{2}'.format(item['type'], item['subtype'], item['name'])

print(dfam)

#doms = doms.map(genelist=gencode_sliced, key='enst') # Only keep those with a known gene structure. Why?

def init_res_store():
    return {
        'TSS': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'MID': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'TTS': {'ES': 0, 'notES': 0, 'GENCODE': 0}}

def init_res_store_totals():
    return {'ES': 0, 'notES': 0, 'GENCODE': 0}

res = {'LTR':
        {'TSS': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'MID': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'TTS': {'ES': 0, 'notES': 0, 'GENCODE': 0}},
    'SINE':
        {'TSS': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'MID': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'TTS': {'ES': 0, 'notES': 0, 'GENCODE': 0}},
    'LINE':
        {'TSS': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'MID': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'TTS': {'ES': 0, 'notES': 0, 'GENCODE': 0}},
    'Retroposon':
        {'TSS': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'MID': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'TTS': {'ES': 0, 'notES': 0, 'GENCODE': 0}},
    'DNA':
        {'TSS': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'MID': {'ES': 0, 'notES': 0, 'GENCODE': 0},
        'TTS': {'ES': 0, 'notES': 0, 'GENCODE': 0}},
    }

res_totals = {'LTR': init_res_store_totals(),
    'SINE': init_res_store_totals(),
    'LINE': init_res_store_totals(),
    'Retroposon': init_res_store_totals(),
    'DNA': init_res_store_totals()}

res_by_te = defaultdict(init_res_store)
res_by_te_totals = defaultdict(init_res_store_totals)
res_te_totals = init_res_store_totals()

# preprocss the doms list to remove non-coding genes;
newdoms = []
novel = []

data_to_process = {'ES': doms, # Using all here, for now
    #'notES': Not available yet,
    'GENCODE': gencode}

padding = 400

for dataset in data_to_process:
    print(dataset)
    for idx, gene in enumerate(data_to_process[dataset]):
        if (idx+1) % 10000 == 0:
            print('{:,}'.format(idx+1))

        # You need to make sure you use the correct GENCODE background, otherwise the results will be odd:

        if type == 'pc':
            if ';NC;' in gene['name']: # i.e., remove non-coding
                continue
        elif type == 'ncrna':
            if ';C;' in gene['name']: # i.e., remove coding
                continue
        else:
            print('Type not found!')
            1/0

        # LAter: Split into ES nor notES transcripts
        # Just count them all, at the highest family level:

        tlength = shared.get_transcript_length(gene)

        for dom in gene['doms']: # The DOM locations are always in the + strand orientation.
            #if 'HERVH' in dom['dom']:
            #    print(dom)
            full_name = dfam_lookup[dom['dom']] # get the full TE name from dfam:
            typ = full_name.split(':')[0]

            if typ not in res: # Just the 5 I am mainly interested in
                continue

            res_totals[typ][dataset] += 1
            res_by_te_totals[dom['dom']][dataset] += 1
            res_te_totals[dataset] += 1

            if dom['span'][0] < padding: # leave a small window for TSS inaccuracies
                res_by_te[dom['dom']]['TSS'][dataset] += 1
                res[typ]['TSS'][dataset] += 1

            elif dom['span'][1] > tlength-padding:
                res_by_te[dom['dom']]['TTS'][dataset] += 1
                res[typ]['TTS'][dataset] += 1

            else:
                res_by_te[dom['dom']]['MID'][dataset] += 1
                res[typ]['MID'][dataset] += 1

    # Build them into glbase:

for dataset in ['ES']:# , 'notES']:
    newl = []
    for typ in res_totals.keys():
        enrich = []
        for part in ('TSS', 'MID', 'TTS'):
            obs = res[typ][part][dataset] / (res_totals[typ][dataset] + 1)
            exp = res[typ][part]['GENCODE'] / (res_totals[typ]['GENCODE'] + 1)

            # This version is better, but again, it will normalise away the higher number of TEs/transcript discovered in our assemblies.
            # Still, it gives a second opinion.
            obs = (res_by_te[typ][part][dataset]+1) / (res_te_totals[dataset])
            exp = ((res_by_te[typ][part]['GENCODE']+1) / (res_te_totals['GENCODE']))

            print(typ, obs, exp)

            enrich.append(obs/exp)

        newl.append({'name': typ, 'conditions': enrich})
    expn = expression(loadable_list=newl, cond_names=['TSS', 'MID', 'TTS'])
    expn.saveTSV('freqs_{0}.tsv'.format(dataset))
    expn.heatmap('heatfreqs_{0}.png'.format(dataset), grid=True, heat_hei=0.013*len(expn), heat_wid=0.08, col_cluster=False)

for dataset in ['ES']:# , 'notES']:
    newl = []
    for typ in res_by_te_totals.keys():
        enrich = []
        #if res_by_te_totals[typ]['GENCODE'] < 50:
        #    continue
        #if res_by_te[typ]['MID']['GENCODE'] < 10:
        #    continue
        for part in ('TSS', 'MID', 'TTS'):
            # trim trivial totals:
            '''
            # This is wrong as it would normalise away higher frequency of detection of TEs
            obs = (res_by_te[typ][part][dataset]+1) / (res_by_te_totals[typ][dataset]+1)
            exp = ((res_by_te[typ][part]['GENCODE']+1) / (res_by_te_totals[typ]['GENCODE']+1))
            '''
            # This version is better, but again, it will normalise away the higher number of TEs/transcript discovered in our assemblies.
            # Still, it gives a second opinion.
            obs = (res_by_te[typ][part][dataset]+1) / (res_te_totals[dataset]/1000)
            exp = ((res_by_te[typ][part]['GENCODE']+1) / (res_te_totals['GENCODE']/1000))

            print(part, typ, obs, exp, res_by_te[typ][part][dataset], (res_by_te_totals[typ][dataset]+1), res_by_te[typ][part]['GENCODE'], (res_by_te_totals[typ]['GENCODE']+1))

            enr = obs/exp

            enrich.append(enr)

            #print('{0} \t{1:.2f} {2:.2f} {3:.2f}'.format(typ, obs, exp, enrich[-1]))

        if max(enrich) > 3.0:
            newl.append({'name': dfam_lookup[typ], 'conditions': enrich})

    # split them up by class:
    for cl in ['SINE', 'LINE', 'LTR', 'Retroposon', 'DNA']:
        newl2 = []
        for item in newl:
            if cl in item['name']:
                newl2.append(item)

        expn = expression(loadable_list=newl2, cond_names=['TSS', 'MID', 'TTS'])

        if expn:
            expn.saveTSV('freqs_by_class_{0}-{1}.tsv'.format(cl, dataset))
            expn.heatmap('heatfreqs_by_class_{0}-{1}.png'.format(cl, dataset), size=[6,12],
                row_font_size=6, bracket=[0, 3],
                grid=True, heat_hei=0.0072*len(expn), heat_wid=0.08, col_cluster=False)
