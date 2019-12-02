
'''

Measure the density of TE insertions at splice sites

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

type = 'ncrna' # change me if you want to do ncrna or protein coding;

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # TODO: Add a new key here, PSC and notPSC
gencode = glload('../../te_transcripts/transcript_table_gencode_%s.glb' % type)
print(shared.convert_genocode_to_splice_sites(gencode[0]))
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
# rebuild to a lookup table:
dfam_lookup = {}
for item in dfam:
    dfam_lookup[item['name']] = '{0}:{1}:{2}'.format(item['type'], item['subtype'], item['name'])

#doms = doms.map(genelist=gencode_sliced, key='enst') # Only keep those with a known gene structure. Why?

def init_res_classes():
    return {'ES': 0, 'notES': 0, 'GENCODE': 0}

def init_res_store():
    return {
        'sp+': init_res_classes(),
        'sp-': init_res_classes(),
        }

cond_names = ['sp+', 'sp-']

res = {'LTR': init_res_store(),
    'SINE':init_res_store(),
    'LINE': init_res_store(),
    'Retroposon': init_res_store(),
    'DNA': init_res_store(),
    }

res_totals = {'LTR': init_res_classes(),
    'SINE': init_res_classes(),
    'LINE': init_res_classes(),
    'Retroposon': init_res_classes(),
    'DNA': init_res_classes()}

res_by_te = defaultdict(init_res_store)
res_by_te_totals = defaultdict(init_res_classes)
res_te_totals = init_res_classes()

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

        _, tlength, splice_sites = shared.convert_genocode_to_splice_sites(gene)

        if type == 'pc':
            if ';NC;' in gene['name']: # i.e., remove non-coding
                continue
        elif type == 'ncrna':
            if ';C;' in gene['name']: # i.e., remove coding
                continue
        else:
            print('Type not found!')
            1/0

        for dom in gene['doms']: # The DOM locations are always in the + strand orientation.
            # test ot see if the dom overlaps a splice site:
            for splice_site in splice_sites:
                full_name = dfam_lookup[dom['dom']] # get the full TE name from dfam:
                typ = full_name.split(':')[0]

                if typ not in res: # Just the 5 I am mainly interested in
                    continue

                if splice_site > dom['span'][0] and splice_site < dom['span'][1]:
                    # TE spans a splice site;
                    if dom['strand'] == '+': # same orientation as the transgene
                        res_by_te[typ]['sp+'][dataset] += 1
                    elif dom['strand'] == '-':
                        res_by_te[typ]['sp-'][dataset] += 1

                    res_totals[typ][dataset] += 1
                    res_by_te_totals[dom['dom']][dataset] += 1
                    res_te_totals[dataset] += 1

print(res_by_te)

# Build them into glbase:
for dataset in ['ES']:# , 'notES']:
    newl = []
    for typ in res_totals.keys():
        enrich = []
        for part in ('sp+', 'sp-'):
            obs = res[typ][part][dataset] / (res_totals[typ][dataset] + 1)
            exp = res[typ][part]['GENCODE'] / (res_totals[typ]['GENCODE'] + 1)

            # This version is better, but again, it will normalise away the higher number of TEs/transcript discovered in our assemblies.
            # Still, it gives a second opinion.
            obs = (res_by_te[typ][part][dataset]+1) / (res_te_totals[dataset])
            exp = ((res_by_te[typ][part]['GENCODE']+1) / (res_te_totals['GENCODE']))

            print(typ, obs, exp)

            enrich.append(obs/exp)

        newl.append({'name': typ, 'conditions': enrich})
    expn = expression(loadable_list=newl, cond_names=cond_names)
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
        for part in cond_names:
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

        expn = expression(loadable_list=newl2, cond_names=cond_names)

        if expn:
            expn.saveTSV('freqs_by_class_{0}-{1}.tsv'.format(cl, dataset))
            expn.heatmap('heatfreqs_by_class_{0}-{1}.png'.format(cl, dataset), size=[6,12],
                row_font_size=6, bracket=[0, 3],
                grid=True, heat_hei=0.0072*len(expn), heat_wid=0.08, col_cluster=False)
