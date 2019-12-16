
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

transcriptome_type = 'ncrna' # change me if you want to do ncrna or protein coding;

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # TODO: Add a new key here, PSC and notPSC
print(doms)
gencode = glload('../../te_transcripts/transcript_table_gencode_%s.glb' % transcriptome_type)
print(shared.convert_genocode_to_splice_sites(gencode[0]))
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
# rebuild to a lookup table:
dfam_lookup = {}
for item in dfam:
    dfam_lookup[item['name']] = '{0}:{1}:{2}'.format(item['type'], item['subtype'], item['name'])

#doms = doms.map(genelist=gencode_sliced, key='enst') # Only keep those with a known gene structure. Why?

def init_res_classes():
    return {'ES': 0, 'ES:': 0, 'ES-': 0, 'GENCODE': 0}

def init_res_store():
    return {
        'sp+': init_res_classes(),
        'sp-': init_res_classes(),
        }

cond_names = ['sp+', 'sp-']

# Store TEs by type:
res_tetype = {'LTR': init_res_store(),
    'SINE':init_res_store(),
    'LINE': init_res_store(),
    'Retroposon': init_res_store(),
    'DNA': init_res_store(),
    }
totals_tetype = {'LTR': init_res_classes(),
    'SINE': init_res_classes(),
    'LINE': init_res_classes(),
    'Retroposon': init_res_classes(),
    'DNA': init_res_classes()}

# Store TEs by family:
res_tefamily = defaultdict(init_res_store)
totals_tefamily = defaultdict(init_res_classes)

# Scores for normalisations
totals_all = init_res_classes()

# preprocss the doms list to remove non-coding genes;
newdoms = []
novel = []

# Here, later do ncrna, all and pc as a for:

data_to_process = {'ES+': doms, # Using all here, for now
    'ES:': ,
    'ES-': ,
    'GENCODE': gencode}

for dataset in data_to_process:
    print(dataset)
    for idx, gene in enumerate(data_to_process[dataset]):
        if (idx+1) % 10000 == 0:
            print('{:,}'.format(idx+1))

        _, tlength, splice_sites = shared.convert_genocode_to_splice_sites(gene)

        if transcriptome_type == 'pc': # You need to make sure you use the correct GENCODE background, otherwise the results will be odd:
            if ';NC;' in gene['name']: # i.e., remove non-coding
                continue
        elif transcriptome_type == 'ncrna':
            if ';C;' in gene['name']: # i.e., remove coding
                continue
        else:
            print('Type not found!')
            1/0

        for dom in gene['doms']: # The DOM locations are always in the + strand orientation.
            # test to see if the dom overlaps a splice site:
            for splice_site in splice_sites:
                full_name = dfam_lookup[dom['dom']] # get the full TE name from dfam:
                te_type = full_name.split(':')[0]

                if te_type not in totals_tetype: # Just the 5 I am interested in
                    continue

                if splice_site > dom['span'][0] and splice_site < dom['span'][1]:
                    #print(splice_site, dom)
                    # TE spans a splice site;
                    if dom['strand'] == '+': # same orientation as the transgene
                        res_tefamily[full_name]['sp+'][dataset] += 1
                        res_tetype[te_type]['sp+'][dataset] += 1
                    elif dom['strand'] == '-':
                        res_tefamily[full_name]['sp-'][dataset] += 1
                        res_tetype[te_type]['sp-'][dataset] += 1
                    else:
                        1/0

                    totals_tetype[te_type][dataset] += 1
                    totals_tefamily[full_name][dataset] += 1
                    totals_all[dataset] += 1

# Build them into glbase:
# Analyis by TYPE (LINE, SINE, LTR, etc)
for dataset in ['ES+', 'ES:', 'ES-']:
    newl = []
    for te_type in res_tetype.keys():
        enrich = []
        for part in cond_names:

            #obs = res_tetype[te_family][part][dataset] / (totals_tetype[te_family][dataset] + 1)
            #exp = res_tetype[te_family][part]['GENCODE'] / (totals_tetype[te_family]['GENCODE'] + 1)

            # This version is better, but again, it will normalise away the higher number of TEs/transcript discovered in our assemblies.
            # Still, it gives a second opinion.
            obs = (res_tetype[te_type][part][dataset]+1) / (totals_tetype[te_type][dataset])
            exp = ((res_tetype[te_type][part]['GENCODE']+1) / (totals_tetype[te_type]['GENCODE']))

            enrich.append(obs/exp)

            #print(dataset, part, te_type, obs, exp, enrich[-1], res_tetype[te_type][part][dataset]+1, totals_tetype[te_type][dataset], res_tetype[te_type][part]['GENCODE']+1, totals_tetype[te_type]['GENCODE'])

        newl.append({'name': te_type, 'conditions': enrich})
    expn = expression(loadable_list=newl, cond_names=cond_names)
    expn.saveTSV('freqs_{0}.tsv'.format(dataset))
    expn.heatmap('heatfreqs_{0}.png'.format(dataset), grid=True, heat_hei=0.013*len(expn), heat_wid=0.08, bracket=[0, 3], col_cluster=False)

# Analyiss by TE FAMILY:
for dataset in ['ES']:# , 'notES']:
    newl = []
    for te_family in res_tefamily.keys():
        enrich = []

        for part in cond_names:
            # trim trivial totals:
            if totals_tefamily[te_family][dataset] == 0:
                continue

            obs = (res_tefamily[te_family][part][dataset]) / (totals_tefamily[te_family][dataset])
            exp = ((res_tefamily[te_family][part]['GENCODE']+1) / (totals_tefamily[te_family]['GENCODE']+1))

            #print(te_family, res_tefamily[te_family][part][dataset], (totals_tefamily[te_family][dataset]+1), res_tefamily[te_family][part]['GENCODE'], (totals_tefamily[te_family]['GENCODE']+1))

            enr = obs/(exp)

            enrich.append(enr)

            #print('{0}\t{1:.2f} {2:.2f} {3:.2f}'.format(te_family, obs, exp, enrich[-1]))

        if enrich and max(enrich) > 3.0:
            #newl.append({'name': dfam_lookup[te_family], 'conditions': enrich})
            newl.append({'name': te_family, 'conditions': enrich})

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
                row_font_size=6, bracket=[0, 5],
                grid=True, heat_hei=0.0072*len(expn), heat_wid=0.08, col_cluster=False)
