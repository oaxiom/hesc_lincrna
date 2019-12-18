
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

doms = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # TODO: Add a new key here, PSC and notPSC
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
# rebuild to a lookup table:
dfam_lookup = {item['name']: '{0}:{1}:{2}'.format(item['type'], item['subtype'], item['name']) for item in dfam}

delta = 1000 # bp around the splice site to look for TEs

def draw_pileup(filename,
        enriched_plus, enriched_neg,
        unbiased_plus, unbiased_neg,
        depleted_plus, depleted_neg,
        gencode_plus, gencode_neg):

    print('Drawing {0}'.format(filename))
    fig = plot.figure()

    xs = numpy.arange(0, delta*2)

    datap = {
        'GENCODE     +': gencode_plus,
        'ES-enriched +': enriched_plus,
        'ES-unbiased +': unbiased_plus,
        'ES-depleted +': depleted_plus
        }

    datam = {
        'GENCODE     -': gencode_neg,
        'ES-enriched +': enriched_neg,
        'ES-unbiased +': unbiased_neg,
        'ES-depleted +': depleted_neg
        }

    ax1 = fig.add_subplot(211)
    ymax = -1e6
    for k in datap:
        if datap[k] is not None:
            ax1.plot(xs, datap[k], label=k) # add legend
            if datap[k].max() > ymax: ymax = datap[k].max()

    ax2 = fig.add_subplot(212)
    for k in datam:
        if datam[k] is not None:
            ax2.plot(xs, -datam[k], label=k) # add legend
            if datam[k].max() > ymax: ymax = datam[k].max()

    ax1.legend()
    ax2.legend()

    ax1.axvline(delta, ls=":", color="grey")
    ax2.axvline(delta, ls=":", color="grey")

    ax1.set_xticklabels('')

    ax1.set_ylim([0, ymax])
    ax2.set_ylim([-ymax, 0])

    #ax2.set_xlim([0, delta*2])
    ax2.set_xticks([0, delta//2, delta, (delta//2)*3, delta*2])
    ax2.set_xticklabels([-delta, -delta//2, 0, delta//2, delta])
    fig.savefig(filename)
    plot.close(fig)


def init_res_classes_count():
    return {'ES+': 0, 'ES:': 0, 'ES-':0, 'GENCODE': 0}
def init_res_classes():
    return {'ES+': numpy.zeros(2000), 'ES:': numpy.zeros(2000), 'ES-': numpy.zeros(2000), 'GENCODE': numpy.zeros(2000)}
def init_res_store():
    return {'sp+': init_res_classes(), 'sp-': init_res_classes()}
def init_res_totals_classes():
    return {'ES+': 0, 'ES:': 0, 'ES-':0, 'GENCODE': 0}

for transcriptome_type in ('ncrna', 'pc'): # change me if you want to do ncrna or protein coding;
    gencode = glload('../../te_transcripts/transcript_table_gencode_%s.glb' % transcriptome_type) # needs to be loaded based on the type

    cond_names = ['sp+', 'sp-']

    # Store TEs by type:
    res_tetype = {'LTR': init_res_store(),
        'SINE':init_res_store(),
        'LINE': init_res_store(),
        'Retroposon': init_res_store(),
        'DNA': init_res_store(),
        }
    totals_tetype = {'LTR': init_res_totals_classes(),
        'SINE': init_res_totals_classes(),
        'LINE': init_res_totals_classes(),
        'Retroposon': init_res_totals_classes(),
        'DNA': init_res_totals_classes()}

    # Store TEs by family:
    res_tefamily = defaultdict(init_res_store)
    totals_tefamily = defaultdict(init_res_totals_classes)
    totals_tefamily_counts = defaultdict(init_res_classes_count)

    # Scores for normalisations
    totals_all = init_res_totals_classes()

    # preprocss the doms list to remove non-coding genes;
    newdoms = []
    novel = []

    # Here, later do ncrna, all and pc as a for:

    data_to_process = {'ES+': doms.get(key='expression', value='enriched'), # Using all here, for now
        'ES:': doms.get(key='expression', value='unbiased'),
        'ES-': doms.get(key='expression', value='depleted'),
        'GENCODE': gencode}

    data_lens = {k: len(data_to_process[k]) for k in data_to_process}

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

            for ss in splice_sites:
                for dom in gene['doms']:
                    full_name = dfam_lookup[dom['dom']]
                    te_type = full_name.split(':')[0]
                    if te_type not in res_tetype.keys():
                        continue

                    # get the position relative to the ss:
                    ss_rel_left = dom['span'][0] - ss + delta
                    ss_rel_right = dom['span'][1] - ss + delta

                    if ss_rel_right < 0:
                        continue
                    if ss_rel_left > delta *2:
                        continue

                    ss_rel_left = ss_rel_left if ss_rel_left > 0 else 0
                    ss_rel_right = ss_rel_right if ss_rel_right < delta*2 else delta*2

                    #print(ss, dom['span'], ss_rel_left, ss_rel_right)
                    #print(res_tefamily[full_name]['sp+'][dataset])

                    if dom['strand'] == '+': # same orientation as the transgene
                        res_tefamily[full_name]['sp+'][dataset][ss_rel_left:ss_rel_right] += 1
                        res_tetype[te_type]['sp+'][dataset][ss_rel_left:ss_rel_right] += 1
                    elif dom['strand'] == '-':
                        res_tefamily[full_name]['sp-'][dataset][ss_rel_left:ss_rel_right] += 1
                        res_tetype[te_type]['sp-'][dataset][ss_rel_left:ss_rel_right] += 1

                    totals_tetype[te_type][dataset] += 1
                    totals_tefamily[full_name][dataset] += 1
                    totals_all[dataset] += 1
                    totals_tefamily_counts[full_name][dataset] += 1

    # draw the pileups
    for te_type in res_tetype.keys():
        # normalise to the total number of TEs? Or to the number of transcripts?
        '''
        draw_pileup('pile-type-{0}-{1}.png'.format(te_type, transcriptome_type),
            res_tetype[te_type]['sp+']['ES+'] / totals_all['ES+'],
            res_tetype[te_type]['sp-']['ES+'] / totals_all['ES+'],
            res_tetype[te_type]['sp+']['ES:'] / totals_all['ES:'],
            res_tetype[te_type]['sp-']['ES:'] / totals_all['ES:'],
            res_tetype[te_type]['sp+']['ES-'] / totals_all['ES-'],
            res_tetype[te_type]['sp-']['ES-'] / totals_all['ES-'],
            res_tetype[te_type]['sp+']['GENCODE'] / totals_all['GENCODE'],
            res_tetype[te_type]['sp-']['GENCODE'] / totals_all['GENCODE'],
            )
        '''
        draw_pileup('pile-type-{0}-{1}.png'.format(te_type, transcriptome_type),
            res_tetype[te_type]['sp+']['ES+'] / data_lens['ES+'],
            res_tetype[te_type]['sp-']['ES+'] / data_lens['ES+'],
            res_tetype[te_type]['sp+']['ES:'] / data_lens['ES:'],
            res_tetype[te_type]['sp-']['ES:'] / data_lens['ES:'],
            res_tetype[te_type]['sp+']['ES-'] / data_lens['ES-'],
            res_tetype[te_type]['sp-']['ES-'] / data_lens['ES-'],
            res_tetype[te_type]['sp+']['GENCODE'] / data_lens['GENCODE'],
            res_tetype[te_type]['sp-']['GENCODE'] / data_lens['GENCODE'],
            )

    # draw the pileups
    for te_family in res_tefamily.keys():
        # normalise to the total number of TEs? Or to the number of transcripts?

        if max(totals_tefamily_counts[te_family].values()) < 100:
            continue
        '''
        draw_pileup('pileups/family-{0}-{1}.png'.format(te_family, transcriptome_type),
            res_tefamily[te_family]['sp+']['ES+'] / totals_all['ES+'],
            res_tefamily[te_family]['sp-']['ES+'] / totals_all['ES+'],
            res_tefamily[te_family]['sp+']['ES:'] / totals_all['ES:'],
            res_tefamily[te_family]['sp-']['ES:'] / totals_all['ES:'],
            res_tefamily[te_family]['sp+']['ES-'] / totals_all['ES-'],
            res_tefamily[te_family]['sp-']['ES-'] / totals_all['ES-'],
            res_tefamily[te_family]['sp+']['GENCODE'] / totals_all['GENCODE'],
            res_tefamily[te_family]['sp-']['GENCODE'] / totals_all['GENCODE'],
            )
        '''
        draw_pileup('pileups/family-{0}-{1}.png'.format(te_family, transcriptome_type),
            res_tefamily[te_family]['sp+']['ES+'] / data_lens['ES+'],
            res_tefamily[te_family]['sp-']['ES+'] / data_lens['ES+'],
            res_tefamily[te_family]['sp+']['ES:'] / data_lens['ES:'],
            res_tefamily[te_family]['sp-']['ES:'] / data_lens['ES:'],
            res_tefamily[te_family]['sp+']['ES-'] / data_lens['ES-'],
            res_tefamily[te_family]['sp-']['ES-'] / data_lens['ES-'],
            res_tefamily[te_family]['sp+']['GENCODE'] / data_lens['GENCODE'],
            res_tefamily[te_family]['sp-']['GENCODE'] / data_lens['GENCODE'],
            )
