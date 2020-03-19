import sys, os, numpy, math
from scipy.stats import mannwhitneyu, wilcoxon, ttest_ind, ttest_1samp
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared
sys.path.append('../')
import shared_bundle

# These have the official GENCODE CDS, and the predicted (about ~80% accurate)
all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
tes = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
dfam_dict = {te['name']: '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name']) for te in dfam}

# Broad summary:
def process_bundles(bundle):
    res_fcs = {}
    ps = {}
    gene_with_noTE_and_TE_transcript = 0
    has_no_with_te_transcript = 0
    has_no_nonte_transcript = 0

    # FOR P calc:
    tpms_withTE = {}
    tpms_noTE = {}

    for gene in bundle:
        tpms_for_no_te = []
        tpms_for_with_te = {}
        for transcript in bundle[gene]:
            if transcript['TEs']:
                unq_tes = set([t['dom'] for t in transcript['doms']])
                for te in unq_tes:
                    full_name = dfam_dict[te]
                    if full_name not in tpms_for_with_te:
                        tpms_for_with_te[full_name] = []
                    tpms_for_with_te[full_name].append(transcript['TPM'])
            else: # No TE:
                tpms_for_no_te.append(transcript['TPM'])

        # Get FC:
        # A few ways to do this, take the mean or the max
        if not tpms_for_with_te:
            has_no_with_te_transcript += 1
            continue # No paired
        if not tpms_for_no_te: # There is a few!
            has_no_nonte_transcript += 1
            continue

        gene_with_noTE_and_TE_transcript += 1
        for te in tpms_for_with_te:
            fc = utils.fold_change(max(tpms_for_no_te), max(tpms_for_with_te[te]), pad=0.1) # correct way around
            if 'LTR7' in te and '7' in te[-1]:
                print(te, tpms_for_no_te, tpms_for_with_te[te], fc)
            #fc = utils.fold_change(numpy.mean(tpms_for_no_te), numpy.mean(tpms_for_with_te[te]), pad=0.01)

            # You need to think about this slightly odd way of generating a P value, but it basically keeps all genes in each category
            # that are with or without a specific TE, and then does a MWU against that;
            if te not in tpms_noTE:
                tpms_noTE[te] = []
            if te not in tpms_withTE:
                tpms_withTE[te] = []

            tpms_noTE[te] += tpms_for_no_te
            tpms_withTE[te] += tpms_for_with_te[te]

            if te not in res_fcs:
                res_fcs[te] = []
            res_fcs[te].append(fc)

    # Figure out the P:
    ps = {}
    for te in tpms_withTE:
        #ps[te] = mannwhitneyu(tpms_noTE[te], tpms_withTE[te], alternative='two-sided')[1]
        #ps[te] = ttest_ind(
        #    [math.log2(v) for v in tpms_noTE[te]],
        #    [math.log2(v) for v in tpms_withTE[te]],
        #    equal_var=True)[1]
        #ps[te] = ttest_ind(
        #    tpms_noTE[te],
        #    tpms_withTE[te],
        #    equal_var=False)[1]
        #print(te, ['{0:.2f}'.format(f) for f in sorted(list(res_fcs[te]))])
        ps[te] = ttest_1samp(res_fcs[te], 0)[1] # For FC
        #print(ps[te], te)
    # Q value correct?

    print('{0:,} genes without a non-TE transcript '.format(has_no_nonte_transcript))
    print('{0:,} genes without a TE-containing transcript'.format(has_no_with_te_transcript))
    print('Found {0:,} genes with at least 1 non-TE transcript and 1 TE-containing transcript'.format(gene_with_noTE_and_TE_transcript))
    return res_fcs, ps

coding_bundles = shared_bundle.bundle_up_by_name('coding', all_genes, tes)
noncoding_bundles = shared_bundle.bundle_up_by_name('noncoding', all_genes, tes)
all_bundles = shared_bundle.bundle_up_by_name('all', all_genes, tes)

res_coding, p_coding = process_bundles(coding_bundles)
res_ncrna, p_ncrna  = process_bundles(noncoding_bundles)
res_all, p_all = process_bundles(all_bundles)

coding_tes = [
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

data = {te: res_coding[te] for te in coding_tes}
shared.boxplots('pc.pdf', data, qs=p_coding, no_TE_key=None)

noncoding_tes = data = [
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
    #'LTR:ERV1:HUERS-P3b',
    #'LTR:ERV1:MER110',

    'LTR:ERVK:HERVK',

    'LTR:ERVL-MaLR:MST-int',
    'LTR:ERVL-MaLR:THE1-int',
    #'', res_type['']['pc-all'],

    ]

data = {te: res_ncrna[te] for te in noncoding_tes}
shared.boxplots('ncrna.pdf', data, qs=p_ncrna, no_TE_key=None)

all_tes = list(set(coding_tes + noncoding_tes))
all_tes.sort()
data = {te: res_all[te] for te in all_tes}
shared.boxplots('all.pdf', data, qs=p_all, no_TE_key=None)
