import sys, os, numpy, math, copy
from scipy.stats import mannwhitneyu, wilcoxon, ttest_ind, ttest_1samp
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared
sys.path.append('../')
import shared_bundle

from statsmodels.stats.weightstats import ztest

# These have the official GENCODE CDS, and the predicted (about ~80% accurate)
all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
tes = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
dfam_dict = {te['name']: '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name']) for te in dfam}

dfam_family = {dfam_dict[te['name']]: '{0}:{1}'.format(te['type'], te['subtype']) for te in dfam}

# Broad summary:
def process_bundles(bundle):
    res_fcs = {}
    res_fcs_family = {}
    ps = {}
    ps_family = {}
    gene_with_noTE_and_TE_transcript = 0
    has_no_with_te_transcript = 0
    has_no_nonte_transcript = 0
    for_gl = []

    for gene in bundle:
        tpms_for_no_te = []
        tpms_for_with_te = {}
        temp_for_gl = []
        for transcript in bundle[gene]:
            if transcript['TEs']:
                unq_tes = set([t['dom'] for t in transcript['doms']])
                for te in unq_tes:
                    full_name = dfam_dict[te]
                    if full_name not in tpms_for_with_te:
                        tpms_for_with_te[full_name] = []
                    tpms_for_with_te[full_name].append(transcript['TPM'])
                    transcript['te'] = te
                    temp_for_gl.append(copy.deepcopy(transcript))
            else: # No TE:
                tpms_for_no_te.append(transcript['TPM'])
                transcript['te'] = '-'
                temp_for_gl.append(transcript)

        # Get FC:
        if not tpms_for_with_te:
            has_no_with_te_transcript += 1
            continue # No paired
        if not tpms_for_no_te: # There is a few!
            has_no_nonte_transcript += 1
            continue

        for_gl += temp_for_gl

        gene_with_noTE_and_TE_transcript += 1
        for te in tpms_for_with_te:
            # A few ways to do this, take the mean or the max
            fc = utils.fold_change(max(tpms_for_no_te), max(tpms_for_with_te[te]), pad=0.01) # correct way around
            if 'LTR7' in te and '7' in te[-1]:
                print(te, tpms_for_no_te, tpms_for_with_te[te], fc)
            #fc = utils.fold_change(numpy.mean(tpms_for_no_te), numpy.mean(tpms_for_with_te[te]), pad=0.01)

            if te not in res_fcs:
                res_fcs[te] = []
            res_fcs[te].append(fc)

            family = dfam_family[te]
            if family not in res_fcs_family:
                res_fcs_family[family] = []
            res_fcs_family[family].append(fc)

    # Figure out the P:
    ps = {}
    for te in res_fcs:
        t = ttest_1samp(res_fcs[te], 0)[1] # For FC
        ps[te] = t

    ps_family = {}
    for te in res_fcs_family:
        t = ttest_1samp(res_fcs_family[te], 0)[1] # For FC
        ps_family[te] = t

    gl = genelist()
    gl.load_list(for_gl)

    print('{0:,} genes without a non-TE transcript '.format(has_no_nonte_transcript))
    print('{0:,} genes without a TE-containing transcript'.format(has_no_with_te_transcript))
    print('Found {0:,} genes with at least 1 non-TE transcript and 1 TE-containing transcript'.format(gene_with_noTE_and_TE_transcript))
    return res_fcs, ps, gl, res_fcs_family, ps_family

coding_bundles = shared_bundle.bundle_up_by_name('coding', all_genes, tes)
noncoding_bundles = shared_bundle.bundle_up_by_name('noncoding', all_genes, tes)
all_bundles = shared_bundle.bundle_up_by_name('all', all_genes, tes)

res_coding, p_coding, gl_coding, res_coding_family, ps_coding_family = process_bundles(coding_bundles)
res_ncrna, p_ncrna, gl_ncrna, res_ncrna_family, ps_ncrna_family  = process_bundles(noncoding_bundles)
res_all, p_all, gl_all, res_all_family, ps_all_family = process_bundles(all_bundles)

gl_coding = gl_coding.getColumns(['ensg', 'enst', 'name', 'gene_symbol', 'tags', 'coding', 'expression', 'TPM', 'TEs', 'te'])
gl_coding.saveTSV('coding.tsv')
gl_ncrna = gl_ncrna.getColumns(['ensg', 'enst', 'name', 'gene_symbol', 'tags', 'coding', 'expression', 'TPM', 'TEs', 'te'])
gl_ncrna.saveTSV('ncrna.tsv')
gl_all = gl_all.getColumns(['ensg', 'enst', 'name', 'gene_symbol', 'tags', 'coding', 'expression', 'TPM', 'TEs', 'te'])
gl_all.saveTSV('all.tsv')

data = {}
for te in sorted(res_coding_family):
    if '.' in te:
        continue
    if 'tRNA' in te or 'Satellite' in te:
        continue
    data[te] = res_coding_family[te]

shared.boxplots('family-pc.pdf', data, qs=ps_coding_family, no_TE_key=None)

data = {}
for te in sorted(res_ncrna_family):
    if '.' in te:
        continue
    if 'tRNA' in te or 'Satellite' in te:
        continue
    data[te] = res_ncrna_family[te]

shared.boxplots('family-ncrna.pdf', data, qs=ps_ncrna_family, no_TE_key=None)

coding_tes = [
    'DNA:TcMar-Tigger:Tigger1',
    'DNA:hAT-Charlie:MER117',
    'SINE:MIR:MIR',
    'SINE:Alu:AluSp',
    'SINE:Alu:AluJb',
    'SINE:Alu:FRAM',
    'LINE:L2:L2',
    'LINE:L1:L1M2_orf2',
    'LINE:L1:L1M5_orf2',
    'LINE:L1:L1HS_5end',
    'LINE:L1:L1MB5_3end',
    'LTR:ERVL-MaLR:MLT1B',
    'LTR:ERV1:LTR7',
    'LTR:ERV1:LTR7Y',
    'LTR:ERV1:HERVH',
    'LTR:ERV1:HERV-Fc2',
    'LTR:ERV1:HERVE',

    'LTR:ERVK:HERVK',
    #'', res_type['']['pc-all'],
    ]

data = {te: res_coding[te] for te in coding_tes}
shared.boxplots('pc.pdf', data, qs=p_coding, no_TE_key=None)

noncoding_tes = [
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
    ]

data = {te: res_ncrna[te] for te in noncoding_tes}
ps = {te: p_ncrna[te] for te in noncoding_tes}
shared.boxplots('ncrna.pdf', data, qs=ps, no_TE_key=None)

all_tes = list(set(coding_tes + noncoding_tes))
all_tes.sort()
data = {te: res_all[te] for te in all_tes}
shared.boxplots('all.pdf', data, qs=p_all, no_TE_key=None)
