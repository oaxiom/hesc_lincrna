import sys, os, numpy
from scipy.stats import mannwhitneyu, wilcoxon
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared

# These have the official GENCODE CDS, and the predicted (about ~80% accurate)
all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
tes = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
dfam_dict = {te['name']: '{0}:{1}:{2}'.format(te['type'], te['subtype'], te['name']) for te in dfam}

def bundle_up_by_name(mode):
    # First I need to bundle them up by their name;
    bundles = {}
    for gene in all_genes:
        symbol = gene['name'].split(' ')[0].strip()

        if symbol not in bundles:
            bundles[symbol] = []

        if gene['transcript_id'] in tes:
            gene['doms'] = tes[gene['transcript_id']]['doms']
            gene['TEs'] = True
        else:
            gene['TEs'] = False
            gene['doms'] = []

        # remove the genes that are only coding/non-coding
        if mode != 'all' and gene['coding'] != mode:
           continue
        bundles[symbol].append(gene)

    print(mode)
    print('Found {0:,} genes'.format(len(bundles)))
    bundles = {b: bundles[b] for b in bundles if len(bundles[b]) > 1}
    genes_with_multiple_transcripts = len(bundles)
    print('Found {0:,} genes with >1 transcript'.format(genes_with_multiple_transcripts))

    transcript_variants_per_gene = [len(bundles[gene]) for gene in bundles]
    # limit to 10+
    transcript_variants_per_gene = [min(b, 20) for b in transcript_variants_per_gene]
    # histogram;
    fig = plot.figure(figsize=[1.6,1.1])
    ax = fig.add_subplot(111)
    ax.hist(transcript_variants_per_gene, max(transcript_variants_per_gene)-1, range=(0, 20))
    ax.set_xlim([-0.5, 21.5])

    ax.set_xticks([1.5, 10, 19.5])
    ax.set_xticklabels([2, 10, '>=20'])
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    fig.savefig('transcripts_per_gene-{0}.pdf'.format(mode))

    return bundles

# Broad summary:
def process_bundles(bundle, report_only):
    res_fcs = []

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
                    tpms_for_with_te[full_name].append(transcript)
            else: # No TE:
                tpms_for_no_te.append(transcript['TPM'])


        if not tpms_for_with_te:
            continue # No paired
        if not tpms_for_no_te: # There is a few!
            continue

        for te in tpms_for_with_te:
            if te not in report_only:
                continue
            fc = utils.fold_change(max(tpms_for_no_te), max([t['TPM'] for t in tpms_for_with_te[te]]), pad=0.01) # correct way around

            for transcript in tpms_for_with_te[te]:
                transcript['fc'] = fc
                transcript['te'] = te
                res_fcs.append(transcript)

    gl = genelist()
    gl.load_list(res_fcs)

    return gl

coding_bundles = bundle_up_by_name('coding')
noncoding_bundles = bundle_up_by_name('noncoding')
all_bundles = bundle_up_by_name('all')

coding_tes = [
    'DNA:TcMar-Tigger:Tigger1',
    'SINE:MIR:MIR',
    'SINE:Alu:AluSp',
    'LINE:L2:L2',
    'LINE:L1:L1M2_orf2',
    'LINE:L1:L1M5_orf2',
    'LINE:L1:L1HS_5end',
    'LTR:ERVL-MaLR:MLT1B',
    'LTR:ERV1:LTR7',
    'LTR:ERV1:LTR7Y',
    'LTR:ERV1:HERVH',
    #'', res_type['']['pc-all'],
    ]

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
    'LTR:ERV1:HUERS-P3b',
    'LTR:ERV1:MER110',

    'LTR:ERVK:HERVK',

    'LTR:ERVL-MaLR:MST-int',
    'LTR:ERVL-MaLR:THE1-int',
    #'', res_type['']['pc-all'],
    ]

res_coding = process_bundles(coding_bundles, coding_tes).saveTSV('coding.tsv')
res_ncrna  = process_bundles(noncoding_bundles, noncoding_tes).saveTSV('ncrna.tsv')
res_all  = process_bundles(all_bundles, list(set(coding_tes + noncoding_tes))).saveTSV('all.tsv')
