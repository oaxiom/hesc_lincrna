import sys, os, numpy, math, glob
from scipy.stats import mannwhitneyu, wilcoxon, ttest_ind, ttest_1samp
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared
sys.path.append('../')
import shared_bundle

[os.remove(f) for f in glob.glob('type/*.pdf')]

def get_num_te_in_utr(dataset, TE=None):
    data = {'utr5': [], 'cds': [], 'utr3': []}

    for n, gene in enumerate(dataset):
        pos = gene['cds_local_locs'] # return 0, tlength, cdsl, cdsr, splice_sites

        if pos[0] == pos[1]:
            # Bad CDS, skip this one;
            continue

        utr5_l = 0
        utr5_r = pos[0]

        cds_l = pos[0]
        cds_r = pos[1]
        cds_len  = pos[1] - pos[0]

        utr3_l = pos[1]
        utr3_r = gene['tlength']
        utr3_len = (utr3_r - utr3_l) # in case some utr = 0
        #print(utr5_l, utr5_r, pos, utr3_l, utr3_r, utr3_len)

        add_utr5 = None
        add_cds = None
        add_utr3 = None

        for d in gene['doms']:
            if TE:
                if d['dom'] not in TE:
                    continue

            s = d['span'][0]
            e = d['span'][1]

            if s <= utr5_r: # Inside UTR
                add_utr5 = math.log2(gene['TPM'])
            if e >= cds_l and s <= cds_r: # Inside CDS
                add_cds = math.log2(gene['TPM'])
            if utr3_len > 1 and e > utr3_l: # there are a bunch of messages with UTR3' = 1
                add_utr3 = math.log2(gene['TPM'])

        # Only add the TPM once per transcript;
        if add_utr5: # Inside UTR
            data['utr5'].append(add_utr5)
        if add_cds: # Inside CDS
            data['cds'].append(add_cds)
        if add_utr3: # there are a bunch of messages with UTR3' = 1
            data['utr3'].append(add_utr3)

    return data

all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
contains_te = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
contains_not_te = contains_te.map(genelist=all_genes, key='transcript_id', logic='notright')

cds = glload('../../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
cds = {i['transcript_id']: i for i in cds}

def correct_CDS(data):
    all_doms_to_do = []
    for gene in data: # Correct the CDS
        if gene['coding'] == 'noncoding':
            continue
        if '!' in gene['tags']: # I am not considering these, as they look dubious;
            continue

        if gene['transcript_id'] not in cds:
            print('Warning {0} not found'.format(gene['transcript_id']))
            continue

        gene['cds_local_locs'] = cds[gene['transcript_id']]['cds_local_locs']
        gene['tlength'] = cds[gene['transcript_id']]['tlength']

        if gene['cds_local_locs'][0] == gene['cds_local_locs'][1]: # I am unsure about the cds_loc;
            continue

        all_doms_to_do.append(gene)
    return all_doms_to_do

contains_te = correct_CDS(contains_te)
contains_not_te = correct_CDS(contains_not_te)

types = {
    'DNA': [i['name'] for i in dfam if i['type'] == 'DNA'],
    'LINE': [i['name'] for i in dfam if i['type'] == 'LINE'],
    'SINE': [i['name'] for i in dfam if i['type'] == 'SINE'],
    'LTR': [i['name'] for i in dfam if i['type'] == 'LTR'],
    'Retroposon': [i['name'] for i in dfam if i['type'] == 'Retroposon'],
    }

for T in types:
    te_tpms = get_num_te_in_utr(contains_te, TE=types[T])

    if max([len(te_tpms['utr5']), len(te_tpms['cds']), len(te_tpms['utr3'])]) < 20: # Ignore uncommon/empty
        continue

    res_pc = {
        "3' UTR": te_tpms['utr3'],
        "CDS": te_tpms['cds'],
        "5' UTR": te_tpms['utr5'],
        'no TE': [math.log2(i['TPM']) for i in contains_not_te],
        }

    q_pc = {}
    for k in res_pc:
        q_pc[k] = ttest_ind(res_pc['no TE'], res_pc[k], equal_var=False)[1]

    shared.boxplots('type/num_tes_pc-{0}.pdf'.format(T), res_pc, qs=q_pc, trim_low_samples=10)

