

'''

Work out the score for TEs

# Parse the hmmer TSVs and output summary statistics

'''

import glob, sys, os, gzip
from glbase3 import utils, expression, genelist, glload

all_res_te_count = {}

# load in the FASTA summary statistics:
res = {}
oh = open('lib_sizes.txt', 'r')
for lin in oh:
    if 'name' in lin:
        continue
    t = lin.strip().split('\t')
    res[t[0]] = {'num_fasta': int(t[1]), 'num_bps': int(t[2])}
oh.close()

for f in sorted(glob.glob('trans*.glb')):
    if 'pacbio' in f:
        continue

    sam = os.path.split(f)[1].replace('.glb', '').replace('.gz', '').replace('transcript_table_', '')

    glb = glload(f)

    print(sam)

    res_te_count = {}

    for n, item in enumerate(glb):
        for d in item['doms']:

            if d['dom'] not in res_te_count: # Column 3 in domtblout
                res_te_count[d['dom']] = 0 #
            res_te_count[d['dom']] += 1

        #if n > 10000:
            #break
    all_res_te_count[sam] = res_te_count

# Save summary table
all_te_keys = set(sum([list(all_res_te_count[sam].keys()) for sam in all_res_te_count], []))

all_key = list(all_res_te_count.keys())
header_row = ['name'] +\
    all_key +\
    ['per_entry_%s' % s for s in all_key] +\
    ['per_1e6bp_%s' % s for s in all_key] +\
    ['FC_%s' % s for s in all_key]

oh = open('summary_all.tsv', 'w')
oh.write('%s\n' % ('\t'.join(header_row)))

glct = []
glnorm = []
glperbp = []

for te in sorted(list(all_te_keys)):
    te_row = [te, ]
    expn_norm = {'name': te, 'conditions': []}
    expn_perbp = {'name': te, 'conditions': []}
    expn_counts = {'name': te, 'conditions': []}

    for sam in all_res_te_count:
        if te in all_res_te_count[sam]:
            te_row.append(all_res_te_count[sam][te])
            expn_counts['conditions'].append(all_res_te_count[sam][te])
        else:
            te_row.append(0)
            expn_counts['conditions'].append(0)

    # per_entry:
    for sam in all_res_te_count:
        if te in all_res_te_count[sam]:
            te_row.append(all_res_te_count[sam][te] / res[sam]['num_fasta'])
        else:
            te_row.append(0)
        expn_norm['conditions'].append(te_row[-1] *1e3) # i.e. per 1000k
    # per bp:
    for sam in all_res_te_count:
        if te in all_res_te_count[sam]:
            te_row.append((all_res_te_count[sam][te]  / res[sam]['num_bps'])*1e6)
        else:
            te_row.append(0)
        expn_perbp['conditions'].append(te_row[-1])

    glct.append(expn_counts)
    glnorm.append(expn_norm)
    glperbp.append(expn_perbp)

    oh.write('%s\n' % ('\t'.join([str(i) for i in te_row])))
oh.close()

# Get TE annotations;
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'HMMlength': 1, 'species': 2, 'type': 3, 'subtype': 4})

sam_order = [
    'gencode.v32.transcripts',
    'Homo_sapiens.GRCh38.cdna.all',
    'HSC_SR_PB_merged.transcripts',

    'gencode.v32.pc_transcripts',
    'Homo_sapiens.GRCh38.cds.all',
    'HSC_SR_PB_merged.pc',

    'gencode.v32.lncRNA_transcripts',
    'Homo_sapiens.GRCh38.ncrna',
    'NONCODEv5_human',
    'HSC_SR_PB_merged.ncrna',
    ]

new_names = [
    'gencode.all',
    'GRCh38.all',
    'custom.all',
    'gencode.pc',
    'GRCh38.pc',
    'custom.pc',
    'gencode.ncrna',
    'GRCh38.ncrna',
    'noncode.ncrna',
    'custom.ncrna',
    ]

expn = expression(loadable_list=glct, cond_names=all_res_te_count)
expn = dfam.map(genelist=expn, key='name')
expn = expn.sliceConditions(sam_order)
expn.setConditionNames(new_names)
print(expn)
expn.save('te_counts.glb')
expn.saveTSV('te_counts.tsv')

expn = expression(loadable_list=glnorm, cond_names=all_res_te_count)
expn = dfam.map(genelist=expn, key='name')
expn = expn.sliceConditions(sam_order)
expn.setConditionNames(new_names)
print(expn)
expn.save('te_norm_counts.glb')
expn.saveTSV('te_norm_counts.tsv')

expn = expression(loadable_list=glperbp, cond_names=all_res_te_count)
expn = dfam.map(genelist=expn, key='name')
expn = expn.sliceConditions(sam_order)
expn.setConditionNames(new_names)
print(expn)
expn.save('te_per1Mbp_seq.glb')
expn.saveTSV('te_per1Mbp_seq.tsv')

# Do the FCs here:
