
import sys, os
from glbase3 import *

sys.path.append('../../')
import shared

gencode = glload('../../gencode/hg38_gencode_v32.glb')
gencode_map = {gene['enst']:gene for gene in gencode}
cds = glload('coding_genes_with_local_CDS-predicted.glb')
annotations = glload('../packed/all_genes.glb')
cds = cds.map(genelist=annotations, key='transcript_id')

res = {'perfect': 0,
    'start_correct': 0,
    'end_correct': 0,
    'incorrect': 0,
    'not-tested': 0,
    'not-found': 0}

for idx, gene in enumerate(cds):
    #print(gene)
    if '=' not in gene['name']: # only keep genocde
        continue

    # get the GENCODE transcript
    enst = gene['enst'].split('.')[0]
    if enst not in gencode_map:
        print('Warning: {0} not found in GENCODE v32'.format(enst))
        res['not-found'] += 1
        continue
    gencode_gene = gencode_map[enst]

    #print(gene)
    #print(gencode_gene)

    gene_tlength = shared.get_transcript_length(gene)
    _, gencode_tlength, actual_cdsl, actual_cdsr, splice_sites = shared.convert_genocode_to_local(gencode_gene)

    if len(gencode_gene['cds_loc']) == 0: # gencode doesn't know the CDS
        res['not-tested'] += 1
        continue

    predicted_cdsl = gene['cds_local_locs'][0]
    predicted_cdsr = gene['cds_local_locs'][1]
    #print(gene['enst'], gene['loc'], gencode_gene['loc'], gencode_gene['cds_loc'])
    #print(predicted_cdsl, actual_cdsl, ':', predicted_cdsr, actual_cdsr, gene['strand'], ':', gencode_tlength, gene_tlength)
    #print((predicted_cdsr-predicted_cdsl)/3.0, (actual_cdsr-actual_cdsl)/3.0)

    if gencode_tlength != gene_tlength:
        # the assembly is truncated/expanded in some way relative to the GENCODE
        # Isaac allows ~25% loss of overlap at the 5' and 3' ends
        # I can still go ahead if the TSS are the same:
        if gene['strand'] == '+':
            tss_gene = gene['loc']['left']
            tss_gencode = gencode_gene['loc']['left']
        elif gene['strand'] == '-':
            tss_gene = gene['loc']['right']
            tss_gencode = gencode_gene['loc']['right']
        if tss_gene != tss_gencode:
            res['not-tested'] += 1
            continue

    predicted_cdsl = gene['cds_local_locs'][0]
    predicted_cdsr = gene['cds_local_locs'][1]

    # The positions are often very slightly wrong, due to a problem I can't chase down
    # It's okay in the domain plots as the differences are always tiny,
    # but here we need to be a little fuzzy on the ends particularly.
    if abs(predicted_cdsl-actual_cdsl) < 3 and abs(predicted_cdsr-actual_cdsr) < 3:
        res['perfect'] += 1
    elif abs(predicted_cdsl-actual_cdsl) < 3:
        res['start_correct'] += 1
    elif abs(predicted_cdsr-actual_cdsr) < 3:
        res['end_correct'] += 1
    else:
        res['incorrect'] += 1

for k in res:
    print('{0}:\t {1}\t{2:.1f}%'.format(k, res[k], res[k]/sum(res.values())*100))
