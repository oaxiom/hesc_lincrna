
import sys, os
from glbase3 import *

sys.path.append('../../')
import shared

gencode = glload('../../gencode/hg38_gencode_v32.glb')
gencode_map = {gene['enst']:gene for gene in gencode}
cds = glload('coding_genes_with_local_CDS.glb')
annotations = glload('../packed/all_genes.glb')
cds = cds.map(genelist=annotations, key='transcript_id')

res = {'perfect': 0,
    'start_correct': 0,
    'end_correct': 0,
    'incorrect': 0,
    'not-tested': 0}

for gene in cds:
    #print(gene)
    if '=' not in gene['name']: # only keep genocde
        continue

    # get the GENCODE transcript
    enst = gene['enst'].split('.')[0]
    if enst not in gencode_map:
        print('Warning: {0} not found in GENCODE v32'.format(enst))
        res['not-tested'] += 1
        continue
    gencode_gene = gencode_map[enst]

    print(gene)
    print(gencode_gene)

    gene_tlength = shared.get_transcript_length(gene)
    _, gencode_tlength, actual_cdsl, actual_cdsr, splice_sites = shared.convert_genocode_to_local(gencode_gene)

    if gencode_tlength != gene_tlength:
        # the assembly is truncated/expanded in some way relative to the GENCODE
        # Isaac allows ~25% loss of overlap at the 5' and 3' ends
        # TODO:
        print('Warning: lengths are not the same, positions may be incorrect; {0} {1}'.format(gencode_tlength, gene_tlength))
        res['not-tested'] += 1

        # I need to correct the

    else: # just deal with the naive case first:
        predicted_cdsl = gene['cds_local_locs'][0]
        predicted_cdsr = gene['cds_local_locs'][1]

        if predicted_cdsl == actual_cdsl and predicted_cdsr == actual_cdsr:
            res['perfect'] += 1
        elif predicted_cdsl == actual_cdsl:
            res['start_correct'] += 1
        elif predicted_cdsr == actual_cdsr:
            res['end_correct'] += 1
        else:
            res['incorrect'] += 1

        print(predicted_cdsl, actual_cdsl, ':', predicted_cdsr, actual_cdsr, gene['strand'], ':', gencode_tlength, gene_tlength)


    1/0
print(res)
