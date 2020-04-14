import sys, os, numpy, math
from scipy.stats import mannwhitneyu, wilcoxon, ttest_ind, ttest_1samp
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared
sys.path.append('../')
import shared_bundle

all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
dfam = genelist('../../dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
contains_te = glload('../../te_transcripts/transcript_table_merged.mapped.glb')
contains_not_te = contains_te.map(genelist=all_genes, key='transcript_id', logic='notright')

print(contains_not_te)

res_pc = {10-k: [] for k in range(0, 11)}
res_ncrna = {10-k: [] for k in range(0, 11)}

for gene in contains_not_te:
    if gene['coding'] == 'coding':
        res_pc[0].append(math.log(gene['TPM'], 2))
    else:
        res_ncrna[0].append(math.log(gene['TPM'], 2))

for gene in contains_te:
    num_tes = min(len(gene['doms']), 10)
    if gene['coding'] == 'coding':
        res_pc[num_tes].append(math.log(gene['TPM'], 2))
    else:
        res_ncrna[num_tes].append(math.log(gene['TPM'], 2))

q_pc = []
for k in res_pc:
    q_pc.append(ttest_ind(res_pc[0], res_pc[k], equal_var=False)[1])

q_ncrna = []
for k in res_ncrna:
    q_ncrna.append(ttest_ind(res_ncrna[0], res_ncrna[k], equal_var=False)[1])
# The q's have to be inverted as the drawing is bottom to top;
q_pc.reverse()
q_ncrna.reverse()

print(res_pc)
shared.boxplots('num_tes_pc.pdf', res_pc, qs=q_pc, no_TE_key=None)

print(q_ncrna)
shared.boxplots('num_tes_ncrna.pdf', res_ncrna, qs=q_ncrna, no_TE_key=None)
