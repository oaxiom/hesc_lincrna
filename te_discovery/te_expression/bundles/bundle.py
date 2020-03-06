import sys, os
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared

# These have the official GENCODE CDS, and the predicted (about ~80% accurate)
all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
tes = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}

def bundle_up_by_name(mode):
    # First I need to bundle them up by their name;
    bundles = {}
    for gene in all_genes:
        symbol = gene['name'].split(' ')[0].strip()

        if symbol not in bundles:
            bundles[symbol] = []

        if gene['enst'] in tes:
            gene['doms'] = tes[gene['enst']]['doms']
            gene['TEs'] = True
        else:
            gene['TEs'] = False
            gene['doms'] = []

        bundles[symbol].append(gene)

    # remove the genes that are only non-coding
    newbundles = {}
    for gene in bundles:
        for transcript in bundles[gene]:
            if mode == 'all' or transcript['coding'] == mode:
                newbundles[gene] = bundles[gene]
                break
    bundles = newbundles

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

coding_bundles = bundle_up_by_name('coding')
noncoding_bundles = bundle_up_by_name('noncoding')
noncoding_bundles = bundle_up_by_name('all')

print(coding_bundles)

# Broad summary:
for g in bundles:
