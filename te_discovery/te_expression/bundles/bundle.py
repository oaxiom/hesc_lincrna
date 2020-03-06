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

        if gene['transcript_id'] in tes:
            gene['doms'] = tes[gene['transcript_id']]['doms']
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

# Broad summary:
def process_bundles(bundle):
    res_fcs = {}
    gene_with_noTE_and_TE_transcript = 0
    has_no_with_te_transcript = 0
    has_no_nonte_transcript = 0

    for gene in bundle:
        tpms_for_no_te = []
        tpms_for_with_te = {}
        for transcript in bundle[gene]:
            if transcript['TEs']:
                for te in transcript['doms']:
                    if te['dom'] not in tpms_for_with_te:
                        tpms_for_with_te[te['dom']] = []
                    tpms_for_with_te[te['dom']].append(transcript['TPM'])
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
            fc = utils.fold_change(max(tpms_for_no_te), max(tpms_for_with_te[te]), pad=0.01)
            if te not in res_fcs:
                res_fcs[te] = []
            res_fcs[te].append(fc)
    print('{0:,} genes without a non-TE transcript '.format(has_no_nonte_transcript))
    print('{0:,} genes without a TE-containing transcript'.format(has_no_with_te_transcript))
    print('Found {0:,} genes with at least 1 non-TE transcript and 1 TE-containing transcript'.format(gene_with_noTE_and_TE_transcript))
    return res_fcs

coding_bundles = bundle_up_by_name('coding')
#noncoding_bundles = bundle_up_by_name('noncoding')
#noncoding_bundles = bundle_up_by_name('all')

res_coding = process_bundles(coding_bundles)


