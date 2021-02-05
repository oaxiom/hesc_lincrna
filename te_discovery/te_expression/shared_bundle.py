import sys, os, numpy
from scipy.stats import mannwhitneyu, wilcoxon
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared

def bundle_up_by_name(mode, all_genes, tes, _draw_hist=True):
    # First I need to bundle them up by their name;
    bundles = {}
    for gene in all_genes:
        #if gene['expression'] == 'depleted':
        #    continue

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
    if _draw_hist:
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
            fc = utils.fold_change(max(tpms_for_no_te), max(tpms_for_with_te[te]), pad=0.01) # correct way around
            #print(te, max(tpms_for_no_te), max(tpms_for_with_te[te]), fc)
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
        ps[te] = mannwhitneyu(tpms_noTE[te], tpms_withTE[te], alternative='two-sided')[1]
    # Q value correct?

    print('{0:,} genes without a non-TE transcript '.format(has_no_nonte_transcript))
    print('{0:,} genes without a TE-containing transcript'.format(has_no_with_te_transcript))
    print('Found {0:,} genes with at least 1 non-TE transcript and 1 TE-containing transcript'.format(gene_with_noTE_and_TE_transcript))
    return res_fcs, ps
