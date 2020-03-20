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

    return bundles

coding_bundles = bundle_up_by_name('coding')
noncoding_bundles = bundle_up_by_name('noncoding')
#all_bundles = bundle_up_by_name('all')

# Barplots for these genes:
def barplot(bundle, filename, gene_name, TE):
    all_TEs = bundle[gene_name]

    vals = [i['TPM'] for i in all_TEs]
    cols = []
    labs = [i['transcript_id'] for i in all_TEs]
    for i in all_TEs:
        if TE in [te['dom'] for te in i['doms']]:
            cols.append('red') # This TE
        elif i['doms']:
            cols.append('orange') # A differnet TE
        else:
            cols.append('blue') # No TE

    ys = numpy.arange(len(vals))

    mmheat_hei = 0.1+(0.04*len(vals))
    fig = plot.figure(figsize=[2,2])
    fig.subplots_adjust(left=0.6, top=mmheat_hei, right=0.95, bottom=0.1)
    ax = fig.add_subplot(111)
    ax.barh(ys, vals, 0.5, color=cols)
    ax.set_yticks(ys)
    ax.set_yticklabels(labs)
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    fig.savefig(filename)

# Noncoding:
genes = {
    'LTR7': ['HDAC2-AS2'],
    'HERVK': ['PDCL3P4', 'AC108519.1'],
    }

for te in genes:
    for g in genes[te]:
        barplot(noncoding_bundles, 'ncrna-{0}-{1}.pdf'.format(te, g), g, te)

genes = {
    'HERVH': ['DISP1', 'CCDC141', 'TRMT44', 'KLKB1', 'PRKG1', 'GUCY2C', 'ABHD12B', 'PCDH11X'],

    'L1M5_orf2': [ 'CTPS2', 'ZNF542P', 'ZNF780B','ZNF69', ],

    'L1M2_orf2': [ 'CTPS1', 'UTY', 'ADGRL3', 'CCDC191',
        'R3HCC1L', 'ZNF91',  'LRRTM4', 'NCOA2', 'SHANK2', 'MED21', 'RAB21',
        ],

    'FRAM': ['CDH2', 'CNDP2', 'BACH2', 'FBXL4', 'GINS4', 'ZNF250', 'RBM14'],
    }

for te in genes:
    for g in genes[te]:
        barplot(coding_bundles, 'pc-{0}-{1}.pdf'.format(te, g), g, te)

