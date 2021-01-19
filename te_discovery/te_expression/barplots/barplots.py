import sys, os, numpy
from scipy.stats import mannwhitneyu, wilcoxon
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared
sys.path.append('../')
import shared_bundle

# These have the official GENCODE CDS, and the predicted (about ~80% accurate)
all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
tes = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}

coding_bundles = shared_bundle.bundle_up_by_name('coding', all_genes, tes)
noncoding_bundles = shared_bundle.bundle_up_by_name('noncoding', all_genes, tes)
all_bundles = shared_bundle.bundle_up_by_name('all', all_genes, tes)

# Barplots for these genes:
def mini_barplot(bundle, filename, gene_name, TE):
    all_TEs = bundle[gene_name]

    vals = [i['TPM'] for i in all_TEs]
    cols = []
    labs = ['{} {} {}'.format(i['name'], i['enst'], i['transcript_id']) for i in all_TEs]
    for i in all_TEs:
        if TE in [te['dom'] for te in i['doms']]:
            cols.append('red') # This TE
        elif i['doms']:
            cols.append('orange') # A differnet TE
        else:
            cols.append('blue') # No TE

    ys = numpy.arange(len(vals))

    mmheat_hei = 0.1+(0.04*len(vals))
    fig = plot.figure(figsize=[4,2])
    fig.subplots_adjust(left=0.7, top=mmheat_hei, right=0.95, bottom=0.1)
    ax = fig.add_subplot(111)
    ax.barh(ys, vals, 0.5, color=cols)
    ax.set_yticks(ys)
    ax.set_yticklabels(labs)
    [t.set_fontsize(6) for t in ax.get_yticklabels()]
    [t.set_fontsize(6) for t in ax.get_xticklabels()]
    fig.savefig(filename)
    plot.close(fig)

# Noncoding:
genes = {
    'LTR7': ['HDAC2-AS2', 'AC027288.1', 'FIRRE',],
    'HERVH': ['LINC01108', 'AL590705.1', 'AP002856.2', ],
    'L1M2': ['AC012213.1', ],
    'FRAM': ['FIRRE',],
    'HERVK': ['LINC02018', 'PDCL3P4',
        'AC108519.1',
        'AF228730.5',
        'AC068587.4'],

    'HERVK11': ['LINC00467', 'LINC00467',
        'LANCL1-AS1',
        'AC116562.4',
        'AF228730.5',
        'AC068587.4',],

    'HERVK14C': ['AC116562.4', 'AC068587.4', 'ZNF433-AS1', 'ZNF577'],
    'HERVK3': ['ZNF577', 'ZNF702P'],
    'HERVK9': ['DAP3', 'SDHA', 'SEPTIN7P2', 'ZNF433-AS1'],
    'HERVKC4': ['DAP3', 'ZNF577'],

    }

print(noncoding_bundles)

for te in genes:
    for g in genes[te]:
        mini_barplot(noncoding_bundles, 'ncrna-{0}-{1}.pdf'.format(te, g), g, te)

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
        mini_barplot(coding_bundles, 'pc-{0}-{1}.pdf'.format(te, g), g, te)

