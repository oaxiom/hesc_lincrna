import sys, os
from glbase3 import *
import matplotlib.pyplot as plot
os.path.append('../../')
import shared

# These have the official GENCODE CDS, and the predicted (about ~80% accurate)
all_genes = glload('../../transcript_assembly/packed/all_genes.glb')
cds = glload('../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-predicted.glb')
cds = {gene['transcript_id']: gene for gene in cds}
tes = glload('../te_transcripts/transcript_table_gencode_all.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}

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

    # Get hte CDS info:
    if gene['transcript_id'] in cds:
        gene['cds_info'] = True # and presumably 'coding' == True?
        gene['cds_local_locs'] = cds[gene['transcript_id']]['cds_local_locs']
    else:
        gene['cds_info'] = False
        gene['cds_local_locs'] = (-1,-1)

    bundles[symbol].append(gene)

print('Found {0:,} genes'.format(len(bundles)))
bundles = {b: bundles[b] for b in bundles if len(bundles[b]) > 1}
print('Found {0:,} genes with >1 transcript'.format(len(bundles)))
transcript_variants_per_gene = [len(bundles[gene]) for gene in bundles]
# limit to 10+
transcript_variants_per_gene = [min(b, 11) for b in transcript_variants_per_gene]

# histogram;
fig = plot.figure()
ax = fig.add_subplot(111)
n, bins, patches = ax.hist(transcript_variants_per_gene, max(transcript_variants_per_gene), facecolor='blue', alpha=0.5)
ax.set_xlim([0, 11])
fig.savefig('trasncripts_per_gene.png')


# Okay, now we check each bundle has at least 1 CDS:
newbundles = {}
for b in bundles:
    has_coding = False
    has_noncoding = False
    for gene in bundles[b]:
        if gene['coding'] == 'coding':
            has_coding = True
        if gene['coding'] == 'noncoding':
            has_noncoding = True

        if has_coding and has_noncoding:
            newbundles[b] = bundles[b]
            break
bundles = newbundles

print('Found {0:,} bundles of genes with at least 1 coding variant, and at least 1 non-coding variant'.format(len(bundles)))

