import sys, os
from glbase3 import *
import matplotlib.pyplot as plot
sys.path.append('../../../')
import shared

# These have the official GENCODE CDS, and the predicted (about ~80% accurate)
all_genes = glload('../../../transcript_assembly/packed/all_genes.glb')
tes = glload('../../te_transcripts/transcript_table_merged.mapped.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}

# First I need to bundle them up by their name;
bundles = {}
for gene in all_genes:
    symbol = gene['name'].split(' ')[0].strip()
    
    if gene['expression'] != 'unbiased':
        continue
    
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
        if transcript['coding'] == 'coding':
            newbundles[gene] = bundles[gene]
            break
bundles = newbundles

all_coding_genes = len(bundles)
print('Found {0:,} genes'.format(all_coding_genes))
bundles = {b: bundles[b] for b in bundles if len(bundles[b]) > 1}
genes_with_multiple_transcripts = len(bundles)
print('Found {0:,} genes with >1 transcript'.format(genes_with_multiple_transcripts))
transcript_variants_per_gene = [len(bundles[gene]) for gene in bundles]
# limit to 10+
transcript_variants_per_gene = [min(b, 20) for b in transcript_variants_per_gene]
# histogram;
fig = plot.figure(figsize=[3,2])
ax = fig.add_subplot(111)
ax.hist(transcript_variants_per_gene, max(transcript_variants_per_gene)-1)
ax.set_xlim([-0.5, 21.5])
fig.savefig('transcripts_per_gene.png')
fig.savefig('transcripts_per_gene.pdf')

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

genes_one_coding_and_noncoding_variant = len(bundles)
print('Found {0:,} bundles of genes with at least 1 coding variant, and at least 1 non-coding variant'.format(genes_one_coding_and_noncoding_variant))

# Now, for each gene in each bundle, see if the non-coding contains a TE
hasTE = 0
noTE = 0
for g in bundles:
    for transcript in bundles[g]:
        if transcript['coding'] == 'noncoding':
            # see if it has a TE:
            if transcript['transcript_id'] in tes:
                hasTE += 1
                break
noTE = len(bundles) - hasTE

print('Has a TE non-coding variant: {0}'.format(hasTE))
print('Non-coding variant no TEs  : {0}'.format(noTE))

pies = {
    'Coding, single transcript': all_coding_genes - genes_with_multiple_transcripts,
    'Coding, no non-coding variant': genes_with_multiple_transcripts - genes_one_coding_and_noncoding_variant,
    'Coding, with a non-coding variant': genes_one_coding_and_noncoding_variant
    }

shared.pie('coding_to_non-coding.pdf', pies.values(), pies.keys())

pies = {
    'no TE': noTE,
    'Contains TE': hasTE,
    }

shared.pie('coding_to_non-coding_TE.pdf', pies.values(), pies.keys())