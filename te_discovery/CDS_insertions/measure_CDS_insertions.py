import sys, os
from glbase3 import *

# These have the official GENCODE CDS, and the predicted (about ~80% accurate)
#canonical = glload('../../gencode/hg38_gencode_v32.pc.glb')
#gencode_cds = glload('../../transcript_assembly/get_CDS/gencode_cds.glb')
#canonical = canonical.map(genelist=gencode_cds, key='enst')
all_genes = glload('../../transcript_assembly/packed/all_genes.glb')
cds = glload('../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
cds = {gene['transcript_id']: gene for gene in cds}
tes = glload('../te_transcripts/transcript_table_merged.mapped.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}

# First I need to bundle them up by their name;
bundles = {}
for gene in all_genes:
    symbol = gene['name'].split(' ')[0]

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

print('Found {0:,} bundles of genes'.format(len(bundles)))

# Okay, now we check each bundle has at least 1 CDS:
newbundles = {}
for b in bundles:
    for gene in bundles[b]:
        if gene['coding'] == 'coding': # Make sure there's at least 1 coding in there
            newbundles[b] = bundles[b]
            break
bundles = newbundles

print('Found {0:,} bundles of genes with at least 1 coding gene'.format(len(bundles)))

res = {'insertion_no_disruption': 0, # numbers are per-transcript;
    'insertion_truncation': 0,
    'disrupt_coding': 0,
    'insertion_alternate_cds': 0,
    'no_disruption_5prime': 0,
    'no_disruption_3prime': 0,
    'no_variants': 0,
    'total': 0}

for idx, gene_name in enumerate(bundles):
    res['total'] += len(bundles[gene_name])
    if len(bundles[gene_name]) == 1:
        # TODO: Check if it's a ~ and get the GENCODE canonical one;
        continue

    all_types = [i['tags'][-1] for i in bundles[gene_name]]

    if '~' not in all_types: # No variants here
        res['no_variants'] += len(all_types)
        continue

    for transcript in bundles[gene_name]:
        te = None
        if transcript['transcript_id'] in tes:
            te = tes[transcript['transcript_id']]

        if te:
            if transcript['coding'] == 'noncoding':
                res['disrupt_coding'] += 1

            # find out if the TE is inside the CDS:
            for t in te:


    if idx > 2000:
        break

print()
for k in res:
    print(k, res[k])

