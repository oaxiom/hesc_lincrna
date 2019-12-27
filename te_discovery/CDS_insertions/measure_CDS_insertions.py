import sys, os
from glbase3 import *

def qcollide(al, ar, bl, br):
    return ar >= bl and al <= br

def contained(al, ar, bl, br): # Is A in B?
    return al >= bl and ar <= br

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

res = {'inframe_insertion': [], # Class 1 # numbers are per-transcript;
    'insertion_truncation': [], # Class 2
    'new_ATG': [],
    'disrupt_coding': [], # Class 3
    'insertion_alternate_cds': [], # Class 4
    'no_disruption_5prime': [], # Class 5
    'no_disruption_3prime': [], # Class 6
    'class_not_found': []}

total = 0
no_variants = 0

for idx, gene_name in enumerate(bundles):
    total += len(bundles[gene_name])

    all_types = [i['tags'][-1] for i in bundles[gene_name]]

    if '~' not in all_types: # No variants here
        no_variants += len(bundles[gene_name])
        continue

    if len(bundles[gene_name]) == 1 and all_types[0] == '=':
        # TODO: Check if it's a ~ and get the GENCODE canonical one;
        continue

    if '=' not in all_types:
        # TODO: Ugh, don't have the canonical transcript, need to add it from GENCODE
        continue

    # Add the doms key to all transcripts;
    for transcript in bundles[gene_name]:
        if transcript['transcript_id'] in tes:
            te = tes[transcript['transcript_id']]
            transcript['doms'] = te['doms']
        else:
            transcript['doms'] = []
    # divide into known and novel:
    known = [t for t in bundles[gene_name] if '=' in t['tags']]
    novel = [t for t in bundles[gene_name] if '~' in t['tags']]

    #print(known, novel)

    for transcript in novel:
        te = None
        if transcript['transcript_id'] in tes:
            te = tes[transcript['transcript_id']]
        if te:
            if transcript['coding'] == 'noncoding':
                res['disrupt_coding'].append(transcript)
                continue

            # find out if a TE overlaps the CDS:
            for t in te['doms']:
                if qcollide(t['span'][0], t['span'][1], transcript['cds_local_locs'][0], transcript['cds_local_locs'][1]):
                    # See if the TE is entirely contained:
                    if contained(t['span'][0], t['span'][1], transcript['cds_local_locs'][0], transcript['cds_local_locs'][1]):
                        te_length = (t['span'][1] - t['span'][0])
                        expected_cds_length_if_in_frame = (transcript['cds_local_locs'][1] - transcript['cds_local_locs'][0]) + te_length -1
                        cds_lengths = [i['cds_local_locs'][1]-i['cds_local_locs'][0]+te_length for i in known if i['coding'] == 'coding']
                        if expected_cds_length_if_in_frame in cds_lengths:
                            res['inframe_insertion'].append(transcript) # what about multiple TE insertions?
                            break
                    else: # It is flapping over the edge of the CDS;
                        te_length = (t['span'][1] - t['span'][0])
                        # See if it's at the end:
                        if t['span'][1] > transcript['cds_local_locs'][1]: # It's stopping the CDS
                            res['insertion_truncation'].append(transcript)
                            break
                        elif t['span'][0] < transcript['cds_local_locs'][0]: # It's at the START;
                            res['new_ATG'].append(transcript)
                            break
                        1/0 # should not be possible to get here;
                else: # No collision; check it's 5' or 3':
                    # Check that it still contains a CDS of the correct length:
                    expected_cds_length_if_in_frame = ((transcript['cds_local_locs'][1] - transcript['cds_local_locs'][0])-1)
                    cds_lengths = [i['cds_local_locs'][1]-i['cds_local_locs'][0] for i in known if i['coding'] == 'coding']
                    #print(expected_cds_length_if_in_frame, cds_lengths)
                    if expected_cds_length_if_in_frame in cds_lengths: # It's a simple insertion 5' or 3':
                        # I know it's not a collision, so just test the edge:
                        if t['span'][1] < transcript['cds_local_locs'][0]: # 5'
                            res['no_disruption_5prime'].append(transcript)
                            break
                        elif t['span'][0] > transcript['cds_local_locs'][1]:
                            res['no_disruption_3prime'].append(transcript)
                            break
                    else:
                        res['insertion_alternate_cds'].append(transcript)
                        break

                # If you find something you break, this should be 0 if it's all working;
                res['class_not_found'].append(transcript)

    #if idx > 5000:
    #    break

for k in res:
    gl = genelist()
    if res[k]:
        gl.load_list(res[k])
        gl.saveTSV('table_{0}.tsv'.format(k))
        gl.save('table_{0}.glb'.format(k))

print()
for k in res:
    print(k, len(res[k]))
print('No variants', no_variants)
print('Total', total)
