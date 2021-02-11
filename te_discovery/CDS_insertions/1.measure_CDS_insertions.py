import sys, os, glob
from glbase3 import *
sys.path.append('../../')
import shared

[os.remove(f) for f in glob.glob('*.tsv')]
[os.remove(f) for f in glob.glob('*.glb')]

def qcollide(al, ar, bl, br):
    return ar >= bl and al <= br # return(self.loc["right"] >= loc.loc["left"] and self.loc["left"] <= loc.loc["right"]) # nice one-liner

def contained(al, ar, bl, br): # Is A in B?
    return al > bl and ar < br

# These have the official GENCODE CDS, and the predicted (about ~80% accurate)
canonical = glload('../../gencode/hg38_gencode_v32.pc.glb')
gencode_cds = glload('../../transcript_assembly/get_CDS/gencode_cds.glb')
print(gencode_cds)
canonical_all = canonical.map(genelist=gencode_cds, key='enst')

canonical = {} # convert to a quick look up for speed
for gene in canonical_all:
    if gene['name'] not in canonical:
        canonical[gene['name']] = []
    gene['cds_info'] = True
    gene['cds_gencode_loc'] = gene['cds_loc']
    gene['cds_local_to_genome'] = gene['cds_loc']
    canonical[gene['name']].append(gene)

# get sequences for extra string literal matching;
gencode_peptide_fastas = genelist('../../transcript_assembly/get_CDS/gencode.v32.pc_translations.fa.gz', format=format.fasta, gzip=True)
#print(gencode_peptide_fastas)
gencode_peptide_fastas_lookup = {}
for gene in gencode_peptide_fastas:
    name = gene['name'].split('|')[6]
    if name not in gencode_peptide_fastas_lookup:
        gencode_peptide_fastas_lookup[name] = []
    gencode_peptide_fastas_lookup[name].append(gene['seq'])

all_transcripts = glload('../../transcript_assembly/packed/all_genes.glb')
cds = glload('../../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
cds = {gene['transcript_id']: gene for gene in cds}
#print(cds)
tes = glload('../te_transcripts/transcript_table_merged.mapped.glb') # yes, all as I am measuring PC -> ncRNA as well;
tes = {gene['transcript_id']: gene for gene in tes}

def check_sequence_versus_gencode(seq, gencode_peptide_fastas_lookup):
    if name in gencode_peptide_fastas_lookup:
        for seq in gencode_peptide_fastas_lookup[gene['name'].split(' ')[0]]:
            if seq == aa:
                print(name)
                print(seq)
                print(aa)
                print()
                found = True
                return True
    return False

# First I need to bundle them up by their name;
bundles = {}
for transcript in all_transcripts:
    symbol = transcript['name'].split(' ')[0]

    if symbol not in bundles:
        bundles[symbol] = []

    if transcript['transcript_id'] in tes:
        transcript['doms'] = tes[transcript['transcript_id']]['doms']
        transcript['TEs'] = True
    else:
        transcript['TEs'] = False
        transcript['doms'] = []

    # Get hte CDS info:
    if transcript['transcript_id'] in cds:
        transcript['cds_info'] = True # and presumably 'coding' == True?
        transcript['cds_local_locs'] = cds[transcript['transcript_id']]['cds_local_locs']
        transcript['cds_gencode_loc'] = cds[transcript['transcript_id']]['cds_gencode_loc']
        transcript['cds_local_to_genome'] = cds[transcript['transcript_id']]['cds_local_to_genome']
    else:
        # Probably non-coding;
        transcript['cds_info'] = False
        transcript['cds_local_locs'] = (-1,-1)
        transcript['cds_gencode_loc'] = None
        transcript['cds_local_to_genome'] = None

    bundles[symbol].append(transcript)

print('Found {:,} bundles of genes'.format(len(bundles)))

res = {'inframe_insertion': [], # Class 1 # numbers are per-transcript;
    'frameshift_insertion': [],
    'noncoding_to_coding_withTE': [],
    'noncoding_to_coding_noTE': [],
    'new_STOP': [], # Class 2
    'new_ATG': [],
    'coding_to_noncoding': [], # Class 3
    'insertion_alternate_cds': [], # Class 4
    'no_disruption_5prime': [], # Class 5
    'no_disruption_3prime': [], # Class 6
    'no_disruption_5_3prime': [],
    'class_not_found': [],
    'variant_coding_but_noTE': [],
    'no_coding': []}

total = 0
no_variants = 0
canonical_not_found = 0

for idx, gene_name in enumerate(bundles):
    if '-' in gene_name:
        continue # Skip these

    total += len(bundles[gene_name])

    all_types = [i['tags'][-1] for i in bundles[gene_name]]

    if '~' not in all_types: # No variants here
        no_variants += len(bundles[gene_name])
        continue

    if len(bundles[gene_name]) == 1 and all_types[0] == '=': # Only the canonical one was found, skip;
        continue

    # Always add the canonical transcripts from GENCODE, otherwise you just end up guessing existing CDS that have a slightly different transcript
    if gene_name in canonical:
        can = None
        can = canonical[gene_name]

        #can = canonical.getRowsByKey(key='name', values=gene_name, silent=True)
        if can:
            for i in can:
                i['tags'] = '='
                i['coding'] = 'coding'
                bundles[gene_name].append(i)
            #print(bundles[gene_name])
        else:
            #print(gene_name)
            canonical_not_found += 1 # probably non-coding
            continue
        # Don't worry about duplicate removal as we only care about the ~ transcripts anyways

    # check there is at least 1 coding in there, coming from either the GENCODE canonical, or internally
    if 'coding' not in [i['coding'] for i in bundles[gene_name]]:
        res['no_coding'] += bundles[gene_name]
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

    # The positions are often slightly off by -2, -1, 1 and 2 bp;
    cds_lengths = []
    cds_lengths = [(i['cds_local_locs'][1]-i['cds_local_locs'][0]) for i in known if i['coding'] == 'coding']
    cds_inframe = set(cds_lengths)
    cds_lengths += [l-2 for l in cds_lengths]
    cds_lengths += [l-1 for l in cds_lengths]
    cds_lengths += [l+1 for l in cds_lengths]
    cds_lengths += [l+2 for l in cds_lengths]
    cds_lengths = set(cds_lengths)

    canonical_cds_edges_genome = []
    for t in known:
        #print('\n',t)
        if t['coding'] != 'coding': continue
        if not t['cds_info']: continue
        if 'cds_gencode_loc' not in t: continue

        if transcript['strand'] == '+':
            canonical_cds_edges_genome = []
            if t['cds_gencode_loc']:
                edge = t['cds_gencode_loc']['left']
                canonical_cds_edges_genome += [edge, edge-1, edge+1, edge+2]
            edge = t['cds_local_to_genome']['left']
            canonical_cds_edges_genome += [edge, edge-1, edge+1, edge+2]
        else:
            canonical_cds_edges_genome = []
            if t['cds_gencode_loc']:
                edge = t['cds_gencode_loc']['right']
                canonical_cds_edges_genome += [edge, edge-1, edge+1, edge+2]
            edge = t['cds_local_to_genome']['right']
            canonical_cds_edges_genome += [edge, edge-1, edge+1, edge+2]

    canonical_cds_edges_genome = set(canonical_cds_edges_genome)

    for transcript in novel:
        # calls:
        inframe_insertion = False
        frameshift_insertion = False
        insertion_alternate_cds = False
        new_STOP = False
        new_ATG = False
        no_disruption_5prime = False
        no_disruption_3prime = False
        noncoding_to_coding_withTE = False
        noncoding_to_coding_noTE = False
        variant_coding_but_noTE = False

        # noncoding to coding
        known_coding_status = set([t['coding'] for t in bundles[gene_name] if '=' in t['tags']])
        novel_coding_status = transcript['coding']

        te = None
        if transcript['transcript_id'] in tes:
            te = tes[transcript['transcript_id']]

        if transcript['coding'] == 'noncoding':
            res['coding_to_noncoding'].append(transcript)
            continue

        if te:
            # find out if a TE overlaps the CDS:
            for t in te['doms']:
                # Collect the parameters:
                #print(t, '\n', transcript)
                colliding = qcollide(t['span'][0], t['span'][1], transcript['cds_local_locs'][0], transcript['cds_local_locs'][1])
                enclosed = contained(t['span'][0], t['span'][1], transcript['cds_local_locs'][0], transcript['cds_local_locs'][1])
                te_length = (t['span'][1] - t['span'][0])
                expected_cds_length = (transcript['cds_local_locs'][1] - transcript['cds_local_locs'][0])

                # cut the te for partially translated TEs:
                te_edges = (max(t['span'][0], i['cds_local_locs'][0]), min(t['span'][1], i['cds_local_locs'][1]))
                te_span = te_edges[1] - te_edges[0]
                cds_lengths_plus_te = [l+te_span for l in cds_inframe]
                cds_edges = [i['cds_local_locs'][0] for i in known if i['coding'] == 'coding']
                cds_edges += [i['cds_local_locs'][1] for i in known if i['coding'] == 'coding']

                if colliding: # with this CDS;
                    if 'coding' in novel_coding_status and 'coding' not in known_coding_status:
                        noncoding_to_coding_withTE = True # This will override all classes

                    if enclosed: # TE is entirely contained in the transcript
                        if expected_cds_length in cds_lengths_plus_te: # TE is most likely in frame inserted:
                            inframe_insertion = True
                        elif expected_cds_length in cds_lengths: # It's probably already annotated as contained, and has no effect on the CDS
                            frameshift_insertion = True# This will almost certainly get trimmed in the BLAST step so it's safe to leave it in this class;
                            #pass # I suppose it could get here by chance, but chances are less than 1 in 1e5 assuming ~1000 bp for both transcript and TE.
                        else: # probably a new CDS
                            frameshift_insertion = True

                    else: # It is colliding, but extends past the CDS, i.e. the STOP is inside the TE.
                        if t['span'][1] >= transcript['cds_local_locs'][1]: # It's STOP the CDS
                            if expected_cds_length in cds_lengths: # It's probably already annotated as contained, and has no effect on the CDS
                                pass
                            else: # probably a novel truncation
                                new_STOP = True

                        elif t['span'][0] <= transcript['cds_local_locs'][0]: # It's at the START;
                            if expected_cds_length in cds_lengths: # It's probably already annotated as contained, and has no effect on the CDS
                                pass
                            else: # probably a novel truncation
                                new_ATG = True

                        else: # I can't ID it;
                            1/0 # Never gets here;
                            frameshift_insertion = True

                else: # No collision with this CDS; check it's 5' or 3':
                    # Check it against the canonical CDSs;
                    if expected_cds_length in cds_lengths: # It's a simple insertion 5' or 3':
                        # I know it's not a collision, so just test the edge:
                        if t['span'][1] < transcript['cds_local_locs'][0]: # 5'
                            no_disruption_5prime = True
                        elif t['span'][0] > transcript['cds_local_locs'][1]:
                            no_disruption_3prime = True
                        else:
                            frameshift_insertion = True # It would just get trimmed in the BLAST step;

                    elif transcript['cds_local_to_genome']['left'] in canonical_cds_edges_genome and transcript['cds_local_to_genome']['right'] in canonical_cds_edges_genome:
                        no_disruption_5prime = True # No disruption to CDS
                        no_disruption_3prime = True

                    else:
                        if 'coding' in novel_coding_status and 'coding' not in known_coding_status:
                            noncoding_to_coding_noTE = True # No TE IN CDS! This will override all classes

                        #print('\n', transcript)
                        # see if one of the cds edges perfectly matches a canonical edge: Most likey a mid_CDS_insertion, that results in a STOP before the TE (hence no collision)'
                        if transcript['cds_local_locs'][0] in cds_edges or transcript['cds_local_locs'][1] in cds_edges:
                            insertion_alternate_cds = True

                        if transcript['cds_local_to_genome']['left'] in canonical_cds_edges_genome or transcript['cds_local_to_genome']['right'] in canonical_cds_edges_genome:
                            insertion_alternate_cds = True

                        else:
                            insertion_alternate_cds = True
                            # I find this category to be a bit dubious, and seems to have too many False+

        else: # No TE
            if 'coding' in novel_coding_status and 'coding' not in known_coding_status:
                noncoding_to_coding_noTE = True # This will override all classes
            else:
                variant_coding_but_noTE = True

        # transcripts only get called once. Add it here based ona hierarchy:
        # Use the calls above to assign to the preferred classes:
        if noncoding_to_coding_withTE:                       res['noncoding_to_coding_withTE'].append(transcript)
        elif noncoding_to_coding_noTE:                       res['noncoding_to_coding_noTE'].append(transcript)
        elif inframe_insertion:                              res['inframe_insertion'].append(transcript)
        elif new_ATG:                                        res['new_ATG'].append(transcript)
        elif new_STOP:                                       res['new_STOP'].append(transcript)
        elif frameshift_insertion:                           res['frameshift_insertion'].append(transcript)
        elif no_disruption_5prime and no_disruption_3prime:  res['no_disruption_5_3prime'].append(transcript)
        elif no_disruption_5prime:                           res['no_disruption_5prime'].append(transcript)
        elif no_disruption_3prime:                           res['no_disruption_3prime'].append(transcript)
        elif insertion_alternate_cds:                        res['insertion_alternate_cds'].append(transcript)
        elif variant_coding_but_noTE:                        res['variant_coding_but_noTE'].append(transcript)
        else:
            res['class_not_found'].append(transcript)
            #print(transcript)

    #if idx > 5000:
    #    break

for k in res:
    gl = genelist()
    if res[k]:
        gl.load_list(res[k])
        gl.saveTSV('table_{}.tsv'.format(k))
        gl.save('table_{}.glb'.format(k))

print()
for k in res:
    print(k, len(res[k]))
print('No variants', no_variants)
print('Total', total)
print('canonical_not_found', canonical_not_found)
