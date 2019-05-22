

'''

# Parse the hmmer TSVs and output the domain locations on a per-transcript basis;
Before [(280, 573), (280, 574), (280, 574), (280, 574), (280, 572), (280, 573), (280, 574), (280, 574), (280, 575), (281, 361), (417, 515), (547, 689), (671, 995), (671, 995), (671, 995), (672, 995), (673, 995), (673, 995), (673, 995), (673, 995), (673, 995), (913, 994), (1061, 1349), (1061, 1349), (1061, 1349), (1061, 1349), (1061, 1349), (1061, 1349), (1061, 1349), (1061, 1349), (1061, 1349), (1191, 1288), (1917, 2155), (1917, 2163), (1917, 2163), (1917, 2162), (1917, 2162), (1917, 2155), (1917, 2155), (1918, 2155), (1919, 2155), (1985, 2083), (2325, 2557), (2577, 2717), (2578, 2716), (2578, 2717), (2578, 2716), (2578, 2717), (2578, 2717), (2578, 2717), (2578, 2717), (2578, 2717), (2833, 2973), (2834, 2973), (2834, 2973), (2834, 2973), (2834, 2973), (2834, 2973), (2834, 2973), (2834, 2973), (2834, 2973), (3362, 3616), (3674, 3976), (3674, 3975), (3674, 3976), (3674, 3976), (3674, 3975), (3674, 3976), (3674, 3976), (3674, 3976), (3674, 3975), (3674, 3976), (3809, 3905), (3990, 4171)]

After [(280, 574), (281, 361), (547, 689), (671, 995), (913, 994), (1061, 1349), (1917, 2155), (2325, 2557), (2578, 2716), (2834, 2973), (3362, 3616), (3674, 3975), (3990, 4171)]

'''

import glob, sys, os, gzip
from operator import itemgetter
from glbase3 import utils, expression, genelist, glload

dfam = genelist('dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'HMMlength': 1, 'species': 2, 'type': 3, 'subtype': 4})

reject_species = set([
    'Chrysochloris_asiatica', # Golden mole;
    'Caenorhabditis_elegans',
    'Danio',
    'Danio_rerio',
    'Drosophila_fruit_fly_genus',
    'Ficedula_albicollis',
    'Halyomorpha_halys',
    'Muridae',
    'Murinae',
    'Mus_mouse_genus',
    'Mus_musculus',
    'Protostomia', # Humans are Deuterostomia
    'Rodentia',
    'Teleostei',
    'Testudines',
    ])

for f in sorted(glob.glob('raw_data/tblout*.tsv.gz')):
    skipped = 0
    transcript_table = {}
    sam = os.path.split(f)[1].replace('.tsv', '').replace('.gz', '').replace('domtblout.', '').replace('tblout.', '')
    print('Working %s' % sam)
    oh = gzip.open(f, 'rt')

    res_te_count = {}

    for n, line in enumerate(oh):
        if line[0] == '#':
            continue
        line = [l for l in line.strip().split(' ') if l]

        dom_name = line[2]

        species = dfam.get(key='name', value=dom_name)['species'][0]
        #print(species)

        if species in reject_species: # Only keep Hominid TEs
            skipped += 1
            continue

        if line[0] not in transcript_table: # Column is transcript ID
            transcript_table[line[0]] = {'transcript_id': line[0], 'doms': []}

        span = (min((int(line[6]), int(line[7]))), max((int(line[6]), int(line[7])))) # spans on the - strand go in the 'wrong' direction;

        transcript_table[line[0]]['doms'].append({'dom': dom_name, 'span': span, 'E': float(line[12]), 'strand': line[11]})
        #print(transcript_table[line[0]])

        if n % 1000000 == 0:
            print('Processed: {:,} domains'.format(n))
            #break
    oh.close()

    print('Skipped {:,} non-Hominid TEs'.format(skipped))

    # Save summary table
    # Can use a glbase, although saving it as a TSV would be a mess

    # Why does this part take so long?
    print('Removing duplicate domains in {:,} transcripts'.format(len(transcript_table)))
    # remove duplicate, overlapping doms; take the + strand, best scoring domains if there are any overlaps;
    for idx, transcript_id in enumerate(transcript_table):
        transcript = transcript_table[transcript_id]

        if len(transcript['doms']) <= 1:
            continue
        #all_spans = [d for d in transcript['doms']]
        marked_for_deletion = []
        for idx_dom1, dom1 in enumerate(transcript['doms']):
            if dom1 in marked_for_deletion: # if already marked for deletion; don't do it again, in case of ties
                continue
            for idx_dom2, dom2 in enumerate(transcript['doms']):
                if idx_dom1 >= idx_dom2: # do a triangular match, so that you only compare each pair once;
                    continue

                if dom2 in marked_for_deletion:# if already marked for deletion; don't do it again, in case of ties
                    continue

                dom1_span = dom1['span']
                dom2_span = dom2['span']
                # see if there is an overlap
                # Check for the simple case where two domains are identical:
                if dom1_span[0] == dom2_span[0] and dom1_span[1] == dom2_span[1]:
                    if dom1['E'] <= dom2['E']:
                        marked_for_deletion.append(dom2)
                        continue
                    else:
                        marked_for_deletion.append(dom1)
                        continue

                # The second simple case: One domain is entirely contained in another domain
                if (dom1_span[0] > dom2_span[0] and dom1_span[1] < dom2_span[1]): # dom1 is inside dom2:
                    if dom1['E'] <= dom2['E']:
                        marked_for_deletion.append(dom2)
                        continue
                    else:
                        marked_for_deletion.append(dom1)
                        continue

                if (dom2_span[0] > dom1_span[0] and dom2_span[1] < dom1_span[1]): # dom2 is inside dom1:
                    if dom1['E'] <= dom2['E']:
                        marked_for_deletion.append(dom2)
                        continue
                    else:
                        marked_for_deletion.append(dom1)
                        continue

                # The third complex case, check the percentage overlap:
                if (dom1_span[1] > dom2_span[0] and dom1_span[0] < dom2_span[1]): # Has an overlap;
                    # get number of base pairs overlapping:
                    bp_overlap = min(abs(dom1_span[1] - dom2_span[0]), abs(dom1_span[0] - dom2_span[1]))

                    perc_overlap = bp_overlap / ((dom1_span[1] - dom1_span[0]) + (dom2_span[1] - dom2_span[0]) - bp_overlap + 1)

                    #print(dom1['span'], dom2['span'], bp_overlap, perc_overlap)
                    if perc_overlap > 0.6:
                        # Whichever one has the lower E value:
                        if dom1['E'] <= dom2['E']:
                            marked_for_deletion.append(dom2)
                            continue
                        else:
                            marked_for_deletion.append(dom1)
                            continue

                        marked_for_deletion.append(dom2) # choose a random one;

        if marked_for_deletion: # Sometimes none;
            #print('Before', sorted([d['span'] for d in transcript['doms']], key=itemgetter(0)))
            marked_for_deletion = list(map(dict, frozenset(frozenset(i.items()) for i in marked_for_deletion)))
            for m in marked_for_deletion:
                transcript_table[transcript_id]['doms'].remove(m)

            #print()
            #print('After', sorted([d['span'] for d in transcript['doms']], key=itemgetter(0)))
            #print()
            #break

        if (idx+1) % 1000 == 0:
            print('Processed: {:,} transcripts'.format(idx+1))

    #print(transcript_table)

    gl = genelist()
    gl.load_list([transcript_table[i] for i in sorted(transcript_table)])
    gl.save('transcript_table_%s.glb' % sam)
    gl.saveTSV('transcript_table_%s.tsv' % sam) # A mess, just for checking
