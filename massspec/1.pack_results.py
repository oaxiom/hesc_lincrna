
import glob
from glbase3 import *
import gzip as gzipfile

res_peps = {} # Peptides that pass
res_genes = [] # Genes that pass

delete_chars = set('0123456789.+-')
# So I can get hte locations;
all_fastas = genelist('2.blast_searches/all_masked_peptides.fa', format=format.fasta) # Just for getting the positions;
all_fastas = {i['name']: i['seq'] for i in all_fastas}
dfam = genelist('../te_discovery/dfam/dfam_annotation.tsv', format={'force_tsv': True, 'name': 0, 'type': 3, 'subtype': 4})
all_te_transcripts = glload('../te_discovery/te_transcripts/transcript_table_merged.mapped.glb')

CDSs = glload('../transcript_assembly/get_CDS/coding_genes_with_local_CDS-corrected.glb')
CDSs = {i['transcript_id']: i for i in CDSs}
all_te_transcripts = {i['transcript_id']: i['doms'] for i in all_te_transcripts}
#all_data = all_matches.map(genelist=all_te_transcripts, key='transcript_id')

for filename in glob.glob('3.hipsci_results/PT*.tsv.gz'):
    oh = gzipfile.open(filename, 'rt')
    print(filename)

    for idx, line in enumerate(oh):
        if line[0] == '#':
            continue

        line = line.strip().split('\t')

        #if 'table_noncoding_to_coding_noTE' in line[9]:
        #    # new coding, but I'm not interested in them here as they have no TE;
        #    continue
        #elif 'table_variant_coding_but_noTE' in line[9]:
        #    continue # same as above

        peptide = line[8]
        matches = line[9]

        # Only take high quality matches
        e = float(line[14]) # QValue column
        if e > 0.05: # FDR if you use q
            continue

        # Much more relaxed:
        #e = float(line[13]) # EValue column
        #if e > 0.001:
        #    continue

        # delete peptides matching to two peptides or more with the same name;
        matches = matches.split(');')
        hsc_names = [m.split('|')[2] for m in matches]
        enst_names = [m.split('|')[3].split('(')[0] for m in matches]
        symbol_names = [m.split('|')[1] for m in matches]
        class_names = [m.split('|')[0] for m in matches]

        if (idx+1) % 1000 == 0:
            print('{:,}'.format(idx))

        if peptide not in res_peps:
            res_peps[peptide] = 0
        res_peps[peptide] += 1

        for hsc, enst, symbol, class_ in zip(hsc_names, enst_names, symbol_names, class_names):
            peptide_string = ''.join([i for i in peptide if i not in delete_chars])
            #print(peptide_string, peptide)

            # Find out where it is in the CDS, and if it's in a TE:
            # Get the position in the Peptide_fasta
            fasta_name = '|'.join([class_, symbol.replace(' ', ''), hsc, enst])
            aa_seq = all_fastas[fasta_name]
            left = aa_seq.index(peptide_string)
            rite = left+len(peptide_string)

            # see if it's in a TE domain:
            fullname = 'No'
            if 'noTE' in line[9]: # No TE peptide;
                te_doms = []
            else:
                te_doms = all_te_transcripts[hsc]

            # the left/rite values don't include the UTRs. I need to add them as the doms are in mRNA +UTRs positions;
            mrna_left = CDSs[hsc]['cds_local_locs'][0] + (left*3)
            mrna_right = CDSs[hsc]['cds_local_locs'][0] + (rite*3)

            for d in te_doms:
                if d['span'][1] >= mrna_left and d['span'][0] <= mrna_right:
                    te = dfam.get(key='name', value=d['dom'])[0]
                    fullname = '{0}:{1}:{2}'.format(te['type'], te['subtype'], d['dom'])

            res_genes.append({'transcript_id': hsc,
                'enst': enst,
                'name': symbol,
                'class': class_,
                'peptide': peptide,
                'E': e,
                'peptide_string': peptide_string,
                'peptide_left': left,
                'peptide_right': rite,
                'mrna_left': mrna_left,
                'mrna_right': mrna_right,
                'insideTE': fullname,
                })

gl = genelist()
gl.load_list(res_genes)
gl.sort('name')
gl.sort('class')
gl.saveTSV('results_gene.tsv', key_order=['transcript_id', 'enst', 'name', 'class', 'E', 'peptide', 'peptide_string', 'insideTE'])
gl.save('results_gene.glb')

print('{0} peptides passed'.format(len(res_peps)))
print('{0} genes found'.format(len(set(gl['transcript_id']))))

