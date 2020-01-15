
import glob
from glbase3 import *
import gzip as gzipfile

res_peps = {} # Peptides that pass
res_genes = [] # Genes that pass

for filename in glob.glob('PT*.tsv.gz'):
    oh = gzipfile.open(filename, 'rt')
    print(filename)

    for idx, line in enumerate(oh):
        if line[0] == '#':
            continue

        line = line.strip().split('\t')

        peptide = line[8]
        matches = line[9]
        e = float(line[13])

        # Only take high quality matches
        if e > 0.1:
            continue

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
            res_genes.append({'transcript_id': hsc,
                'enst': enst,
                'name': symbol,
                'class': class_,
                'peptide': peptide,
                'e': e,
                })

gl = genelist()
gl.load_list(res_genes)
gl.sort('name')
gl.sort('class')
gl.saveTSV('results_gene.tsv', key_order=['transcript_id', 'enst', 'name'])
gl.save('results_gene.glb')

print('{0} peptides passed'.format(len(res_peps)))
print('{0} genes found'.format(len(set(gl['transcript_id']))))

